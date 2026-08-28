r"""
Extensive form games.

This module provides :class:`ExtensiveFormGame`, a thin wrapper around the
extensive form (game tree) ``Game`` of gambit ([Gambit]_).  An extensive form
game [NN2007]_ is described by a game tree; rather than re-implementing that
tree in Sage, this class simply holds a ``pygambit`` tree game and delegates
to it, in the same spirit as
:class:`~sage.game_theory.normal_form_game.NormalFormGame`.

A game can either be built in Sage with the tree-building methods
(:meth:`~ExtensiveFormGame.add_player`,
:meth:`~ExtensiveFormGame.append_move`, ...) or wrapped from an existing
``pygambit`` game.  Following gambit, tree nodes carry no labels: a node is
referred to by the ``pygambit`` ``Node`` object reached by navigating action
labels down from the root, ``g.root``.  For instance, in the game below
``battle.root`` is the root and ``battle.root.children['game']`` is the child
reached by Amy's action ``'game'``.  The moves themselves take the **action**
labels (a required argument).

Its other main purpose is to provide a Sage entry point to gambit's extensive
form games: convert back to the underlying gambit game with
:meth:`~ExtensiveFormGame._gambit_`, save and load games in gambit's ``.efg``
format with :meth:`~ExtensiveFormGame.save_efg` /
:meth:`~ExtensiveFormGame.load_efg`, load games from the literature out of the
gambit catalog with :meth:`~ExtensiveFormGame.load_from_gambit_catalog` (two of
them are worked through at the end of this page), and compute Nash equilibria
with :meth:`~ExtensiveFormGame.obtain_nash`.

Game trees can be drawn with :meth:`~ExtensiveFormGame.plot`, by default with
Sage's own graph plotting, which needs nothing beyond Sage, or -- for
publication-quality TikZ pictures like the two below -- with ``backend='gtdraw'``
and the optional `gtdraw
<https://www.gambit-project.org/gtdraw/>`_ package.  The latter gives a
:class:`GameTreeTikzPicture`; displaying one compiles it with LaTeX and opens
the picture, so seeing it needs a LaTeX installation on top of the package
itself (see the `gtdraw installation guide
<https://www.gambit-project.org/gtdraw/installation/>`_), while reading its
TikZ source does not.

EXAMPLES:

A game is built entirely in Sage, move by move, and the tree records not only
who moves when but also what each player knows when they move.  Take the Battle
of the Sexes: Amy and Bob want to spend the evening together but disagree on
how, Amy preferring video games and Bob a movie.  Suppose first that Amy
chooses and that Bob then chooses, having seen what she picked::

    sage: # optional - pygambit
    sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
    sage: battle = ExtensiveFormGame(players=['Amy', 'Bob'])
    sage: battle.append_move(battle.root, 'Amy', ['game', 'movie'])
    sage: for amy in ['game', 'movie']:
    ....:     battle.append_move(battle.root.children[amy], 'Bob',
    ....:                        ['game', 'movie'])
    sage: payoffs = {('game', 'game'): [3, 2], ('game', 'movie'): [1, 1],
    ....:            ('movie', 'game'): [0, 0], ('movie', 'movie'): [2, 3]}
    sage: for (amy, bob), payoff in payoffs.items():
    ....:     battle.set_outcome(battle.root.children[amy].children[bob],
    ....:                        '{},{}'.format(amy, bob), payoff)
    sage: battle
    An extensive form game with 2 players

Drawing the game with ``battle.plot(backend='gtdraw')`` shows the tree the
moves have built:
Amy's move at the root, one of Bob's below each of her actions, and the payoffs
to Amy and Bob at the four ends, hers written above his.

.. image:: ../../media/battle-of-the-sexes.png
   :align: center
   :width: 350 px
   :alt: Amy moves at the root and Bob at each of her two children.

This is a game of perfect information: each of its three information sets --
one for Amy's move and one for each of Bob's two moves -- holds a single node,
so Bob always knows where in the tree he is::

    sage: # optional - pygambit
    sage: len(battle.infosets)
    3
    sage: sorted(len(list(s.members)) for s in battle.infosets)
    [1, 1, 1]

Because Bob can react to what he saw, Amy can commit to video games and count
on him to follow, an equilibrium worth 3 to her::

    sage: # optional - pygambit
    sage: amy = battle._gambit_().players['Amy']
    sage: sorted(float(eq.payoff(amy))
    ....:        for eq in battle.obtain_nash(algorithm='enumpure'))
    [2.0, 3.0, 3.0]

Imperfect information -- a player who cannot tell certain nodes of the tree
apart -- is expressed by collecting those nodes into a single information set.
All nodes of an information set carry one and the same move: the same player,
choosing among the same actions.  That is what not being able to tell them
apart means, as a player who could choose differently at two nodes would
thereby distinguish them.  In Sage the first node of a set is given its move
with :meth:`~ExtensiveFormGame.append_move` as above, and every further node is
added to that set with :meth:`~ExtensiveFormGame.append_infoset`, naming a node
whose information set to join.

The Battle of the Sexes is really played simultaneously: neither player sees
the other's choice.  Simultaneity is drawn as Amy moving first with Bob's two
nodes joined into one information set, so that his choice cannot depend on
hers::

    sage: # optional - pygambit
    sage: battle = ExtensiveFormGame(players=['Amy', 'Bob'])
    sage: battle.append_move(battle.root, 'Amy', ['game', 'movie'])
    sage: battle.append_move(battle.root.children['game'], 'Bob',
    ....:                    ['game', 'movie'])
    sage: battle.append_infoset(battle.root.children['movie'],
    ....:                       battle.root.children['game'])
    sage: for (amy, bob), payoff in payoffs.items():
    ....:     battle.set_outcome(battle.root.children[amy].children[bob],
    ....:                        '{},{}'.format(amy, bob), payoff)

The tree is the one drawn above, except that Bob's two nodes are now enclosed
together in the shaded information set carrying his name, which is how a
drawing says that he cannot tell them apart:

.. image:: ../../media/battle-of-the-sexes-simultaneous.png
   :align: center
   :width: 350 px
   :alt: The same tree, with Bob's two nodes enclosed in one information set.

Bob's two nodes now share one information set of two members, so the game has
two information sets rather than three; the shared set carries his single
move::

    sage: # optional - pygambit
    sage: len(battle.infosets)
    2
    sage: bob = battle.root.children['game'].infoset
    sage: battle.root.children['movie'].infoset == bob
    True
    sage: len(list(bob.members))
    2
    sage: bob.player.label
    'Bob'
    sage: [action.label for action in bob.actions]
    ['game', 'movie']

Neither player ever forgets what they knew, so the game still has perfect
recall.  What has changed is that Bob can no longer follow Amy: the only pure
equilibria are the two ways of agreeing in advance, and Amy has no way of
forcing the one she prefers::

    sage: # optional - pygambit
    sage: battle.is_perfect_recall
    True
    sage: amy = battle._gambit_().players['Amy']
    sage: sorted(float(eq.payoff(amy))
    ....:        for eq in battle.obtain_nash(algorithm='enumpure'))
    [2.0, 3.0]

Any simultaneous two-player game is built in just this way: one move for the
first player, then one move and one :meth:`~ExtensiveFormGame.append_infoset`
per further node for the second::

    sage: # optional - pygambit
    sage: def simultaneous(players, actions, payoffs):
    ....:     g = ExtensiveFormGame(players=players)
    ....:     g.append_move(g.root, players[0], actions)
    ....:     first = g.root.children[actions[0]]
    ....:     g.append_move(first, players[1], actions)
    ....:     for action in actions[1:]:
    ....:         g.append_infoset(g.root.children[action], first)
    ....:     for (x, y), payoff in payoffs.items():
    ....:         g.set_outcome(g.root.children[x].children[y],
    ....:                       '{},{}'.format(x, y), payoff)
    ....:     return g

Matching pennies is the constant-sum game in which Alice and Bob each show a
coin, Alice taking both if they match and Bob if they differ.  It has no pure
equilibrium at all: knowing the other's choice would be decisive, and since
neither does, tossing the coin is all either can do::

    sage: # optional - pygambit
    sage: pennies = simultaneous(['Alice', 'Bob'], ['heads', 'tails'],
    ....:                        {('heads', 'heads'): [1, -1],
    ....:                         ('heads', 'tails'): [-1, 1],
    ....:                         ('tails', 'heads'): [-1, 1],
    ....:                         ('tails', 'tails'): [1, -1]})
    sage: pennies.obtain_nash(algorithm='enumpure')
    []
    sage: actions = pennies._gambit_().actions
    sage: [[float(eq[a]) for a in actions]
    ....:  for eq in pennies.obtain_nash(algorithm='lp')]
    [[0.5, 0.5, 0.5, 0.5]]

The games so far have been small enough to build by hand.  Gambit also ships a
catalog of games taken from the literature, which
:meth:`~ExtensiveFormGame.gambit_catalog_games` lists and
:meth:`~ExtensiveFormGame.load_from_gambit_catalog` loads into a Sage game.  The two below are the kind of game the extensive form exists for: what
makes each of them worth studying -- what a player knows, and which parts of
the tree are ever reached -- is precisely what disappears when the game is
flattened into a payoff matrix.

A driver leaving a party has to take the second exit off the motorway.
Exiting at the first is worth 0, continuing and then exiting 4, and driving
past both 1.  The two intersections look exactly alike, and the driver -- this
is Piccione and Rubinstein's absent-minded driver [PR1997]_, in the version
Gilboa put in the catalog -- cannot remember whether one has gone by already.
Both intersections therefore lie in a single information set.  That is not the
imperfect information of the games above, where Bob did not know somebody
else's move: here a player has forgotten their own, which is imperfect
*recall*::

    sage: # optional - pygambit
    sage: driver = ExtensiveFormGame.load_from_gambit_catalog(
    ....:     'journals/geb/gilboa1997/fig1')
    sage: driver
    An extensive form game with 1 player
    sage: driver.is_perfect_recall
    False
    sage: len(driver.infosets)
    1
    sage: sorted(len(list(s.members)) for s in driver.infosets)
    [2]
    sage: driver.root.children['B'].infoset == driver.root.infoset
    True

One move is made at both intersections, so a plan is a single probability `p`
of driving on (the action ``'B'``, against ``'E'`` for exiting), worth
`4p(1 - p) + p^2 = 4p - 3p^2`.  Gambit's equilibrium solvers assume perfect
recall, so rather than calling :meth:`~ExtensiveFormGame.obtain_nash` we hand
gambit the plans themselves and ask what each is worth::

    sage: # optional - pygambit
    sage: from fractions import Fraction
    sage: game = driver._gambit_()
    sage: def value(p):
    ....:     plan = game.mixed_behavior_profile(rational=True)
    ....:     plan[game.actions['B']] = Fraction(str(p))
    ....:     plan[game.actions['E']] = Fraction(str(1 - p))
    ....:     return QQ(plan.payoff(game.players['Player 1']))
    sage: [value(p) for p in [0, 1/3, 1/2, 2/3, 1]]
    [0, 1, 5/4, 4/3, 1]

Driving on with probability `2/3` is worth `4/3`, strictly more than either
way of deciding in advance, which are worth 0 and 1.  There is no opponent
here to keep guessing: tossing a coin buys the absent-minded driver something
that no deterministic plan can, and that is what forgetting does to a game.

The other game is Selten's horse, from the paper [Selten1975]_ that introduced
trembling-hand perfection.  Player 1 either passes across to Player 2 or drops
down to Player 3; Player 2, if reached, either ends the game at `(1, 1, 1)` or
passes down to Player 3 as well.  Player 3 cannot tell which of the two routes
led to them, and because that one information set straddles both branches, no
node below the root starts a subgame of its own::

    sage: # optional - pygambit
    sage: horse = ExtensiveFormGame.load_from_gambit_catalog(
    ....:     'journals/ijgt/selten1975/fig1')
    sage: horse
    An extensive form game with 3 players
    sage: [n.is_subgame_root for n in [horse.root, horse.root.children['R'],
    ....:                              horse.root.children['L']]]
    [True, False, False]

The game has two equilibria in pure strategies.  Each of them leaves one
player's information set off the equilibrium path -- it is reached with
probability 0 -- so what that player would have done is never put to the
test.  ``'enumpure'`` enumerates the pure strategies, so it answers with mixed
strategy profiles; ``as_behavior`` turns one into the behavior profile that
gives the tree its probabilities::

    sage: # optional - pygambit
    sage: eqs = horse.obtain_nash(algorithm='enumpure')
    sage: [[QQ(eq.payoff(p)) for p in horse.players] for eq in eqs]
    [[3, 2, 2], [1, 1, 1]]
    sage: game = horse._gambit_()
    sage: sets = [list(game.players[n].infosets)[0]
    ....:         for n in ['Player 1', 'Player 2', 'Player 3']]
    sage: [[QQ(eq.as_behavior().infoset_prob(s)) for s in sets] for eq in eqs]
    [[1, 0, 1], [1, 1, 0]]

Take the first one, worth `(3, 2, 2)`.  Player 1 goes down, Player 2 is never
asked to move, and the plan the equilibrium credits to Player 2 is to end the
game for 1 -- which is exactly what keeps Player 1 from passing across.  Give
Player 1's move to ``'R'`` so that Player 2's node is reached, leave the other
two players as they are, and ask what Player 2's two actions are then worth::

    sage: # optional - pygambit
    sage: plan = game.mixed_behavior_profile(rational=True)
    sage: for s in sets:
    ....:     plan[s.actions['L']] = Fraction('0')
    ....:     plan[s.actions['R']] = Fraction('1')
    sage: [(a.label, QQ(plan.action_value(a))) for a in sets[1].actions]
    [('R', 1), ('L', 4)]

Player 2 would pass down and take 4.  The profile is a Nash equilibrium all
the same, because being one only requires plans to be optimal where they are
actually carried out, and it survives subgame perfection too: with no proper
subgame there is nothing for that refinement to check.  Ruling it out is what
Selten's trembling hands are for, and the logit tracing procedure, which
approaches an equilibrium along a path of slightly trembling play, keeps only
the other one::

    sage: # optional - pygambit
    sage: [[QQ(eq.payoff(p)) for p in horse.players]
    ....:  for eq in horse.obtain_nash(algorithm='logit')]
    [[1, 1, 1]]

REFERENCES:

- [NN2007]_

- [PR1997]_

- [Selten1975]_

- [Gambit]_

AUTHORS:

- Amelie Kleber (06-2026): Initial version
"""

# ****************************************************************************
#       Copyright (C) 2026 Amelie Kleber
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

import os
import warnings
from tempfile import TemporaryDirectory

from sage.structure.sage_object import SageObject
from sage.misc.temporary_file import atomic_write
from sage.misc.latex_standalone import TikzPicture

from sage.features.gambit import pygambit, gtdraw
from sage.misc.lazy_import import lazy_import
lazy_import('pygambit', ['Game', 'read_efg', 'catalog'], feature=pygambit())
lazy_import('pygambit', 'nash', 'gambit_nash', feature=pygambit())
lazy_import('gtdraw', 'tikz', 'gtdraw_tikz', feature=gtdraw())

# The solvers of :meth:`ExtensiveFormGame.obtain_nash` that compute agent
# equilibria: they work on the extensive form and on it only, and their answer
# is verified as an agent rather than as a Nash equilibrium.
_AGENT_ALGORITHMS = ('enumpure_agent', 'liap_agent')

# The colors :meth:`ExtensiveFormGame.plot` gives the players of a game, in the
# order the game lists them: the first player is red, the second blue, and a
# game with more players than this cycles through the list again.
_PLAYER_COLORS = ['red', 'blue', 'green', 'orange', 'purple', 'brown']

# How :meth:`ExtensiveFormGame.plot` sizes a tree: the vertical distance
# between two levels, as a multiple of the horizontal distance between two
# nodes, how many inches of figure a node is given, and the size in inches past
# which the figure does not grow any further.
_LEVEL_SPACING = 2
_INCHES_PER_NODE = 1
_MAX_FIGSIZE = 30


class GameTreeTikzPicture(TikzPicture):
    r"""
    The TikZ picture of a game tree, as drawn by gtdraw.

    This is the :class:`~sage.misc.latex_standalone.TikzPicture` returned by
    :meth:`ExtensiveFormGame.plot` with ``backend='gtdraw'``, and it inherits
    that class' whole interface: :meth:`~sage.misc.latex_standalone.Standalone.pdf`,
    :meth:`~sage.misc.latex_standalone.Standalone.png`,
    :meth:`~sage.misc.latex_standalone.Standalone.svg`,
    :meth:`~sage.misc.latex_standalone.Standalone.tex` and
    :meth:`~sage.misc.latex_standalone.Standalone.save` write the picture to a
    file, :meth:`~sage.misc.latex_standalone.Standalone.content` returns the
    TikZ source and ``str`` of it the whole standalone LaTeX document.

    It differs from a plain ``TikzPicture`` only in how it displays itself.  A
    ``TikzPicture`` renders inline in a Jupyter notebook but prints its LaTeX
    source at the Sage command line; this class renders in both, compiling the
    picture to a PDF in a temporary file and handing it to the platform viewer,
    the way :meth:`~sage.plot.graphics.Graphics.show` does for an ordinary Sage
    plot.  Displaying it therefore needs a LaTeX installation, while merely
    building it -- and reading its :meth:`content` -- does not.

    EXAMPLES::

        sage: # optional - pygambit gtdraw
        sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
        sage: g = ExtensiveFormGame(players=['Alice', 'Bob'])
        sage: g.append_move(g.root, 'Alice', ['L', 'R'])
        sage: g.set_outcome(g.root.children['L'], 'L', [2, 5])
        sage: g.set_outcome(g.root.children['R'], 'R', [3, 1])
        sage: t = g.plot(backend='gtdraw')
        sage: t
        Game tree TikZ picture
        sage: r'\begin{tikzpicture}' in t.content()
        True
    """
    def _repr_(self):
        r"""
        Return a short string representation.

        Unlike :meth:`sage.misc.latex_standalone.Standalone._repr_`, which
        prints the LaTeX document, this stays on one line: it is the text the
        command line shows next to the launched viewer.  Use ``print(self)``
        for the document and :meth:`content` for the TikZ source.

        EXAMPLES::

            sage: # optional - pygambit gtdraw
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame(players=['Alice', 'Bob'])
            sage: g.append_move(g.root, 'Alice', ['L', 'R'])
            sage: g.set_outcome(g.root.children['L'], 'L', [2, 5])
            sage: g.set_outcome(g.root.children['R'], 'R', [3, 1])
            sage: g.plot(backend='gtdraw')._repr_()
            'Game tree TikZ picture'
        """
        return "Game tree TikZ picture"

    def _rich_repr_(self, display_manager, **kwds):
        r"""
        Return a rich output container for this picture.

        At the Sage command line this compiles the picture to a PDF in a
        temporary file and returns it as an
        :class:`~sage.repl.rich_output.output_graphics.OutputImagePdf`, which
        the command line backend opens in the platform's PDF viewer -- the same
        route an ordinary Sage plot takes.  In a Jupyter notebook the picture is
        rendered inline by
        :meth:`sage.misc.latex_standalone.Standalone._rich_repr_`.  Anywhere
        else, notably in a doctest or a plain script, this returns ``None`` and
        nothing is compiled.

        See :mod:`sage.repl.rich_output` for details.

        INPUT:

        - ``display_manager`` -- the display manager

        - ``**kwds`` -- passed on to
          :meth:`~sage.misc.latex_standalone.Standalone.pdf` (command line) or
          to :meth:`sage.misc.latex_standalone.Standalone._rich_repr_`

        OUTPUT: a rich output container, or ``None``

        EXAMPLES:

        Doctests neither run in a terminal nor in a notebook, so nothing is
        compiled and no LaTeX is needed::

            sage: # optional - pygambit gtdraw
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: from sage.repl.rich_output import get_display_manager
            sage: g = ExtensiveFormGame(players=['Alice', 'Bob'])
            sage: g.append_move(g.root, 'Alice', ['L', 'R'])
            sage: g.set_outcome(g.root.children['L'], 'L', [2, 5])
            sage: g.set_outcome(g.root.children['R'], 'R', [3, 1])
            sage: dm = get_display_manager()
            sage: dm.is_in_terminal()
            False
            sage: g.plot(backend='gtdraw')._rich_repr_(dm) is None
            True
        """
        if display_manager.is_in_terminal():
            types = display_manager.types
            if types.OutputImagePdf not in display_manager.supported_output():
                return None
            from sage.repl.rich_output.buffer import OutputBuffer
            filename = self.pdf(view=False, **kwds)
            return types.OutputImagePdf(OutputBuffer.from_file(filename))
        return super()._rich_repr_(display_manager, **kwds)


class ExtensiveFormGame(SageObject):
    r"""
    An extensive form game [NN2007]_.

    This is a thin wrapper around gambit's extensive form (game tree) ``Game``
    ([Gambit]_).  It stores a ``pygambit`` tree game internally and delegates
    to it, in the same spirit as
    :class:`~sage.game_theory.normal_form_game.NormalFormGame`.

    A game can be built from scratch with the tree-building methods (see below)
    or by wrapping an existing ``pygambit`` game passed to the constructor; use
    :meth:`_gambit_` to recover the underlying gambit game.

    Following gambit, tree nodes carry no labels: a node is the ``pygambit``
    ``Node`` reached by navigating action labels down from the root
    (:meth:`root`).  For example ``g.root.children['L']`` is the child of the
    root reached by the action labeled ``'L'``, and
    ``g.root.children['L'].children['l']`` the grandchild after ``'l'``.  The
    tree-building methods take such nodes and the required **action** labels;
    they are :meth:`add_player`, :meth:`append_move`, :meth:`append_chance_move`,
    :meth:`append_infoset`, :meth:`set_outcome`, :meth:`set_chance_probs`,
    :meth:`insert_move` and :meth:`delete_tree`.

    INPUT:

    - ``generator`` -- the game to wrap; either

      * a ``pygambit`` extensive form (tree) ``Game`` (requires the optional
        gambit package), or

      * ``None`` (default), giving an empty tree ready to build on.

    - ``players`` -- (default: ``None``) list of strategic player labels for a
      new game (a shortcut for calling :meth:`add_player` for each); only
      allowed when ``generator`` is ``None``.

    EXAMPLES:

    Build a two-player game in Sage.  Alice chooses between ``'L'`` and
    ``'R'``, ending the game with payoffs ``[2, 5]`` or ``[3, 1]``::

        sage: # optional - pygambit
        sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
        sage: g = ExtensiveFormGame(players=['Alice', 'Bob'])
        sage: g.append_move(g.root, 'Alice', ['L', 'R'])
        sage: g.set_outcome(g.root.children['L'], 'L', [2, 5])
        sage: g.set_outcome(g.root.children['R'], 'R', [3, 1])
        sage: g
        An extensive form game with 2 players
        sage: sorted(p.label for p in g.players)
        ['Alice', 'Bob']

    The familiar textbook games are built the same way.  In the Prisoner's
    Dilemma, two suspects are questioned separately, each staying silent or
    confessing without knowing what the other does.  Choosing simultaneously
    means, in a tree, that Alice is drawn as moving first while Bob's two nodes
    are joined by :meth:`append_infoset` into a single information set, so that
    he cannot tell them apart.  Gambit maximizes payoffs, so the sentences are
    scored as utilities (higher is better) rather than as the years in prison
    of the :mod:`~sage.game_theory.normal_form_game` example::

        sage: # optional - pygambit
        sage: prisoners = ExtensiveFormGame(players=['Alice', 'Bob'])
        sage: prisoners.append_move(prisoners.root, 'Alice',
        ....:                       ['silent', 'confess'])
        sage: prisoners.append_move(prisoners.root.children['silent'], 'Bob',
        ....:                       ['silent', 'confess'])
        sage: prisoners.append_infoset(prisoners.root.children['confess'],
        ....:                          prisoners.root.children['silent'])
        sage: payoffs = {('silent', 'silent'): [3, 3],
        ....:            ('silent', 'confess'): [0, 5],
        ....:            ('confess', 'silent'): [5, 0],
        ....:            ('confess', 'confess'): [1, 1]}
        sage: for (alice, bob), payoff in payoffs.items():
        ....:     node = prisoners.root.children[alice].children[bob]
        ....:     prisoners.set_outcome(node, '{},{}'.format(alice, bob), payoff)

    As in the normal form, confessing dominates and the unique equilibrium is
    the one both players would rather avoid::

        sage: # optional - pygambit
        sage: actions = prisoners._gambit_().actions
        sage: [[float(eq[a]) for a in actions] for eq in prisoners.obtain_nash()]
        [[0.0, 1.0, 0.0, 1.0]]

    A game built directly in gambit can also be wrapped::

        sage: # optional - pygambit
        sage: from pygambit import Game
        sage: gt = Game.new_tree(players=['Alice', 'Bob'])
        sage: gt.append_move(gt.root, gt.players['Alice'], ['L', 'R'])
        sage: for leaf, label, (a, b) in zip(gt.root.children, ['L', 'R'],
        ....:                                [[2, 5], [3, 1]]):
        ....:     gt.set_outcome(leaf, gt.add_outcome(label, [a, b]))
        sage: ExtensiveFormGame(gt)
        An extensive form game with 2 players

    REFERENCES:

    - [NN2007]_

    - [Gambit]_
    """

    def __init__(self, generator=None, players=None):
        r"""
        Initialize an extensive form game.

        TESTS::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: ExtensiveFormGame()
            An extensive form game with 0 players

        The strategic players can be given up front, as with gambit's
        ``Game.new_tree(players=...)``::

            sage: # optional - pygambit
            sage: g = ExtensiveFormGame(players=['Alice', 'Bob'])
            sage: sorted(p.label for p in g.players)
            ['Alice', 'Bob']
            sage: g.root.is_terminal
            True

        Passing both a game to wrap and ``players`` is an error, as is a
        ``generator`` that is not a gambit game::

            sage: # optional - pygambit
            sage: from pygambit import Game
            sage: ExtensiveFormGame(Game.new_tree(), players=['Alice'])
            Traceback (most recent call last):
            ...
            ValueError: players cannot be given when wrapping an existing game
            sage: ExtensiveFormGame(42)
            Traceback (most recent call last):
            ...
            TypeError: generator must be a gambit extensive form game or None

        A strategic-form (non-tree) gambit game is rejected::

            sage: # optional - pygambit
            sage: import numpy as np
            sage: from pygambit import Game
            sage: tbl = Game.from_arrays(np.array([[1, 0], [0, 1]]),
            ....:                        np.array([[1, 0], [0, 1]]))
            sage: ExtensiveFormGame(tbl)
            Traceback (most recent call last):
            ...
            ValueError: gambit game is not an extensive form (tree) game
        """
        pygambit().require()
        if generator is None:
            self._game = Game.new_tree(players=list(players or []))
        elif Game is not None and isinstance(generator, Game):
            if players is not None:
                raise ValueError("players cannot be given when wrapping an "
                                 "existing game")
            self._gambit_game(generator)
        else:
            raise TypeError("generator must be a gambit extensive form game "
                            "or None")

    def _check_actions(self, actions):
        r"""
        Validate ``actions`` as the action labels of a single move.

        There must be at least one, and each must be a nonempty string,
        distinct from the others.  Distinctness matters because nodes are
        navigated by action label (``node.children['<label>']``); ``pygambit``
        itself only warns on duplicates.  Returns ``list(actions)``.

        TESTS::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame(players=['Alice'])
            sage: g.append_move(g.root, 'Alice', ['L', 'L'])
            Traceback (most recent call last):
            ...
            ValueError: duplicate action label 'L'
            sage: g.append_move(g.root, 'Alice', ['L', ''])
            Traceback (most recent call last):
            ...
            ValueError: action label must be a nonempty string, not ''
            sage: g.append_move(g.root, 'Alice', [])
            Traceback (most recent call last):
            ...
            ValueError: a move needs at least one action
        """
        actions = list(actions)
        if not actions:
            raise ValueError("a move needs at least one action")
        seen = set()
        for action in actions:
            if not isinstance(action, str) or not action:
                raise ValueError("action label must be a nonempty string, "
                                 "not {0!r}".format(action))
            if action in seen:
                raise ValueError("duplicate action label {0!r}".format(action))
            seen.add(action)
        return actions

    def _gambit_game(self, game):
        r"""
        Populate this game from a gambit extensive form ``Game``, in place.

        This stores the gambit tree game, replacing any game already wrapped.
        It is the inverse of :meth:`_gambit_`.

        This requires the optional gambit package.

        TESTS::

            sage: # optional - pygambit
            sage: from pygambit import Game
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: gt = Game.new_tree(players=['Alice', 'Bob'])
            sage: gt.append_move(gt.root, gt.players['Alice'], ['L', 'R'])
            sage: g = ExtensiveFormGame()
            sage: g._gambit_game(gt); g
            An extensive form game with 2 players
        """
        pygambit().require()
        if not game.is_tree:
            raise ValueError("gambit game is not an extensive form (tree) game")
        self._game = game

    def _gambit_(self):
        r"""
        Return the underlying gambit extensive form ``Game``.

        This requires the optional gambit package.

        EXAMPLES::

            sage: # optional - pygambit
            sage: from pygambit import Game
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: gt = Game.new_tree(players=['Alice', 'Bob'])
            sage: gt.append_move(gt.root, gt.players['Alice'], ['L', 'R'])
            sage: g = ExtensiveFormGame(gt)
            sage: g._gambit_() is gt
            True

        A freshly created game wraps an empty tree::

            sage: # optional - pygambit
            sage: ExtensiveFormGame()._gambit_().root.is_terminal
            True
        """
        pygambit().require()
        return self._game

    def _repr_(self) -> str:
        r"""
        Return a concise description of the game.

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame(players=['Alice'])
            sage: g.append_move(g.root, 'Alice', ['L', 'R'])
            sage: g
            An extensive form game with 1 player
            sage: ExtensiveFormGame()
            An extensive form game with 0 players
        """
        n = len(self._game.players)
        return "An extensive form game with {0} player{1}".format(
            n, "" if n == 1 else "s")

    @property
    def root(self):
        r"""
        The root node of the tree, i.e. where to start building and navigating.

        This is the underlying gambit ``Node``; reach any other node from it by
        indexing children with action labels, e.g. ``g.root.children['L']``.

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: ExtensiveFormGame().root.is_terminal
            True
            sage: g = ExtensiveFormGame(players=['Alice'])
            sage: g.append_move(g.root, 'Alice', ['L', 'R'])
            sage: [child.is_terminal for child in g.root.children]
            [True, True]
        """
        return self._gambit_().root

    @property
    def players(self):
        r"""
        The (strategic) players of the game, as gambit players.

        The chance player is not included (gambit keeps it separate).

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame(players=['Alice', 'Bob'])
            sage: sorted(p.label for p in g.players)
            ['Alice', 'Bob']
        """
        return self._gambit_().players

    @property
    def infosets(self):
        r"""
        The information sets of the game, as gambit information sets.

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame(players=['Alice'])
            sage: g.append_move(g.root, 'Alice', ['L', 'R'])
            sage: len(g.infosets)
            1
        """
        return self._gambit_().infosets

    @property
    def outcomes(self):
        r"""
        The outcomes of the game, as gambit outcomes.

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame(players=['Alice', 'Bob'])
            sage: g.append_move(g.root, 'Alice', ['L', 'R'])
            sage: g.set_outcome(g.root.children['L'], 'L', [2, 5])
            sage: g.set_outcome(g.root.children['R'], 'R', [3, 1])
            sage: len(g.outcomes)
            2
        """
        return self._gambit_().outcomes

    @property
    def is_perfect_recall(self):
        r"""
        Whether the game is a game of perfect recall.

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame(players=['Alice'])
            sage: g.append_move(g.root, 'Alice', ['L', 'R'])
            sage: g.is_perfect_recall
            True
        """
        return self._gambit_().is_perfect_recall

    def add_player(self, label):
        r"""
        Add a (strategic) player to the game.

        INPUT:

        - ``label`` -- string; a label identifying the player.  Following
          gambit it must be nonempty and distinct from the labels of the
          game's other players.

        OUTPUT: the label of the new player

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame()
            sage: g.add_player('Alice')
            'Alice'
            sage: [p.label for p in g.players]
            ['Alice']

        TESTS:

        Gambit rejects an empty label, and a label already used by another
        player of the game::

            sage: # optional - pygambit
            sage: g.add_player('')
            Traceback (most recent call last):
            ...
            ValueError: Player label must not be empty
            sage: g.add_player('Alice')
            Traceback (most recent call last):
            ...
            ValueError: Player label must be unique within the game
        """
        pygambit().require()
        return self._gambit_().add_player(label).label

    def append_move(self, node, player, actions):
        r"""
        Append a move for ``player`` at the terminal node ``node``.

        This makes ``node`` a decision node at which ``player`` chooses among
        ``actions``, creating one child per action.  The child reached by an
        action is then ``node.children['<action>']``.

        INPUT:

        - ``node`` -- the (terminal) node to move at (see :meth:`root`)

        - ``player`` -- string; the label of the player who moves

        - ``actions`` -- list of strings; the labels of the available actions,
          which must be distinct and nonempty

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame(players=['Alice', 'Bob'])
            sage: g.append_move(g.root, 'Alice', ['L', 'R'])
            sage: g.append_move(g.root.children['L'], 'Bob', ['l', 'r'])
            sage: g.root.children['L'].children['l'].is_terminal
            True
        """
        pygambit().require()
        actions = self._check_actions(actions)
        self._gambit_().append_move(node, player, actions)

    def append_chance_move(self, node, actions, probs=None):
        r"""
        Append a chance (nature) move at the terminal node ``node``.

        This makes ``node`` a chance node with one branch per action; if
        ``probs`` is given it sets the probabilities of those branches.

        INPUT:

        - ``node`` -- the (terminal) node to move at (see :meth:`root`)

        - ``actions`` -- list of strings; the labels of the chance branches,
          which must be distinct and nonempty

        - ``probs`` -- (default: ``None``) list of branch probabilities; each
          may be a string such as ``'1/2'`` or a Sage rational

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame(players=['Alice'])
            sage: g.append_chance_move(g.root, ['H', 'T'], probs=['1/2', '1/2'])
            sage: g.root.infoset.is_chance
            True
        """
        pygambit().require()
        game = self._gambit_()
        actions = self._check_actions(actions)
        game.append_move(node, game.players.chance, actions)
        if probs is not None:
            game.set_chance_probs(node.infoset, probs)

    def append_infoset(self, node, like):
        r"""
        Append a move at ``node`` in the same information set as ``like``.

        This makes ``node`` a decision node belonging to the information set of
        the node ``like``: the same player moves with the same actions, and the
        two nodes cannot be told apart (imperfect information).

        INPUT:

        - ``node`` -- the (terminal) node to move at (see :meth:`root`)

        - ``like`` -- a decision node whose information set to join

        EXAMPLES:

        Bob moves after Alice without observing her choice::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame(players=['Alice', 'Bob'])
            sage: g.append_move(g.root, 'Alice', ['L', 'R'])
            sage: g.append_move(g.root.children['L'], 'Bob', ['l', 'r'])
            sage: g.append_infoset(g.root.children['R'], g.root.children['L'])
            sage: len(g.infosets)
            2
        """
        pygambit().require()
        self._gambit_().append_infoset(node, like.infoset)

    def set_outcome(self, node, outcome, payoffs=None):
        r"""
        Set the outcome awarded at the terminal node ``node``.

        Given ``payoffs``, a new outcome awarding them is created; given none,
        ``outcome`` is an outcome the game already has, so that the same
        outcome can be awarded at several nodes.

        INPUT:

        - ``node`` -- a terminal node (see :meth:`root`)

        - ``outcome`` -- when ``payoffs`` is given, a string labelling the new
          outcome; following gambit it must be nonempty and distinct from the
          labels of the game's other outcomes.  Otherwise an outcome the game
          already has, either as a gambit outcome (see :meth:`outcomes`) or by
          its label, or ``None`` to award no outcome at ``node``.

        - ``payoffs`` -- (default: ``None``) a list of payoffs, one per player
          (in the order of :meth:`players`), awarded by the new outcome
          labelled ``outcome``

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame(players=['Alice', 'Bob'])
            sage: g.append_move(g.root, 'Alice', ['L', 'R'])
            sage: g.set_outcome(g.root.children['L'], 'L', [2, 5])
            sage: float(g.root.children['L'].outcome['Alice'])
            2.0

        The same outcome can be awarded at several nodes, naming it either by
        its label or as an outcome of the game.  Alice below wins the same 1
        whichever way the coin lands, and the game has a single outcome rather
        than one per leaf::

            sage: # optional - pygambit
            sage: coin = ExtensiveFormGame(players=['Alice', 'Bob'])
            sage: coin.append_chance_move(coin.root, ['heads', 'tails'])
            sage: heads, tails = coin.root.children
            sage: coin.set_outcome(heads, 'Alice wins', [1, -1])
            sage: coin.set_outcome(tails, 'Alice wins')
            sage: coin.set_outcome(tails, heads.outcome)
            sage: len(coin.outcomes)
            1

        Passing ``None`` awards no outcome at ``node``::

            sage: # optional - pygambit
            sage: coin.set_outcome(tails, None)
            sage: bool(tails.outcome)
            False

        TESTS:

        Gambit rejects an empty label, and a label already used by another
        outcome of the game::

            sage: # optional - pygambit
            sage: g.set_outcome(g.root.children['R'], '', [3, 1])
            Traceback (most recent call last):
            ...
            ValueError: Outcome label must not be empty
            sage: g.set_outcome(g.root.children['R'], 'L', [3, 1])
            Traceback (most recent call last):
            ...
            ValueError: Outcome label must be unique within the game

        Awarding an existing outcome requires one the game actually has::

            sage: # optional - pygambit
            sage: g.set_outcome(g.root.children['R'], 'nope')
            Traceback (most recent call last):
            ...
            KeyError: "set_outcome(): no outcome with label 'nope'"
        """
        pygambit().require()
        game = self._gambit_()
        if payoffs is not None:
            outcome = game.add_outcome(outcome, list(payoffs))
        game.set_outcome(node, outcome)

    def set_chance_probs(self, node, probs):
        r"""
        Set the branch probabilities of the chance move at ``node``.

        INPUT:

        - ``node`` -- a node at which a chance move sits (see
          :meth:`append_chance_move`)

        - ``probs`` -- list of branch probabilities; each may be a string such
          as ``'1/3'`` or a Sage rational

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame(players=['Alice'])
            sage: g.append_chance_move(g.root, ['H', 'T'])
            sage: g.set_chance_probs(g.root, ['1/3', '2/3'])
            sage: float(list(g.root.infoset.actions)[0].prob)
            0.3333333333333333
        """
        pygambit().require()
        self._gambit_().set_chance_probs(node.infoset, probs)

    def insert_move(self, node, player, actions):
        r"""
        Insert a new move for ``player`` immediately above ``node``.

        A new decision node with the given ``actions`` is inserted immediately
        above ``node``; ``node`` (keeping its subtree) becomes the first child
        of the new node, reached by the first action, and the remaining actions
        lead to fresh terminal siblings.

        INPUT:

        - ``node`` -- the node to insert the move above (see :meth:`root`).  It
          becomes the first child of the new move.

        - ``player`` -- string; the label of the player who moves

        - ``actions`` -- list of strings; the labels of the new move's actions,
          which must be distinct and nonempty

        .. NOTE::

            Inserting above the root makes the new node the root of the tree,
            so :meth:`root` then returns it and the former root becomes its
            first child.

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame(players=['Alice'])
            sage: g.append_move(g.root, 'Alice', ['L', 'R'])
            sage: former_root = g.root
            sage: g.insert_move(g.root, 'Alice', ['a', 'b', 'c'])
            sage: len(g.root.children)   # a new 3-action move above the old root
            3
            sage: g.root.children['a'] == former_root   # old root reached by 'a'
            True
            sage: g.root.children['b'].is_terminal       # a fresh sibling
            True
        """
        pygambit().require()
        actions = self._check_actions(actions)
        self._gambit_().insert_move(node, player, len(actions))
        for action, label in zip(node.parent.infoset.actions, actions):
            action.label = label

    def delete_tree(self, node):
        r"""
        Delete the subtree below ``node``, making it a terminal node.

        INPUT:

        - ``node`` -- the node to prune (see :meth:`root`)

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame(players=['Alice'])
            sage: g.append_move(g.root, 'Alice', ['L', 'R'])
            sage: g.delete_tree(g.root)
            sage: g.root.is_terminal
            True
        """
        pygambit().require()
        self._gambit_().delete_tree(node)

    def to_efg(self):
        r"""
        Return the game as a string in gambit's extensive-form ``.efg`` format.

        This is a thin passthrough to gambit's serializer.  It is handy for
        printing the full structure of the tree, which the concise
        :meth:`_repr_` does not show; use :meth:`save_efg` to write it to a
        file instead.

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame(players=['Alice', 'Bob'])
            sage: g.append_move(g.root, 'Alice', ['L', 'R'])
            sage: g.set_outcome(g.root.children['L'], 'L', [2, 5])
            sage: g.set_outcome(g.root.children['R'], 'R', [3, 1])
            sage: print(g.to_efg().splitlines()[0])
            EFG 2 R "Untitled extensive game" { "Alice" "Bob" }
        """
        pygambit().require()
        return self._gambit_().to_efg()

    def plot(self, backend='sage', **kwargs):
        r"""
        Plot the game tree.

        Two drawing backends are available, selected with ``backend``.

        The ``'gtdraw'`` backend hands the underlying gambit game to
        `gtdraw <https://www.gambit-project.org/gtdraw/>`_, the game tree
        drawing tool of the gambit project, which produces publication-quality
        TikZ pictures.  The result is a :class:`GameTreeTikzPicture`, a
        :class:`~sage.misc.latex_standalone.TikzPicture` like the ones
        :meth:`~sage.graphs.generic_graph.GenericGraph.tikz` returns for
        graphs.  Displaying it compiles the picture to a PDF in a temporary
        file and opens it in the platform's viewer, just as showing a Sage
        plot does, and renders it inline in a Jupyter notebook; that step needs
        a LaTeX installation providing ``pdflatex`` and TikZ, see the `gtdraw
        installation guide
        <https://www.gambit-project.org/gtdraw/installation/>`_.  Building the
        picture and reading its TikZ source with
        :meth:`~sage.misc.latex_standalone.Standalone.content` need no LaTeX,
        only the optional gtdraw package.

        The default ``'sage'`` backend draws the tree with Sage's own graph
        plotting and needs nothing beyond Sage.  Labels are written beside the
        nodes: the moving player above a decision node, ``'Chance'`` above a
        chance node, the payoffs below a terminal node, and the action labels
        along the branches, a branch out of a chance node also carrying its
        probability, as it does in gtdraw.  A node is colored by whose move it
        is -- the first player of the game red, the second blue, chance gray
        and a terminal node white -- and, when the player has more than one
        information set, the number after their name says which of them the
        node belongs to, so imperfect information is visible: the nodes a
        player cannot tell apart carry the same number.

        The figure grows with the tree, up to a point, so that the nodes keep
        their spacing however large the game is; past that size the labels
        shrink instead.  Both can be overridden with the ``figsize`` and
        ``label_fontsize`` keywords.  The result is a
        :class:`~sage.plot.graphics.Graphics` object, so it composes with the
        rest of Sage's plotting.

        INPUT:

        - ``backend`` -- string (default: ``'sage'``); which backend to draw
          with, one of

          * ``'sage'`` -- draw the tree with Sage's graph plotting

          * ``'gtdraw'`` -- draw the tree with gtdraw (this requires the
            optional gtdraw package)

        - ``**kwargs`` -- passed on to the selected backend: to
          :meth:`~sage.graphs.generic_graph.GenericGraph.plot` for ``'sage'``
          (for instance ``figsize`` or ``vertex_size``), and to ``gtdraw.tikz``
          for ``'gtdraw'``.  Useful gtdraw options are ``horizontal``,
          ``mirror``, ``scale_factor``, ``node_size``, ``hide_action_labels``,
          ``iset_fill``, ``font_size``, ``legend_position`` and
          ``color_scheme`` (one of ``'default'``, ``'gambit'``,
          ``'distinctipy'`` and ``'colorblind'``); see the `gtdraw
          documentation <https://www.gambit-project.org/gtdraw/>`_ for the
          full list.

        OUTPUT:

        For ``backend='sage'``, a :class:`~sage.plot.graphics.Graphics` object.
        For ``backend='gtdraw'``, a :class:`GameTreeTikzPicture`.

        EXAMPLES:

        Take a small game of imperfect information: Alice chooses ``'L'`` or
        ``'R'``, and Bob then chooses ``'l'`` or ``'r'`` without having seen
        her choice, so his two nodes share an information set::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame(players=['Alice', 'Bob'])
            sage: g.append_move(g.root, 'Alice', ['L', 'R'])
            sage: g.append_move(g.root.children['L'], 'Bob', ['l', 'r'])
            sage: g.append_infoset(g.root.children['R'], g.root.children['L'])
            sage: for a, payoff in zip(['l', 'r'], [[3, 1], [0, 0]]):
            ....:     g.set_outcome(g.root.children['L'].children[a],
            ....:                   'L,{}'.format(a), payoff)
            sage: for a, payoff in zip(['l', 'r'], [[0, 0], [1, 3]]):
            ....:     g.set_outcome(g.root.children['R'].children[a],
            ....:                   'R,{}'.format(a), payoff)

        The default ``'sage'`` backend needs nothing beyond Sage and gives a
        :class:`~sage.plot.graphics.Graphics` object, which can be shown, saved
        or combined with other graphics.  The ``'gtdraw'`` backend, shown
        further below, needs the optional gtdraw package and gives a
        :class:`GameTreeTikzPicture` instead::

            sage: # optional - pygambit
            sage: from sage.plot.graphics import Graphics
            sage: isinstance(g.plot(backend='sage'), Graphics)
            True

        Being the default, it is what plotting without a ``backend`` gives, so
        that drawing a tree never requires an optional package::

            sage: # optional - pygambit
            sage: isinstance(g.plot(), Graphics)
            True

        Keyword arguments reach Sage's graph plotting::

            sage: # optional - pygambit
            sage: isinstance(g.plot(backend='sage', figsize=8, vertex_size=400),
            ....:            Graphics)
            True

        The nodes are colored by whose move it is: Alice's node is red, Bob's
        two are blue, and the four terminal nodes are white::

            sage: # optional - pygambit
            sage: from sage.plot.scatter_plot import ScatterPlot
            sage: nodes = [p for p in g.plot(backend='sage')
            ....:          if isinstance(p, ScatterPlot)][0]
            sage: sorted(set(nodes.options()['facecolor']))
            ['blue', 'red', 'white']
            sage: nodes.options()['facecolor'].count('blue')
            2

        Bob has a single information set here, so his nodes are labeled with
        his name alone::

            sage: # optional - pygambit
            sage: drawn = [p.string for p in g.plot(backend='sage')
            ....:          if hasattr(p, 'string')]
            sage: sorted(set(s for s in drawn if 'Bob' in s))
            ['Bob']

        A branch out of a chance node is labeled with its probability as well
        as its action, so a fair coin flip that Alice then bets on reads
        ``'H (1/2)'`` and ``'T (1/2)'``::

            sage: # optional - pygambit
            sage: coin = ExtensiveFormGame(players=['Alice'])
            sage: coin.append_chance_move(coin.root, ['H', 'T'],
            ....:                         probs=['1/2', '1/2'])
            sage: for side in ['H', 'T']:
            ....:     coin.append_move(coin.root.children[side], 'Alice',
            ....:                      ['bet', 'fold'])
            sage: drawn = [p.string for p in coin.plot(backend='sage')
            ....:          if hasattr(p, 'string')]
            sage: sorted(s for s in drawn if '/' in s)
            ['H (1/2)', 'T (1/2)']
            sage: sorted(set(s for s in drawn if s in ['bet', 'fold']))
            ['bet', 'fold']

        The chance node is gray, and Alice, who moves at two nodes she can tell
        apart, has her two information sets numbered::

            sage: # optional - pygambit
            sage: nodes = [p for p in coin.plot(backend='sage')
            ....:          if isinstance(p, ScatterPlot)][0]
            sage: sorted(set(nodes.options()['facecolor']))
            ['lightgray', 'red', 'white']
            sage: sorted(s for s in drawn if s.startswith('Alice'))
            ['Alice 1', 'Alice 2']

        The gtdraw backend draws the same tree as TikZ.  Its
        :meth:`~sage.misc.latex_standalone.Standalone.content` is the source,
        which mentions the players and the tree itself::

            sage: # optional - pygambit gtdraw
            sage: picture = g.plot(backend='gtdraw')
            sage: picture
            Game tree TikZ picture
            sage: tikz = picture.content()
            sage: print(tikz.splitlines()[0])
            % TikZ code with built-in styling for game trees
            sage: r'\begin{tikzpicture}' in tikz
            True
            sage: r'\def\playerone{Alice}' in tikz
            True
            sage: r'\def\playertwo{Bob}' in tikz
            True

        Its many options control the look of the picture.  Growing the tree
        left to right instead of top to bottom, shading the information sets
        and picking a colorblind-safe palette all change the source::

            sage: # optional - pygambit gtdraw
            sage: g.plot(backend='gtdraw', horizontal=True).content() != tikz
            True
            sage: g.plot(backend='gtdraw', iset_fill=True).content() != tikz
            True
            sage: colorblind = g.plot(backend='gtdraw', color_scheme='colorblind')
            sage: print([line for line in colorblind.content().splitlines()
            ....:        if 'definecolor{p1rgb}' in line][0])
            \definecolor{p1rgb}{RGB}{0,127,255}

        Suppressing the action labels drops ``'L'`` from the picture::

            sage: # optional - pygambit gtdraw
            sage: r'{L\strut}' in tikz
            True
            sage: hidden = g.plot(backend='gtdraw', hide_action_labels=True)
            sage: r'{L\strut}' in hidden.content()
            False

        Options can of course be combined::

            sage: # optional - pygambit gtdraw
            sage: fancy = g.plot(backend='gtdraw', horizontal=True,
            ....:                color_scheme='gambit', scale_factor=1.5,
            ....:                node_size=2.0, font_size='small')
            sage: r'\begin{tikzpicture}' in fancy.content()
            True

        Games from the gambit catalog can be drawn the same way.  Catalog slugs
        are full paths, which :meth:`gambit_catalog_games` lists::

            sage: # optional - pygambit
            sage: 'books/myerson1991/fig2_1' in list(
            ....:     ExtensiveFormGame.gambit_catalog_games()['Game'])
            True

        Myerson's simple poker game is a two-player tree with three information
        sets, and it draws with either backend::

            sage: # optional - pygambit
            sage: g = ExtensiveFormGame.load_from_gambit_catalog(
            ....:     'books/myerson1991/fig2_1')
            sage: g
            An extensive form game with 2 players
            sage: len(g.infosets)
            3
            sage: isinstance(g.plot(backend='sage'), Graphics)
            True

        ::

            sage: # optional - pygambit gtdraw
            sage: r'\begin{tikzpicture}' in g.plot(backend='gtdraw').content()
            True

        So does the stripped-down poker game of Reiley et al., here drawn
        sideways::

            sage: # optional - pygambit
            sage: h = ExtensiveFormGame.load_from_gambit_catalog(
            ....:     'journals/other/reiley2008/fig1')
            sage: sorted(p.label for p in h.players)
            ['Professor', 'Student']

        ::

            sage: # optional - pygambit gtdraw
            sage: sideways = h.plot(backend='gtdraw', horizontal=True)
            sage: r'\begin{tikzpicture}' in sideways.content()
            True

        Larger trees are where gtdraw earns its keep.  Kuhn poker is not in the
        catalog, but it is quick to build: chance deals Alice and Bob one card
        each from ``J``, ``Q``, ``K``, and each player then sees only their own
        card.  We label a chance branch by the pair of cards dealt, so that
        ``'JQ'`` means Alice holds the jack and Bob the queen::

            sage: # optional - pygambit
            sage: kuhn = ExtensiveFormGame(players=['Alice', 'Bob'])
            sage: deals = ['JQ', 'JK', 'QJ', 'QK', 'KJ', 'KQ']
            sage: kuhn.append_chance_move(kuhn.root, deals, probs=['1/6'] * 6)

        Alice opens by checking or betting.  She knows only her own card, so
        the two deals giving her the same card are in one information set::

            sage: # optional - pygambit
            sage: opening = {}
            sage: for d in deals:
            ....:     node = kuhn.root.children[d]
            ....:     if d[0] in opening:
            ....:         kuhn.append_infoset(node, opening[d[0]])
            ....:     else:
            ....:         kuhn.append_move(node, 'Alice', ['check', 'bet'])
            ....:         opening[d[0]] = node

        Bob replies, seeing his own card and Alice's action but not her card::

            sage: # optional - pygambit
            sage: reply = {}
            sage: for d in deals:
            ....:     for a in ['check', 'bet']:
            ....:         node = kuhn.root.children[d].children[a]
            ....:         if (d[1], a) in reply:
            ....:             kuhn.append_infoset(node, reply[(d[1], a)])
            ....:         else:
            ....:             acts = ['check', 'bet'] if a == 'check' else ['fold', 'call']
            ....:             kuhn.append_move(node, 'Bob', acts)
            ....:             reply[(d[1], a)] = node

        Alice, having checked, may still fold or call when Bob bets::

            sage: # optional - pygambit
            sage: facing_bet = {}
            sage: for d in deals:
            ....:     node = kuhn.root.children[d].children['check'].children['bet']
            ....:     if d[0] in facing_bet:
            ....:         kuhn.append_infoset(node, facing_bet[d[0]])
            ....:     else:
            ....:         kuhn.append_move(node, 'Alice', ['fold', 'call'])
            ....:         facing_bet[d[0]] = node

        Finally the payoffs: at a showdown the higher card takes the pot, and
        a fold hands one unit to the other player::

            sage: # optional - pygambit
            sage: for d in deals:
            ....:     w = 1 if 'JQK'.index(d[0]) > 'JQK'.index(d[1]) else -1
            ....:     n = kuhn.root.children[d]
            ....:     kuhn.set_outcome(n.children['check'].children['check'],
            ....:                      '{} check,check'.format(d), [w, -w])
            ....:     kuhn.set_outcome(n.children['check'].children['bet'].children['fold'],
            ....:                      '{} check,bet,fold'.format(d), [-1, 1])
            ....:     kuhn.set_outcome(n.children['check'].children['bet'].children['call'],
            ....:                      '{} check,bet,call'.format(d), [2 * w, -2 * w])
            ....:     kuhn.set_outcome(n.children['bet'].children['fold'],
            ....:                      '{} bet,fold'.format(d), [1, -1])
            ....:     kuhn.set_outcome(n.children['bet'].children['call'],
            ....:                      '{} bet,call'.format(d), [2 * w, -2 * w])

        Each player ends up with six information sets, three for each card they
        might hold, and the game has perfect recall::

            sage: # optional - pygambit
            sage: len(kuhn.infosets)
            12
            sage: len(kuhn._gambit_().players['Alice'].infosets)
            6
            sage: kuhn.is_perfect_recall
            True

        Its unique equilibrium gives Alice the classical value `-1/18`::

            sage: # optional - pygambit
            sage: eq = kuhn.obtain_nash()[0]
            sage: QQ(eq.payoff(kuhn._gambit_().players['Alice']))
            -1/18

        Both backends draw it; a tree this wide is much more readable drawn
        sideways, and shading the information sets makes the deal structure
        stand out::

            sage: # optional - pygambit
            sage: isinstance(kuhn.plot(backend='sage'), Graphics)
            True

        Its figure is much wider than that of the small game above, which is
        what keeps the fifty-odd nodes of this one apart::

            sage: # optional - pygambit
            sage: kuhn.plot(backend='sage').get_axes_range()['xmax'] > 20
            True
            sage: (kuhn.plot(backend='sage')._extra_kwds['figsize'][0]
            ....:  > g.plot(backend='sage')._extra_kwds['figsize'][0])
            True
            sage: g.plot(backend='sage', figsize=8)._extra_kwds['figsize']
            8

        Each of Bob's six information sets is numbered, and the two deals
        leaving him the same card after the same move of Alice's -- which he
        cannot tell apart -- carry the same number::

            sage: # optional - pygambit
            sage: drawn = [p.string for p in kuhn.plot(backend='sage')
            ....:          if hasattr(p, 'string')]
            sage: sorted(set(s for s in drawn if s.startswith('Bob')))
            ['Bob 1', 'Bob 2', 'Bob 3', 'Bob 4', 'Bob 5', 'Bob 6']
            sage: sorted(s for s in drawn if s.startswith('Bob')).count('Bob 1')
            2

        ::

            sage: # optional - pygambit gtdraw
            sage: wide = kuhn.plot(backend='gtdraw', horizontal=True,
            ....:                  iset_fill=True, color_scheme='colorblind')
            sage: r'\begin{tikzpicture}' in wide.content()
            True

        To write a picture to a file rather than display it, use the
        :class:`~sage.misc.latex_standalone.TikzPicture` methods.  These
        compile the TikZ source, so they need a LaTeX installation, plus
        ImageMagick for ``png`` and ``pdftocairo`` for ``svg``::

            sage: _ = wide.pdf(view=False)   # long time (2s), optional - pygambit gtdraw latex
            sage: wide.pdf('kuhn.pdf')                   # not tested
            sage: wide.png('kuhn.png', density=300)      # not tested
            sage: wide.svg('kuhn.svg')                   # not tested
            sage: wide.tex('kuhn.tex')                   # not tested

        TESTS:

        A game with no moves at all is a single terminal node, and draws as
        one::

            sage: # optional - pygambit
            sage: isinstance(ExtensiveFormGame().plot(backend='sage'), Graphics)
            True

        An unknown backend is rejected, and is caught before the optional
        gtdraw package is looked for::

            sage: # optional - pygambit
            sage: g.plot(backend='bogus')
            Traceback (most recent call last):
            ...
            ValueError: unknown backend 'bogus'; must be 'sage' or 'gtdraw'

        The gtdraw backend returns a :class:`GameTreeTikzPicture`, and building
        it compiles nothing, so it needs no LaTeX::

            sage: # optional - pygambit gtdraw
            sage: from sage.game_theory.extensive_form_game import GameTreeTikzPicture
            sage: isinstance(g.plot(backend='gtdraw'), GameTreeTikzPicture)
            True
            sage: str(g.plot(backend='gtdraw')).count(r'\begin{document}')
            1

        Before drawing the tree, gtdraw writes the game out in its own ``.ef``
        format.  That file is a temporary one, so plotting leaves nothing
        behind in the current directory::

            sage: # optional - pygambit gtdraw
            sage: import os
            sage: from tempfile import TemporaryDirectory
            sage: t = ExtensiveFormGame(players=['Alice', 'Bob'])
            sage: t.append_move(t.root, 'Alice', ['L', 'R'])
            sage: t.set_outcome(t.root.children['L'], 'L', [2, 5])
            sage: t.set_outcome(t.root.children['R'], 'R', [3, 1])
            sage: with TemporaryDirectory() as d:
            ....:     cwd = os.getcwd()
            ....:     os.chdir(d)
            ....:     try:
            ....:         tikz = t.plot(backend='gtdraw').content()
            ....:     finally:
            ....:         os.chdir(cwd)
            ....:     os.listdir(d)
            []

        The source names that file the way gtdraw itself would, so no
        temporary path appears in the picture::

            sage: # optional - pygambit gtdraw
            sage: [line for line in tikz.splitlines()
            ....:  if line.startswith('% Game tree content')]
            ['% Game tree content from Untitled extensive game.ef']

        Asking gtdraw to keep the file, with its own ``save_to`` option, still
        works::

            sage: # optional - pygambit gtdraw
            sage: path = tmp_filename(ext='.ef')
            sage: _ = t.plot(backend='gtdraw', save_to=path)
            sage: os.path.getsize(path) > 0
            True
        """
        pygambit().require()
        if backend == 'gtdraw':
            gtdraw().require()
            game = self._gambit_()
            save_to = kwargs.pop('save_to', None)
            if save_to is None:
                # gtdraw first writes the game out in its own .ef format, and
                # without a path it drops that file in the current directory.
                # Hand it a temporary one instead, then put the name it would
                # have used back into the source, where it appears in a
                # comment, so that no temporary path leaks into the picture.
                with TemporaryDirectory() as tmp_dir:
                    ef_file = os.path.join(tmp_dir, 'game.ef')
                    tikz = gtdraw_tikz(game, save_to=ef_file, **kwargs)
                tikz = tikz.replace(ef_file, game.title + '.ef')
            else:
                tikz = gtdraw_tikz(game, save_to=save_to, **kwargs)
            # The preamble mirrors the standalone document gtdraw itself wraps
            # its TikZ code in, so that the picture matches gtdraw's own output.
            return GameTreeTikzPicture(
                tikz,
                standalone_config=['border=10pt'],
                usepackage=['graphicx'],
                usetikzlibrary=['shapes', 'arrows.meta'],
                macros=[r'\IfFileExists{newpxtext.sty}'
                        r'{\usepackage{newpxtext,newpxmath}}{}',
                        r'\linespread{1.10}'])
        if backend != 'sage':
            raise ValueError("unknown backend {0!r}; must be 'sage' or "
                             "'gtdraw'".format(backend))

        from sage.graphs.digraph import DiGraph
        from sage.plot.text import text

        game = self._gambit_()
        graph = DiGraph()
        labels = {}
        colors = {}
        parent_of = {}
        counter = [0]

        color_of_player = {p.label: _PLAYER_COLORS[i % len(_PLAYER_COLORS)]
                           for i, p in enumerate(game.players)}

        # Number the information sets of each player, so that nodes the player
        # cannot tell apart carry the same number.  A player with a single
        # information set needs no number.  Gambit hands back a fresh wrapper
        # object on every access, but two wrappers of one information set
        # compare and hash alike, so a dictionary keyed by them is sound.
        numbers = {}
        for player in list(game.players) + [game.players.chance]:
            player_infosets = list(player.infosets)
            if len(player_infosets) > 1:
                for number, infoset in enumerate(player_infosets, start=1):
                    numbers[infoset] = number

        def visit(node):
            index = counter[0]
            counter[0] += 1
            if node.is_terminal:
                # Following gambit, ``node.outcome`` is a predicate evaluated
                # on demand: it is falsy, rather than ``None``, at a terminal
                # node carrying no outcome.
                if node.outcome:
                    payoffs = ", ".join(str(node.outcome[p]) for p in game.players)
                    labels[index] = "({0})".format(payoffs)
                else:
                    labels[index] = ""
                colors[index] = "white"
            else:
                infoset = node.infoset
                if infoset.is_chance:
                    labels[index] = "Chance"
                    colors[index] = "lightgray"
                else:
                    labels[index] = infoset.player.label
                    colors[index] = color_of_player[infoset.player.label]
                if infoset in numbers:
                    labels[index] += " {0}".format(numbers[infoset])
                # The tree layout lays the children out in the reverse of the
                # order they are added in, so add them reversed to draw the
                # moves left to right in the order the game lists them.
                branches = list(zip(node.children, infoset.actions))
                for child, action in reversed(branches):
                    child_index = visit(child)
                    parent_of[child_index] = index
                    if infoset.is_chance:
                        # As gtdraw does, a branch out of a chance node carries
                        # its probability next to the action label.
                        edge_label = "{0} ({1})".format(action.label, action.prob)
                    else:
                        edge_label = action.label
                    graph.add_edge(index, child_index, edge_label)
            return index

        root_index = visit(game.root)
        graph.add_vertex(root_index)   # a game with no moves is a lone node

        vertex_colors = {}
        for index, color in colors.items():
            vertex_colors.setdefault(color, []).append(index)

        # Lay the tree out by hand, rather than with ``layout='tree'``, so that
        # the picture can be sized from it.  The levels are spread out
        # vertically to leave room for the action labels along the branches.
        pos = graph.layout_tree(tree_root=root_index, tree_orientation='down')
        pos = {v: (x, _LEVEL_SPACING * y) for v, (x, y) in pos.items()}
        xs = [x for x, _ in pos.values()]
        ys = [y for _, y in pos.values()]

        # Graph plotting draws with an aspect ratio of 1, so the figure has to
        # grow with the tree in both directions at the same rate if the nodes
        # are to keep their spacing however large the tree is.  Past the cap
        # the picture can only get denser, and the labels shrink with it.
        width = _INCHES_PER_NODE * (max(xs) - min(xs) + 2)
        height = _INCHES_PER_NODE * (max(ys) - min(ys) + 2)
        shrink = min(1, _MAX_FIGSIZE / max(width, height))
        kwargs.setdefault('figsize', (width * shrink, height * shrink))
        kwargs.setdefault('label_fontsize', max(5, 10 * shrink))
        fontsize = kwargs['label_fontsize']

        # The labels are drawn beside the nodes rather than on them: a player
        # just above their node and the payoffs below their leaf.  The player
        # goes on the side away from the branch coming into the node, which is
        # also the side away from the rest of the tree.
        plot = graph.plot(pos=pos, edge_labels=True, vertex_labels=False,
                          vertex_colors=vertex_colors, **kwargs)
        for index, (x, y) in pos.items():
            if not labels[index]:
                continue
            if colors[index] == "white":
                plot += text(labels[index], (x, y - 0.45), color='black',
                             fontsize=fontsize, vertical_alignment='top',
                             zorder=8)
                continue
            parent = parent_of.get(index)
            if parent is not None and pos[parent][0] > x:
                side, alignment = -0.2, 'right'
            else:
                side, alignment = 0.2, 'left'
            plot += text(labels[index], (x + side, y + 0.3), color='black',
                         fontsize=fontsize, horizontal_alignment=alignment,
                         vertical_alignment='bottom', zorder=8)
        # Adding the labels to the graph plot brings the axes back and leaves
        # the outermost of them hanging over the edge of the picture.
        plot.axes(kwargs.get('axes', False))
        plot.set_aspect_ratio(1)
        plot.set_axes_range(min(xs) - 1.5, max(xs) + 1.5,
                            min(ys) - 1.2, max(ys) + 0.8)
        return plot

    def save_efg(self, path):
        r"""
        Save the game to ``path`` in gambit's extensive-form ``.efg`` format.

        The underlying gambit game (see :meth:`_gambit_`) is serialised with
        gambit's writer and written atomically with
        :func:`~sage.misc.temporary_file.atomic_write` so that a partially
        written file is never left behind.  This is the inverse of
        :meth:`load_efg`.

        INPUT:

        - ``path`` -- string; the file path to write the game to

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame(players=['Alice', 'Bob'])
            sage: g.append_move(g.root, 'Alice', ['L', 'R'])
            sage: path = tmp_filename(ext='.efg')
            sage: g.save_efg(path)
            sage: with open(path) as f:
            ....:     print(f.read()[:5])
            EFG 2
        """
        pygambit().require()
        g = self._gambit_()
        with atomic_write(path) as f:   # text mode by default (binary=False)
            f.write(g.to_efg())

    @classmethod
    def load_efg(cls, path):
        r"""
        Read a game from a gambit extensive-form ``.efg`` file.

        The file at ``path`` is read with gambit's ``read_efg`` reader and the
        resulting gambit game is wrapped in a new
        :class:`ExtensiveFormGame` (see :meth:`_gambit_game`).  This is the
        inverse of :meth:`save_efg`.

        INPUT:

        - ``path`` -- string; the path of an ``.efg`` file to read

        OUTPUT: a new :class:`ExtensiveFormGame` wrapping the game in the file

        EXAMPLES:

        A game can be saved and then read back in::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame(players=['Alice', 'Bob'])
            sage: g.append_move(g.root, 'Alice', ['L', 'R'])
            sage: g.set_outcome(g.root.children['L'], 'L', [2, 5])
            sage: g.set_outcome(g.root.children['R'], 'R', [3, 1])
            sage: path = tmp_filename(ext='.efg')
            sage: g.save_efg(path)
            sage: h = ExtensiveFormGame.load_efg(path); h
            An extensive form game with 2 players

        A file that gambit cannot parse is reported as such::

            sage: # optional - pygambit
            sage: path = tmp_filename(ext='.efg')
            sage: with open(path, 'w') as f:
            ....:     _ = f.write('not an efg file')
            sage: ExtensiveFormGame.load_efg(path)
            Traceback (most recent call last):
            ...
            ValueError: Parse error in game file: ...
        """
        pygambit().require()
        return cls(read_efg(path))

    def obtain_nash(self, algorithm=None, use_strategic=False, rational=True,
                    tolerance=1e-4, stop_after='auto'):
        r"""
        Compute the Nash equilibria of the game.

        This delegates to a solver of gambit's ``pygambit.nash`` module (see
        the `pygambit Nash documentation
        <https://gambitproject.readthedocs.io/en/stable/pygambit.api.html#module-pygambit.nash>`_).
        Most of the algorithms operate directly on the extensive form and return
        *behavior-strategy* equilibria; with ``use_strategic`` they work on the
        reduced strategic form instead and return *mixed strategy* equilibria.
        See the OUTPUT section below, as the two are different kinds of object,
        are not indexed in the same way, and which one an algorithm gives is not
        the caller's choice throughout.

        INPUT:

        - ``algorithm`` -- (default: ``None``) the solver to use; one of

          * ``'lcp'`` -- linear complementarity (two-player games), the
            default for games with at most two players

          * ``'lp'`` -- linear programming (two-player *constant-sum* games
            only)

          * ``'enumpure'`` -- enumeration of the pure-strategy equilibria; works
            on the strategic form only

          * ``'enumpoly'`` -- enumeration via systems of polynomial equations
            (any number of players), the default otherwise

          * ``'logit'`` -- the logit quantal response tracing procedure

          * ``'liap'`` -- minimisation of the Lyapunov function, starting from
            the centroid; works on the strategic form only

          * ``'enumpure_agent'`` -- enumeration of the pure-strategy *agent*
            equilibria; works on the extensive form only

          * ``'liap_agent'`` -- minimisation of the *agent* Lyapunov function,
            starting from the centroid of the extensive form; works on the
            extensive form only

          When ``None`` the default is ``'lcp'`` for games with at most two
          players and ``'enumpoly'`` for more.  The name is matched without
          regard to case, so ``'LCP'`` and ``'lcp'`` name the same solver.

          The two ``_agent`` solvers compute a different solution concept.  An
          *agent equilibrium* treats every information set as a player of its
          own -- an "agent" of the player who moves there -- and asks only that
          no agent gain by deviating on its own.  A Nash equilibrium is immune
          to more than that, namely to a player changing what they do at several
          of their information sets at once, so every Nash equilibrium is an
          agent equilibrium but not conversely: a profile can have an agent
          maximum regret of zero and a positive maximum regret.  Gambit's
          `tutorial on the two regrets
          <https://gambitproject.readthedocs.io/en/stable/tutorials/advanced_tutorials/agent_versus_non_agent_regret.html>`_
          works an example out; one is computed below.  Agent equilibria matter
          chiefly as a step towards the refinements -- sequential equilibrium
          and the like -- that ask what a player would do at an information set
          the equilibrium never reaches.  ``'enumpure_agent'`` also works on a
          game of imperfect recall, where a player's pure strategies are not
          something gambit will enumerate, so that ``'enumpure'`` and ``'liap'``
          raise a :class:`RuntimeError` there.

        - ``use_strategic`` -- boolean (default: ``False``); when ``False`` the
          equilibria are computed on the extensive form, when ``True`` on the
          reduced strategic form.  This changes what the method returns, see
          OUTPUT.

          Only ``'lcp'``, ``'lp'``, ``'enumpoly'`` and ``'logit'`` can do both,
          and for them this is passed straight to the gambit solver.  Gambit
          computes pure-strategy and Lyapunov equilibria with one solver per
          form instead of one solver taking an argument, and the two forms do
          not answer the same question there: ``'enumpure'`` and ``'liap'``
          always work on the strategic form and ignore this flag, while
          ``'enumpure_agent'`` and ``'liap_agent'`` always work on the extensive
          form and warn that they ignore it when it is ``True``.

        - ``rational`` -- boolean (default: ``True``); whether to answer with
          rational probabilities rather than floating point ones, which is done
          in whichever of two ways the algorithm allows.  ``'lcp'`` and
          ``'lp'`` are asked to compute exactly throughout; ``'enumpoly'``,
          ``'logit'``, ``'liap'`` and ``'liap_agent'`` have no exact mode, so
          they compute in floating point and their answer is then rounded to the
          exact equilibrium it approximates, as described under ``tolerance``.
          ``'enumpure'`` and ``'enumpure_agent'`` are always exact and ignore
          the flag.

        - ``tolerance`` -- a positive number (default: ``1e-4``); how far a
          probability computed by one of the numerical algorithms may be moved
          in order to round it to a rational.  A rounded profile is returned
          only once gambit has confirmed, in exact arithmetic, that it is an
          equilibrium; when it is not -- as happens when the equilibrium is
          genuinely irrational -- the floating point answer is returned
          instead.  The larger the tolerance the simpler the rationals that are
          tried, so a value that is too small can round an equilibrium to an
          unenlightening exact form rather than fail to round it; the default
          matches the accuracy the gambit solvers themselves aim for.  Has no
          effect when ``rational`` is ``False``, nor on a game whose own
          payoffs -- including the probabilities of its chance moves -- are
          inexact.

        - ``stop_after`` -- (default: ``'auto'``) how many equilibria to compute
          before stopping; ``None`` computes all of them, a positive integer at
          most that many.  Only ``'enumpoly'`` and ``'lcp'`` can stop early, and
          giving this argument with any other algorithm is an error.

          The default, ``'auto'``, is ``1`` for ``'enumpoly'`` and ``None`` for
          ``'lcp'``.  ``'enumpoly'`` works through the supports of the game --
          every choice of which actions are played with positive probability --
          and solves a system of polynomial equations for each one, so its cost
          climbs steeply with the size of the game; as it is also what a game of
          more than two players is solved with by default, it is asked for a
          single equilibrium unless told otherwise.  Pass ``stop_after=None`` to
          have it enumerate them all.

          Gambit lets ``'lcp'`` stop early on the strategic form only, so
          ``use_strategic=True`` has to be passed along with ``stop_after``
          there.

        OUTPUT:

        A list with one gambit profile per computed equilibrium, of one of two
        kinds.

        A ``MixedBehaviorProfile`` is a dict-like object mapping each action at
        each information set to the probability with which that action is
        played, *conditional on that information set being reached*.  Index it
        by an action, an information set or a player, as in ``eq[action]``.
        This is what the extensive form is solved into, and what
        ``use_strategic=False`` (the default) returns, as do
        ``'enumpure_agent'`` and ``'liap_agent'`` whatever it is set to.

        A ``MixedStrategyProfile`` is a dict-like object mapping each *pure
        strategy* -- a complete contingent plan, choosing one action at every
        information set of that player -- to the probability with which the
        plan is played.  Index it by a strategy or a player, as in
        ``eq[strategy]``.  This is what the strategic form is solved into, and
        what ``use_strategic=True`` returns, as do ``'enumpure'`` and
        ``'liap'`` whatever it is set to.  A profile of either kind converts
        into the other with its ``as_behavior()`` and ``as_strategy()``
        methods.

        The two are not interchangeable: indexing a mixed strategy profile by an
        action, or a mixed behavior profile by a strategy, raises a
        :class:`TypeError`.  Either kind of profile comes in a rational and a
        floating point flavour, and which one is returned follows ``rational``
        and ``tolerance`` above.

        This requires the optional gambit package.

        EXAMPLES:

        Alice chooses between ``'L'`` (payoffs ``[2, 5]``) and ``'R'``
        (payoffs ``[3, 1]``); her unique equilibrium is to play ``'R'``::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame(players=['Alice', 'Bob'])
            sage: g.append_move(g.root, 'Alice', ['L', 'R'])
            sage: g.set_outcome(g.root.children['L'], 'L', [2, 5])
            sage: g.set_outcome(g.root.children['R'], 'R', [3, 1])
            sage: eqs = g.obtain_nash()
            sage: [[float(eq[a]) for a in g._gambit_().actions] for eq in eqs]
            [[0.0, 1.0]]

        A different algorithm can be selected explicitly::

            sage: # optional - pygambit
            sage: eqs = g.obtain_nash(algorithm='enumpoly')
            sage: [[float(eq[a]) for a in g._gambit_().actions] for eq in eqs]
            [[0.0, 1.0]]

        With ``use_strategic`` the same equilibrium comes back as a mixed
        strategy profile, which is indexed by each player's strategies rather
        than by the actions of the tree.  Bob never moves here, so he has the
        single trivial strategy::

            sage: # optional - pygambit
            sage: eqs = g.obtain_nash(use_strategic=True)
            sage: [[[float(eq[s]) for s in p.strategies]
            ....:   for p in g._gambit_().players] for eq in eqs]
            [[[0.0, 1.0], [1.0]]]

        ``'lcp'`` and ``'lp'`` compute exactly by default, so the
        probabilities come back as rationals; ``rational=False`` asks for the
        floating point answer instead.  Exact payoffs survive the trip in the
        first place -- a payoff of ``1/3`` reaches gambit as ``1/3``, not as a
        decimal approximation of it::

            sage: # optional - pygambit
            sage: e = ExtensiveFormGame(players=['Alice', 'Bob'])
            sage: e.append_move(e.root, 'Alice', ['L', 'R'])
            sage: e.set_outcome(e.root.children['L'], 'L', [1/3, 5])
            sage: e.set_outcome(e.root.children['R'], 'R', [3, 1])
            sage: e._gambit_().outcomes['L'][e._gambit_().players['Alice']]
            Rational(1, 3)
            sage: type(e.obtain_nash(algorithm='lcp')[0]).__name__
            'MixedBehaviorProfileRational'
            sage: type(e.obtain_nash(algorithm='lcp', rational=False)[0]).__name__
            'MixedBehaviorProfileDouble'

        ``'enumpure'`` and ``'liap'`` solve the strategic form whatever
        ``use_strategic`` says, so they answer with mixed strategy profiles::

            sage: # optional - pygambit
            sage: eqs = g.obtain_nash(algorithm='liap')
            sage: [[[float(eq[s]) for s in p.strategies]
            ....:   for p in g._gambit_().players] for eq in eqs]
            [[[0.0, 1.0], [1.0]]]
            sage: eqs = g.obtain_nash(algorithm='enumpure', use_strategic=True)
            sage: [[[float(eq[s]) for s in p.strategies]
            ....:   for p in g._gambit_().players] for eq in eqs]
            [[[0.0, 1.0], [1.0]]]

        Their two ``_agent`` counterparts solve the extensive form and answer
        with behavior profiles, but not to the same question.  Take Figure 4.2
        of [Mye1991]_, in which Player 1 moves twice without seeing what Player
        2 did in between.  It has one pure Nash equilibrium and two pure agent
        equilibria::

            sage: # optional - pygambit
            sage: myerson = ExtensiveFormGame.load_from_gambit_catalog(
            ....:     'books/myerson1991/fig4_2')
            sage: nash = myerson.obtain_nash(algorithm='enumpure')
            sage: agent = myerson.obtain_nash(algorithm='enumpure_agent')
            sage: len(nash), len(agent)
            (1, 2)
            sage: actions = myerson._gambit_().actions
            sage: [[float(eq[a]) for a in actions] for eq in agent]
            [[1.0, 0.0, 0.0, 1.0, 0.0, 1.0], [0.0, 1.0, 0.0, 1.0, 1.0, 0.0]]

        The first of the two is the Nash equilibrium, written as a behavior
        profile::

            sage: # optional - pygambit
            sage: [float(nash[0].as_behavior()[a]) for a in actions]
            [1.0, 0.0, 0.0, 1.0, 0.0, 1.0]

        The other one is not: Player 1 gains by moving differently at both of
        their information sets, which is a deviation no single agent of theirs
        can make on its own, so gambit's two regrets disagree on it::

            sage: # optional - pygambit
            sage: odd = agent[1]
            sage: QQ(odd.max_regret()), QQ(odd.agent_max_regret())
            (1, 0)

        The agent solvers have no strategic form to work on, so they warn that
        ``use_strategic`` is being ignored and go on to compute the same agent
        equilibria::

            sage: # optional - pygambit
            sage: import warnings
            sage: with warnings.catch_warnings(record=True) as caught:
            ....:     warnings.simplefilter('always')
            ....:     eqs = myerson.obtain_nash(algorithm='enumpure_agent',
            ....:                               use_strategic=True)
            sage: print(caught[0].message)
            'enumpure_agent' computes agent equilibria of the extensive form;
            ignoring use_strategic=True
            sage: len(eqs)
            2

        ``'enumpoly'``, ``'logit'`` and the two ``'liap'`` solvers have no exact
        mode, but their answer is nonetheless exact: it is rounded back to the
        equilibrium it approximates and returned only once gambit has
        confirmed, in exact arithmetic, that the rounded profile is one.  In
        the asymmetric matching pennies game below Bob cannot tell which way
        Alice's coin fell, so his two nodes share an information set, and
        both players mix::

            sage: # optional - pygambit
            sage: pennies = ExtensiveFormGame(players=['Alice', 'Bob'])
            sage: pennies.append_move(pennies.root, 'Alice', ['H', 'T'])
            sage: pennies.append_move(pennies.root.children['H'], 'Bob', ['h', 't'])
            sage: pennies.append_infoset(pennies.root.children['T'],
            ....:                        pennies.root.children['H'])
            sage: for coin, guess, payoff in [('H', 'h', 2), ('H', 't', -1),
            ....:                             ('T', 'h', -1), ('T', 't', 1)]:
            ....:     pennies.set_outcome(pennies.root.children[coin].children[guess],
            ....:                         coin + guess, [payoff, -payoff])
            sage: eqs = pennies.obtain_nash(algorithm='liap_agent')
            sage: type(eqs[0]).__name__
            'MixedBehaviorProfileRational'
            sage: [eqs[0][a] for a in pennies._gambit_().actions]
            [Rational(2, 5), Rational(3, 5), Rational(2, 5), Rational(3, 5)]

        The same holds of the mixed strategy profiles ``use_strategic``
        returns::

            sage: # optional - pygambit
            sage: eq = pennies.obtain_nash(algorithm='logit', use_strategic=True)[0]
            sage: type(eq).__name__
            'MixedStrategyProfileRational'
            sage: [[eq[s] for s in p.strategies] for p in pennies._gambit_().players]
            [[Rational(2, 5), Rational(3, 5)], [Rational(2, 5), Rational(3, 5)]]

        Rounding is never allowed to pass off an approximation as exact, so
        asking for one so fine that nothing verifies gives the floating point
        answer back unchanged::

            sage: # optional - pygambit
            sage: eqs = pennies.obtain_nash(algorithm='liap_agent', tolerance=1e-15)
            sage: type(eqs[0]).__name__
            'MixedBehaviorProfileDouble'
            sage: [float(eqs[0][a]) for a in pennies._gambit_().actions]  # abs tol 1e-6
            [0.4, 0.6, 0.4, 0.6]

        The two representations need not yield the same number of equilibria.
        Below Bob moves only after Alice has played ``'L'``, which she never
        does.  Solving the strategic form leaves his choice at that unreached
        information set unconstrained, so both of his strategies occur in an
        equilibrium; solving the extensive form requires him to act optimally
        there as well, which pins it down::

            sage: # optional - pygambit
            sage: h = ExtensiveFormGame(players=['Alice', 'Bob'])
            sage: h.append_move(h.root, 'Alice', ['L', 'R'])
            sage: h.append_move(h.root.children['L'], 'Bob', ['a', 'b'])
            sage: h.set_outcome(h.root.children['L'].children['a'], 'la', [2, 5])
            sage: h.set_outcome(h.root.children['L'].children['b'], 'lb', [1, 1])
            sage: h.set_outcome(h.root.children['R'], 'R', [3, 1])
            sage: len(h.obtain_nash(algorithm='enumpoly', stop_after=None))
            1
            sage: len(h.obtain_nash(algorithm='enumpoly', use_strategic=True,
            ....:                   stop_after=None))
            2

        ``stop_after`` is what asks for all of them: ``'enumpoly'`` stops at the
        first equilibrium by default, however many the game has::

            sage: # optional - pygambit
            sage: len(h.obtain_nash(algorithm='enumpoly', use_strategic=True))
            1
            sage: len(h.obtain_nash(algorithm='enumpoly', use_strategic=True,
            ....:                   stop_after=2))
            2

        TESTS:

        Indexing a profile by the wrong kind of object is an error::

            sage: # optional - pygambit
            sage: eq = g.obtain_nash(use_strategic=True)[0]
            sage: eq[list(g._gambit_().actions)[0]]
            Traceback (most recent call last):
            ...
            TypeError: profile index must be Player, Strategy, or str, not Action

        The name of the algorithm is not case-sensitive::

            sage: # optional - pygambit
            sage: g.obtain_nash(algorithm='LCP') == g.obtain_nash(algorithm='lcp')
            True

        An unknown name is quoted back the way it was spelled::

            sage: # optional - pygambit
            sage: g.obtain_nash(algorithm='Bogus')
            Traceback (most recent call last):
            ...
            ValueError: unknown algorithm 'Bogus'; must be one of
            'enumpoly', 'enumpure', 'enumpure_agent', 'lcp', 'liap',
            'liap_agent', 'logit', 'lp'

        ::

            sage: # optional - pygambit
            sage: g.obtain_nash(algorithm='bogus')
            Traceback (most recent call last):
            ...
            ValueError: unknown algorithm 'bogus'; must be one of
            'enumpoly', 'enumpure', 'enumpure_agent', 'lcp', 'liap',
            'liap_agent', 'logit', 'lp'

        There is nothing to round to within a tolerance of zero::

            sage: g.obtain_nash(tolerance=0)                     # optional - pygambit
            Traceback (most recent call last):
            ...
            ValueError: 'tolerance' must be positive; got 0

        Only the two solvers that can stop early accept ``stop_after``, and it
        counts equilibria::

            sage: # optional - pygambit
            sage: len(g.obtain_nash(algorithm='lcp', stop_after=1,
            ....:                     use_strategic=True))
            1
            sage: g.obtain_nash(algorithm='lcp', stop_after=1)
            Traceback (most recent call last):
            ...
            ValueError: 'lcp' can only stop early on the strategic form;
            pass use_strategic=True along with 'stop_after'
            sage: g.obtain_nash(algorithm='liap', stop_after=1)
            Traceback (most recent call last):
            ...
            ValueError: 'stop_after' is only supported by the 'enumpoly' and
            'lcp' algorithms; got algorithm 'liap'
            sage: g.obtain_nash(algorithm='enumpoly', stop_after=0)
            Traceback (most recent call last):
            ...
            ValueError: 'stop_after' must be a positive integer or None; got 0

        A game whose payoffs are inexact keeps its inexact equilibria; the
        probabilities of a chance move are payoff data too, so making them
        inexact is enough::

            sage: # optional - pygambit
            sage: coin = ExtensiveFormGame(players=['Alice'])
            sage: coin.append_chance_move(coin.root, ['H', 'T'], probs=['1/3', '2/3'])
            sage: for side in ['H', 'T']:
            ....:     coin.append_move(coin.root.children[side], 'Alice', ['L', 'R'])
            ....:     for move, payoff in [('L', 1), ('R', 2)]:
            ....:         coin.set_outcome(coin.root.children[side].children[move],
            ....:                          side + move, [payoff])
            sage: type(coin.obtain_nash(algorithm='liap_agent')[0]).__name__
            'MixedBehaviorProfileRational'
            sage: coin.set_chance_probs(coin.root, [0.25, 0.75])
            sage: type(coin.obtain_nash(algorithm='liap_agent')[0]).__name__
            'MixedBehaviorProfileDouble'
        """
        # The two game classes round a solver's floating point answer back to
        # an exact equilibrium in the same way.
        from sage.game_theory.normal_form_game import (
            _gambit_payoffs_are_exact, _rationalize_gambit_profile)

        pygambit().require()
        if tolerance <= 0:
            raise ValueError(f"'tolerance' must be positive; got {tolerance!r}")

        game = self._gambit_()
        # Most solvers take ``use_strategic`` as an argument.  ``enumpure`` and
        # ``liap`` come in two flavours instead: the plain ones only ever work
        # on the strategic form, and the agent ones only on the extensive form.
        # ``liap`` is started from a profile rather than from the game, so it
        # needs the centroid of the right kind.  The lambdas read ``stop_after``
        # when they are called, which is after it has been resolved below.
        solvers = {
            'lcp': lambda: gambit_nash.lcp_solve(game, rational=rational,
                                                 use_strategic=use_strategic,
                                                 stop_after=stop_after),
            'lp': lambda: gambit_nash.lp_solve(game, rational=rational,
                                               use_strategic=use_strategic),
            'enumpoly': lambda: gambit_nash.enumpoly_solve(
                game, use_strategic=use_strategic, stop_after=stop_after),
            'logit': lambda: gambit_nash.logit_solve(game,
                                                     use_strategic=use_strategic),
            'enumpure': lambda: gambit_nash.enumpure_solve(game),
            'enumpure_agent': lambda: gambit_nash.enumpure_agent_solve(game),
            'liap': lambda: gambit_nash.liap_solve(game.mixed_strategy_profile()),
            'liap_agent': lambda: gambit_nash.liap_agent_solve(
                game.mixed_behavior_profile()),
        }
        # The algorithm is matched without regard to case, but the error
        # messages below quote the name back the way the caller spelled it.
        requested = algorithm
        if algorithm is None:
            requested = algorithm = ('lcp' if len(game.players) <= 2
                                     else 'enumpoly')
        else:
            algorithm = algorithm.lower()
        try:
            solver = solvers[algorithm]
        except KeyError:
            names = ", ".join(repr(name) for name in sorted(solvers))
            raise ValueError("unknown algorithm {0!r}; must be one of "
                             "{1}".format(requested, names))

        agent = algorithm in _AGENT_ALGORITHMS
        if agent and use_strategic:
            warnings.warn("{0!r} computes agent equilibria of the extensive "
                          "form; ignoring use_strategic=True".format(requested))

        if stop_after == 'auto':
            # Enumerating the supports of a game costs more with every one of
            # them, and ``'enumpoly'`` is what a game of more than two players
            # is solved with by default, so it is stopped at one equilibrium
            # unless the caller asks for more.
            stop_after = 1 if algorithm == 'enumpoly' else None
        else:
            if algorithm not in ('enumpoly', 'lcp'):
                raise ValueError("'stop_after' is only supported by the "
                                 "'enumpoly' and 'lcp' algorithms; got "
                                 "algorithm {0!r}".format(requested))
            if algorithm == 'lcp' and not use_strategic:
                # gambit's own restriction, raised here for a clearer message.
                raise ValueError("'lcp' can only stop early on the strategic "
                                 "form; pass use_strategic=True along with "
                                 "'stop_after'")
            if stop_after is not None:
                if stop_after != int(stop_after) or stop_after < 1:
                    raise ValueError("'stop_after' must be a positive integer "
                                     "or None; got {0!r}".format(stop_after))
                # gambit counts equilibria with a C integer.
                stop_after = int(stop_after)

        equilibria = list(solver().equilibria)
        if rational and _gambit_payoffs_are_exact(game):
            # ``'enumpoly'``, ``'logit'`` and the two ``'liap'`` solvers answer
            # in floating point whatever the game is made of, so ask for the
            # exact equilibrium their answer renders; ``None`` comes back when
            # there is none, and the approximation is then all there is to
            # return.
            equilibria = [_rationalize_gambit_profile(eq, tolerance,
                                                      agent=agent) or eq
                          for eq in equilibria]
        return equilibria

    @classmethod
    def gambit_catalog_games(cls, **kwargs):
        r"""
        List the extensive form games in the gambit catalog.

        The `gambit catalog
        <https://gambitproject.readthedocs.io/en/stable/catalog.html>`_ ships a
        small collection of example games from the literature.  Only the
        extensive form (tree) games are listed here, as those are the ones
        :meth:`load_from_gambit_catalog` can wrap; the rest of the catalog is
        listed by
        :meth:`~sage.game_theory.normal_form_game.NormalFormGame.gambit_catalog_games`.

        INPUT:

        - ``**kwargs`` -- passed on to ``pygambit.catalog.games``, whose
          keywords filter the listing, for instance ``n_players``,
          ``is_const_sum`` or ``include_descriptions``.  ``is_tree`` is set for
          you and passing it has no effect.

        OUTPUT: a :class:`pandas.DataFrame` with a ``Game`` column holding the
        slugs and a ``Title`` column holding the titles.  Slugs are full paths,
        such as ``'journals/geb/bagwell1995'`` or ``'books/myerson1991/fig2_1'``.

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: games = ExtensiveFormGame.gambit_catalog_games()
            sage: 'journals/geb/bagwell1995' in list(games['Game'])
            True
            sage: list(games.columns)
            ['Game', 'Title']

        The keywords of gambit's own catalog listing filter the table::

            sage: # optional - pygambit
            sage: three = ExtensiveFormGame.gambit_catalog_games(n_players=3)
            sage: len(three) < len(games)
            True
            sage: all(ExtensiveFormGame.load_from_gambit_catalog(slug).players.__len__() == 3
            ....:     for slug in three['Game'])
            True

        Every game listed here can be loaded; the strategic form games in the
        catalog, which cannot, are left out::

            sage: # optional - pygambit
            sage: from sage.game_theory.normal_form_game import NormalFormGame
            sage: set(games['Game']) < set(NormalFormGame.gambit_catalog_games()['Game'])
            True
        """
        pygambit().require()
        kwargs['is_tree'] = True
        return catalog.games(**kwargs)

    @classmethod
    def load_from_gambit_catalog(cls, slug):
        r"""
        Load a game from the gambit catalog.

        The available games are listed by :meth:`gambit_catalog_games`.

        INPUT:

        - ``slug`` -- string; the slug of a catalog game, a full path such as
          ``'journals/geb/bagwell1995'``.  The game must be an extensive form
          (tree) game.

        OUTPUT: a new :class:`ExtensiveFormGame` wrapping the catalog game

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame.load_from_gambit_catalog(
            ....:     'journals/geb/bagwell1995')
            sage: g
            An extensive form game with 2 players

        Once loaded, a catalog game behaves like any other game; it can be
        drawn with :meth:`plot` or solved with :meth:`obtain_nash`::

            sage: # optional - pygambit
            sage: len(g.infosets)
            3
            sage: g.is_perfect_recall
            True

        A slug that is not in the catalog is rejected::

            sage: # optional - pygambit
            sage: ExtensiveFormGame.load_from_gambit_catalog('not_a_real_game')
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
