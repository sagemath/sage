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
labels down from the root, ``g.root``.  For instance ``g.root`` is the root and
``g.root.children['L']`` is the child reached by Alice's action ``'L'``.  The
moves themselves take the **action** labels (a required argument).

Its other main purpose is to provide a Sage entry point to gambit's extensive
form games: convert back to the underlying gambit game with
:meth:`~ExtensiveFormGame._gambit_`, save and load games in gambit's ``.efg``
format with :meth:`~ExtensiveFormGame.save_efg` /
:meth:`~ExtensiveFormGame.load_efg`, load example games from the gambit catalog
with :meth:`~ExtensiveFormGame.load_from_gambit_catalog`, and compute Nash
equilibria with :meth:`~ExtensiveFormGame.obtain_nash`.

Game trees can be drawn with :meth:`~ExtensiveFormGame.plot`, either with
Sage's own graph plotting or -- for publication-quality TikZ pictures -- with
the optional `gtdraw <https://www.gambit-project.org/gtdraw/>`_ package.  The
latter gives a :class:`GameTreeTikzPicture`; displaying one compiles it with
LaTeX and opens the picture, so seeing it needs a LaTeX installation on top of
the package itself (see the `gtdraw installation guide
<https://www.gambit-project.org/gtdraw/installation/>`_), while reading its
TikZ source does not.

EXAMPLES:

A two-player game with imperfect information can be built entirely in Sage.
Alice chooses ``'L'`` or ``'R'``; Bob then chooses ``'l'`` or ``'r'`` without
knowing Alice's choice (his two nodes share an information set)::

    sage: # optional - pygambit
    sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
    sage: g = ExtensiveFormGame(players=['Alice', 'Bob'])
    sage: g.append_move(g.root, 'Alice', ['L', 'R'])
    sage: g.append_move(g.root.children['L'], 'Bob', ['l', 'r'])
    sage: g.append_infoset(g.root.children['R'], g.root.children['L'])
    sage: for a, payoff in zip(['l', 'r'], [[3, 1], [0, 0]]):
    ....:     g.set_outcome(g.root.children['L'].children[a], payoff)
    sage: for a, payoff in zip(['l', 'r'], [[0, 0], [1, 3]]):
    ....:     g.set_outcome(g.root.children['R'].children[a], payoff)
    sage: g
    An extensive form game with 2 players
    sage: g.is_perfect_recall
    True

REFERENCES:

- [NN2007]_

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

from sage.structure.sage_object import SageObject
from sage.misc.temporary_file import atomic_write
from sage.misc.latex_standalone import TikzPicture

from sage.features.gambit import pygambit, gtdraw
from sage.misc.lazy_import import lazy_import
lazy_import('pygambit', ['Game', 'read_efg', 'catalog'], feature=pygambit())
lazy_import('pygambit', 'nash', 'gambit_nash', feature=pygambit())
lazy_import('gtdraw', 'tikz', 'gtdraw_tikz', feature=gtdraw())


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
        sage: g.set_outcome(g.root.children['L'], [2, 5])
        sage: g.set_outcome(g.root.children['R'], [3, 1])
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
            sage: g.set_outcome(g.root.children['L'], [2, 5])
            sage: g.set_outcome(g.root.children['R'], [3, 1])
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
            sage: g.set_outcome(g.root.children['L'], [2, 5])
            sage: g.set_outcome(g.root.children['R'], [3, 1])
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
        sage: g.set_outcome(g.root.children['L'], [2, 5])
        sage: g.set_outcome(g.root.children['R'], [3, 1])
        sage: g
        An extensive form game with 2 players
        sage: sorted(p.label for p in g.players)
        ['Alice', 'Bob']

    The familiar textbook games are built the same way.  In the Battle of the
    Sexes, Amy and Bob want to spend the evening together but disagree on how:
    Amy prefers video games and Bob prefers a movie.  They choose
    simultaneously, which in a tree means that Amy is drawn as moving first
    while Bob's two nodes sit in a single information set, so that he cannot
    tell them apart::

        sage: # optional - pygambit
        sage: battle = ExtensiveFormGame(players=['Amy', 'Bob'])
        sage: battle.append_move(battle.root, 'Amy', ['game', 'movie'])
        sage: battle.append_move(battle.root.children['game'], 'Bob',
        ....:                    ['game', 'movie'])
        sage: battle.append_infoset(battle.root.children['movie'],
        ....:                       battle.root.children['game'])
        sage: payoffs = {('game', 'game'): [3, 2], ('game', 'movie'): [1, 1],
        ....:            ('movie', 'game'): [0, 0], ('movie', 'movie'): [2, 3]}
        sage: for (amy, bob), payoff in payoffs.items():
        ....:     battle.set_outcome(battle.root.children[amy].children[bob],
        ....:                        payoff)
        sage: battle
        An extensive form game with 2 players
        sage: len(battle.infosets)
        2

    Its two pure equilibria are the two ways of spending the evening together,
    both playing video games or both watching the movie::

        sage: # optional - pygambit
        sage: actions = battle._gambit_().actions
        sage: [[float(eq[a]) for a in actions]
        ....:  for eq in battle.obtain_nash(algorithm='enumpure')]
        [[1.0, 0.0, 1.0, 0.0], [0.0, 1.0, 0.0, 1.0]]

    The Prisoner's Dilemma has the same shape: two suspects, questioned
    separately, each stay silent or confess without knowing what the other
    does.  Gambit maximizes payoffs, so the sentences are scored as utilities
    (higher is better) rather than as the years in prison of the
    :mod:`~sage.game_theory.normal_form_game` example::

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
        ....:     prisoners.set_outcome(node, payoff)

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
        sage: for leaf, (a, b) in zip(gt.root.children, [[2, 5], [3, 1]]):
        ....:     gt.set_outcome(leaf, gt.add_outcome([a, b]))
        sage: ExtensiveFormGame(gt)
        An extensive form game with 2 players

    Here is matching pennies built that way: Alice and Bob each show a coin,
    Alice taking both pennies if the coins match and Bob if they differ::

        sage: # optional - pygambit
        sage: gt = Game.new_tree(players=['Alice', 'Bob'])
        sage: gt.append_move(gt.root, gt.players['Alice'], ['heads', 'tails'])
        sage: gt.append_move(gt.root.children['heads'], gt.players['Bob'],
        ....:                ['heads', 'tails'])
        sage: gt.append_infoset(gt.root.children['tails'],
        ....:                   gt.root.children['heads'].infoset)
        sage: for alice in ['heads', 'tails']:
        ....:     for bob in ['heads', 'tails']:
        ....:         win = 1 if alice == bob else -1
        ....:         gt.set_outcome(gt.root.children[alice].children[bob],
        ....:                        gt.add_outcome([win, -win]))
        sage: pennies = ExtensiveFormGame(gt); pennies
        An extensive form game with 2 players

    This is a constant-sum game, so it can be solved by linear programming;
    neither player can do better than tossing their coin::

        sage: # optional - pygambit
        sage: actions = pennies._gambit_().actions
        sage: [[float(eq[a]) for a in actions]
        ....:  for eq in pennies.obtain_nash(algorithm='lp')]
        [[0.5, 0.5, 0.5, 0.5]]

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
            sage: g.set_outcome(g.root.children['L'], [2, 5])
            sage: g.set_outcome(g.root.children['R'], [3, 1])
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

    def add_player(self, label=''):
        r"""
        Add a (strategic) player to the game.

        INPUT:

        - ``label`` -- string (default: ``''``); a label identifying the player

        OUTPUT: the label of the new player

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame()
            sage: g.add_player('Alice')
            'Alice'
            sage: [p.label for p in g.players]
            ['Alice']
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
        self._check_actions(actions)
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
        self._check_actions(actions)
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

    def set_outcome(self, node, payoffs):
        r"""
        Set the payoffs awarded at the terminal node ``node``.

        INPUT:

        - ``node`` -- a terminal node (see :meth:`root`)

        - ``payoffs`` -- a list of payoffs, one per player (in the order of
          :meth:`players`), or ``None`` to clear the outcome

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame(players=['Alice', 'Bob'])
            sage: g.append_move(g.root, 'Alice', ['L', 'R'])
            sage: g.set_outcome(g.root.children['L'], [2, 5])
            sage: float(g.root.children['L'].outcome['Alice'])
            2.0
        """
        pygambit().require()
        game = self._gambit_()
        if payoffs is not None:
            payoffs = game.add_outcome(list(payoffs))
        game.set_outcome(node, payoffs)

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
            sage: g.set_outcome(g.root.children['L'], [2, 5])
            sage: g.set_outcome(g.root.children['R'], [3, 1])
            sage: print(g.to_efg().splitlines()[0])
            EFG 2 R "Untitled extensive game" { "Alice" "Bob" }
        """
        pygambit().require()
        return self._gambit_().to_efg()

    def plot(self, backend='sage', **kwargs):
        r"""
        Plot the game tree.

        Two drawing backends are available, selected with ``backend``.

        The default ``'sage'`` backend draws the tree with Sage's own graph
        plotting and needs nothing beyond Sage.  Decision nodes are labeled
        with the moving player, chance nodes with ``'Chance'``, terminal nodes
        with their payoffs, and edges with the action labels; a branch out of a
        chance node also carries its probability, as it does in gtdraw.  Nodes
        belonging to the same information set share a color, so imperfect
        information is visible.  The result is a
        :class:`~sage.plot.graphics.Graphics` object, so it composes with the
        rest of Sage's plotting.

        The ``'gtdraw'`` backend instead hands the underlying gambit game to
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
            ....:     g.set_outcome(g.root.children['L'].children[a], payoff)
            sage: for a, payoff in zip(['l', 'r'], [[0, 0], [1, 3]]):
            ....:     g.set_outcome(g.root.children['R'].children[a], payoff)

        The default backend gives a Sage :class:`~sage.plot.graphics.Graphics`
        object, which can be shown, saved or combined with other graphics::

            sage: # optional - pygambit
            sage: from sage.plot.graphics import Graphics
            sage: isinstance(g.plot(), Graphics)
            True
            sage: isinstance(g.plot(backend='sage'), Graphics)
            True

        Keyword arguments reach Sage's graph plotting::

            sage: # optional - pygambit
            sage: isinstance(g.plot(figsize=8, vertex_size=400), Graphics)
            True

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
            sage: drawn = [p.string for p in coin.plot() if hasattr(p, 'string')]
            sage: sorted(s for s in drawn if '/' in s)
            ['H (1/2)', 'T (1/2)']
            sage: sorted(set(s for s in drawn if s in ['bet', 'fold']))
            ['bet', 'fold']

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
        are full paths, which :meth:`load_from_gambit_catalog` lists::

            sage: # optional - pygambit
            sage: g = ExtensiveFormGame()
            sage: 'books/myerson1991/fig2_1' in list(g.load_from_gambit_catalog()['Game'])
            True

        Myerson's simple poker game is a two-player tree with three information
        sets, and it draws with either backend::

            sage: # optional - pygambit
            sage: g.load_from_gambit_catalog('books/myerson1991/fig2_1', info=False)
            sage: g
            An extensive form game with 2 players
            sage: len(g.infosets)
            3
            sage: isinstance(g.plot(), Graphics)
            True

        ::

            sage: # optional - pygambit gtdraw
            sage: r'\begin{tikzpicture}' in g.plot(backend='gtdraw').content()
            True

        So does the stripped-down poker game of Reiley et al., here drawn
        sideways::

            sage: # optional - pygambit
            sage: h = ExtensiveFormGame()
            sage: h.load_from_gambit_catalog('journals/other/reiley2008/fig1', info=False)
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
            ....:     kuhn.set_outcome(n.children['check'].children['check'], [w, -w])
            ....:     kuhn.set_outcome(n.children['check'].children['bet'].children['fold'],
            ....:                      [-1, 1])
            ....:     kuhn.set_outcome(n.children['check'].children['bet'].children['call'],
            ....:                      [2 * w, -2 * w])
            ....:     kuhn.set_outcome(n.children['bet'].children['fold'], [1, -1])
            ....:     kuhn.set_outcome(n.children['bet'].children['call'], [2 * w, -2 * w])

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
            sage: isinstance(kuhn.plot(), Graphics)
            True

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
        """
        pygambit().require()
        if backend == 'gtdraw':
            gtdraw().require()
            # The preamble mirrors the standalone document gtdraw itself wraps
            # its TikZ code in, so that the picture matches gtdraw's own output.
            return GameTreeTikzPicture(
                gtdraw_tikz(self._gambit_(), **kwargs),
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
        from sage.plot.colors import rainbow

        game = self._gambit_()
        graph = DiGraph()
        labels = {}
        infoset_of = {}
        counter = [0]

        def visit(node):
            index = counter[0]
            counter[0] += 1
            if node.is_terminal:
                if node.outcome is not None:
                    payoffs = ", ".join(str(node.outcome[p]) for p in game.players)
                    labels[index] = "({0})".format(payoffs)
                else:
                    labels[index] = ""
                infoset_of[index] = None
            else:
                infoset = node.infoset
                if infoset.is_chance:
                    labels[index] = "Chance"
                else:
                    labels[index] = (infoset.player.label
                                     or "Player {0}".format(infoset.player.number + 1))
                infoset_of[index] = infoset
                for child, action in zip(node.children, infoset.actions):
                    child_index = visit(child)
                    if infoset.is_chance:
                        # As gtdraw does, a branch out of a chance node carries
                        # its probability next to the action label.
                        prob = "({0})".format(action.prob)
                        edge_label = ("{0} {1}".format(action.label, prob)
                                      if action.label else prob)
                    else:
                        edge_label = action.label
                    graph.add_edge(index, child_index, edge_label)
            return index

        root_index = visit(game.root)

        infosets = [i for i in {id(s): s for s in infoset_of.values()
                                if s is not None and not s.is_chance}.values()]
        vertex_colors = {}
        for color, infoset in zip(rainbow(max(len(infosets), 1)), infosets):
            vertex_colors[color] = [i for i, s in infoset_of.items()
                                    if s is not None and id(s) == id(infoset)]
        vertex_colors["lightgray"] = [i for i, s in infoset_of.items()
                                      if s is not None and s.is_chance]
        vertex_colors["white"] = [i for i, s in infoset_of.items() if s is None]
        vertex_colors = {c: v for c, v in vertex_colors.items() if v}

        return graph.plot(layout='tree', tree_root=root_index, edge_labels=True,
                          vertex_labels=labels, vertex_colors=vertex_colors,
                          **kwargs)

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

    def load_efg(self, path):
        r"""
        Populate this game from a gambit extensive-form ``.efg`` file.

        The file at ``path`` is read with gambit's ``read_efg`` reader and the
        resulting gambit game is wrapped in place (see :meth:`_gambit_game`),
        replacing any game already wrapped.  This is the inverse of
        :meth:`save_efg`.

        INPUT:

        - ``path`` -- string; the path of an ``.efg`` file to read

        EXAMPLES:

        A game can be saved and then read back in::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame(players=['Alice', 'Bob'])
            sage: g.append_move(g.root, 'Alice', ['L', 'R'])
            sage: g.set_outcome(g.root.children['L'], [2, 5])
            sage: g.set_outcome(g.root.children['R'], [3, 1])
            sage: path = tmp_filename(ext='.efg')
            sage: g.save_efg(path)
            sage: h = ExtensiveFormGame()
            sage: h.load_efg(path); h
            An extensive form game with 2 players
        """
        pygambit().require()
        self._gambit_game(read_efg(path))

    def obtain_nash(self, algorithm=None):
        r"""
        Compute the Nash equilibria of the game.

        This delegates to a solver of gambit's ``pygambit.nash`` module (see
        the `pygambit Nash documentation
        <https://gambitproject.readthedocs.io/en/stable/pygambit.api.html#module-pygambit.nash>`_).
        Every supported algorithm operates directly on the extensive form and
        returns *behavior-strategy* equilibria, i.e. a list of gambit
        ``MixedBehaviorProfile`` objects.

        INPUT:

        - ``algorithm`` -- (default: ``None``) the solver to use; one of

          * ``'lcp'`` -- linear complementarity (two-player games), the
            default for games with at most two players

          * ``'lp'`` -- linear programming (two-player *constant-sum* games
            only)

          * ``'enumpure'`` -- enumeration of pure-strategy (agent) equilibria

          * ``'enumpoly'`` -- enumeration via systems of polynomial equations
            (any number of players), the default otherwise

          * ``'logit'`` -- the logit quantal response tracing procedure

          When ``None`` the default is ``'lcp'`` for games with at most two
          players and ``'enumpoly'`` for more.

        OUTPUT: a list of gambit ``MixedBehaviorProfile`` objects, one per
        computed equilibrium.

        This requires the optional gambit package.

        EXAMPLES:

        Alice chooses between ``'L'`` (payoffs ``[2, 5]``) and ``'R'``
        (payoffs ``[3, 1]``); her unique equilibrium is to play ``'R'``::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame(players=['Alice', 'Bob'])
            sage: g.append_move(g.root, 'Alice', ['L', 'R'])
            sage: g.set_outcome(g.root.children['L'], [2, 5])
            sage: g.set_outcome(g.root.children['R'], [3, 1])
            sage: eqs = g.obtain_nash()
            sage: [[float(eq[a]) for a in g._gambit_().actions] for eq in eqs]
            [[0.0, 1.0]]

        A different algorithm can be selected explicitly::

            sage: # optional - pygambit
            sage: eqs = g.obtain_nash(algorithm='enumpoly')
            sage: [[float(eq[a]) for a in g._gambit_().actions] for eq in eqs]
            [[0.0, 1.0]]

        TESTS::

            sage: # optional - pygambit
            sage: g.obtain_nash(algorithm='bogus')
            Traceback (most recent call last):
            ...
            ValueError: unknown algorithm 'bogus'; must be one of
            'enumpoly', 'enumpure', 'lcp', 'logit', 'lp'
        """
        pygambit().require()
        game = self._gambit_()
        solvers = {
            'lcp': gambit_nash.lcp_solve,
            'lp': gambit_nash.lp_solve,
            'enumpure': gambit_nash.enumpure_agent_solve,
            'enumpoly': gambit_nash.enumpoly_solve,
            'logit': gambit_nash.logit_solve,
        }
        if algorithm is None:
            algorithm = 'lcp' if len(game.players) <= 2 else 'enumpoly'
        try:
            solver = solvers[algorithm]
        except KeyError:
            names = ", ".join(repr(name) for name in sorted(solvers))
            raise ValueError("unknown algorithm {0!r}; must be one of "
                             "{1}".format(algorithm, names))
        return list(solver(game).equilibria)

    def load_from_gambit_catalog(self, game=None, info=True):
        r"""
        List extensive form games in the gambit catalog and/or load one.

        The `gambit catalog
        <https://gambitproject.readthedocs.io/en/stable/catalog.html>`_ ships a
        small collection of example games.  Depending on the arguments this
        method lists the available games, wraps one of them in ``self`` (in
        place, replacing any game already wrapped), or both.

        INPUT:

        - ``game`` -- (default: ``None``) the slug of a catalog game to load.
          Slugs are full paths, such as ``'journals/geb/bagwell1995'`` or
          ``'books/myerson1991/fig2_1'``; the ``Game`` column of the table
          returned by this method lists them all.  When ``None`` no game is
          loaded.  The catalog game must be an extensive form (tree) game.

        - ``info`` -- boolean (default: ``True``); when ``True`` return the
          table of available games (a :class:`pandas.DataFrame` with ``Game``
          slugs and ``Title`` columns).

        OUTPUT: the catalog table when ``info`` is ``True``, otherwise ``None``.

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame()
            sage: 'journals/geb/bagwell1995' in list(g.load_from_gambit_catalog()['Game'])
            True
            sage: g.load_from_gambit_catalog('journals/geb/bagwell1995', info=False)
            sage: g
            An extensive form game with 2 players

        Once loaded, a catalog game behaves like any other game; it can be
        drawn with :meth:`plot` or solved with :meth:`obtain_nash`::

            sage: # optional - pygambit
            sage: len(g.infosets)
            3
            sage: g.is_perfect_recall
            True

        A slug that is not in the catalog is rejected (as is a catalog game
        that is not an extensive form game, via :meth:`_gambit_game`)::

            sage: # optional - pygambit
            sage: g.load_from_gambit_catalog('not_a_real_game', info=False)
            Traceback (most recent call last):
            ...
            ValueError: 'not_a_real_game' is not a game in the gambit catalog; ...
        """
        pygambit().require()
        if game is not None:
            try:
                loaded = catalog.load(game)
            except FileNotFoundError:
                raise ValueError(
                    f"{game!r} is not a game in the gambit catalog; call "
                    "load_from_gambit_catalog() with no argument to see the "
                    "available games"
                )
            self._gambit_game(loaded)
        if info:
            return catalog.games()
