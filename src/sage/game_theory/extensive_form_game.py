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
``pygambit`` game.  In the builder, nodes are referred to by string
**labels** so that the user never has to handle ``pygambit`` objects directly:
each move auto-labels the child nodes it creates by their path in the tree
(the root is ``'root'``, its ``'L'`` child is ``'root.L'``, and so on).

Its other main purpose is to provide a Sage entry point to gambit's extensive
form games: convert back to the underlying gambit game with
:meth:`~ExtensiveFormGame._gambit_`, save and load games in gambit's ``.efg``
format with :meth:`~ExtensiveFormGame.save_efg` /
:meth:`~ExtensiveFormGame.load_efg`, and compute Nash equilibria with
:meth:`~ExtensiveFormGame.obtain_nash`.

EXAMPLES:

A two-player game with imperfect information can be built entirely in Sage.
Alice chooses ``'L'`` or ``'R'``; Bob then chooses ``'l'`` or ``'r'`` without
knowing Alice's choice (his two nodes share an information set)::

    sage: # optional - pygambit
    sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
    sage: g = ExtensiveFormGame()
    sage: g.add_player('Alice'); g.add_player('Bob')
    'Alice'
    'Bob'
    sage: g.append_move('root', 'Alice', ['L', 'R'])
    ['root.L', 'root.R']
    sage: g.append_move('root.L', 'Bob', ['l', 'r'])
    ['root.L.l', 'root.L.r']
    sage: g.append_infoset('root.R', 'root.L')
    ['root.R.l', 'root.R.r']
    sage: payoffs = {'root.L.l': [3, 1], 'root.L.r': [0, 0],
    ....:            'root.R.l': [0, 0], 'root.R.r': [1, 3]}
    sage: for leaf, payoff in payoffs.items():
    ....:     g.set_outcome(leaf, payoff)
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

from sage.features.gambit import pygambit
from sage.misc.lazy_import import lazy_import
lazy_import('pygambit', ['Game', 'read_efg', 'catalog'], feature=pygambit())
lazy_import('pygambit', 'nash', 'gambit_nash', feature=pygambit())


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

    In the builder, nodes are referred to by string **labels**.  A freshly
    created child node is auto-labeled by its path in the tree: the root is
    ``'root'``, the child reached by action ``'L'`` is ``'root.L'``, and so on
    (pass ``children`` to override).  The tree-building methods are
    :meth:`add_player`, :meth:`append_move`, :meth:`append_chance_move`,
    :meth:`append_infoset`, :meth:`set_outcome`, :meth:`set_chance_probs`,
    :meth:`insert_move` and :meth:`delete_tree`.

    INPUT:

    - ``generator`` -- the game to wrap; either

      * a ``pygambit`` extensive form (tree) ``Game`` (requires the optional
        gambit package), or

      * ``None`` (default), giving an empty tree ready to build on.

    EXAMPLES:

    Build a two-player game in Sage.  Alice chooses between ``'L'`` and
    ``'R'``, ending the game with payoffs ``[2, 5]`` or ``[3, 1]``::

        sage: # optional - pygambit
        sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
        sage: g = ExtensiveFormGame()
        sage: g.add_player('Alice'); g.add_player('Bob')
        'Alice'
        'Bob'
        sage: g.append_move('root', 'Alice', ['L', 'R'])
        ['root.L', 'root.R']
        sage: g.set_outcome('root.L', [2, 5])
        sage: g.set_outcome('root.R', [3, 1])
        sage: g
        An extensive form game with 2 players
        sage: sorted(p.label for p in g.players)
        ['Alice', 'Bob']

    A game built directly in gambit can also be wrapped::

        sage: # optional - pygambit
        sage: from pygambit import Game
        sage: gt = Game.new_tree(players=['Alice', 'Bob'])
        sage: gt.append_move(gt.root, gt.players['Alice'], ['L', 'R'])
        sage: for leaf, (a, b) in zip(gt.root.children, [[2, 5], [3, 1]]):
        ....:     gt.set_outcome(leaf, gt.add_outcome([a, b]))
        sage: ExtensiveFormGame(gt)
        An extensive form game with 2 players

    REFERENCES:

    - [NN2007]_

    - [Gambit]_
    """

    def __init__(self, generator=None):
        r"""
        Initialize an extensive form game.

        TESTS::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: ExtensiveFormGame()
            An extensive form game with 0 players
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
            self._game = Game.new_tree()
            self._game.root.label = 'root'
        elif Game is not None and isinstance(generator, Game):
            self._gambit_game(generator)
        else:
            raise TypeError("generator must be a gambit extensive form game "
                            "or None")

    def _node(self, label):
        r"""
        Return the unique node of the tree whose label is ``label``.

        This is the internal resolver that lets the tree-building methods refer
        to nodes by string label instead of by ``pygambit`` object.

        TESTS::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame()
            sage: g.add_player('Alice')
            'Alice'
            sage: g.append_move('root', 'Alice', ['L', 'R'])
            ['root.L', 'root.R']
            sage: g._node('root.L').label
            'root.L'
            sage: g._node('nope')
            Traceback (most recent call last):
            ...
            ValueError: no node labeled 'nope' in the tree; current node
            labels are ['root', 'root.L', 'root.R']
        """
        nodes = list(self._gambit_().nodes)
        matches = [n for n in nodes if n.label == label]
        if not matches:
            available = sorted(n.label for n in nodes if n.label)
            raise ValueError("no node labeled {0!r} in the tree; current node "
                             "labels are {1}".format(label, available))
        if len(matches) > 1:
            raise ValueError("more than one node with label {0!r}".format(label))
        return matches[0]

    def _label_children(self, node, node_label, actions, children):
        r"""
        Label the children of ``node`` and return their labels.

        The children are labeled ``children`` if given, otherwise by their path
        ``"<node_label>.<action>"``.  Each new label is checked not to clash
        with a label already used in the tree.

        TESTS::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame()
            sage: g.add_player('Alice')
            'Alice'
            sage: g.append_move('root', 'Alice', ['L', 'R'], children=['x', 'x'])
            Traceback (most recent call last):
            ...
            ValueError: duplicate node label 'x'
        """
        if children is None:
            labels = ["{0}.{1}".format(node_label, a) for a in actions]
        else:
            labels = list(children)
            if len(labels) != len(actions):
                raise ValueError("children must have one label per action")
        existing = {n.label for n in self._gambit_().nodes if n.label}
        seen = set()
        for label in labels:
            if label in existing or label in seen:
                raise ValueError("duplicate node label {0!r}".format(label))
            seen.add(label)
        for child, label in zip(node.children, labels):
            child.label = label
        return labels

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
            sage: g = ExtensiveFormGame()
            sage: g.add_player('Alice')
            'Alice'
            sage: g.append_move('root', 'Alice', ['L', 'R'])
            ['root.L', 'root.R']
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
        The label of the root node, i.e. where to start building the tree.

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: ExtensiveFormGame().root
            'root'
        """
        return self._gambit_().root.label

    @property
    def players(self):
        r"""
        The (strategic) players of the game, as gambit players.

        The chance player is not included (gambit keeps it separate).

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame()
            sage: g.add_player('Alice'); g.add_player('Bob')
            'Alice'
            'Bob'
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
            sage: g = ExtensiveFormGame()
            sage: g.add_player('Alice')
            'Alice'
            sage: g.append_move('root', 'Alice', ['L', 'R'])
            ['root.L', 'root.R']
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
            sage: g = ExtensiveFormGame()
            sage: g.add_player('Alice'); g.add_player('Bob')
            'Alice'
            'Bob'
            sage: g.append_move('root', 'Alice', ['L', 'R'])
            ['root.L', 'root.R']
            sage: g.set_outcome('root.L', [2, 5]); g.set_outcome('root.R', [3, 1])
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
            sage: g = ExtensiveFormGame()
            sage: g.add_player('Alice')
            'Alice'
            sage: g.append_move('root', 'Alice', ['L', 'R'])
            ['root.L', 'root.R']
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

    def append_move(self, node, player, actions, children=None):
        r"""
        Append a move for ``player`` at the terminal node ``node``.

        This makes ``node`` a decision node at which ``player`` chooses among
        ``actions``, creating one child per action.

        INPUT:

        - ``node`` -- string; the label of the (terminal) node to move at

        - ``player`` -- string; the label of the player who moves

        - ``actions`` -- list of strings; the labels of the available actions

        - ``children`` -- (default: ``None``) list of strings; labels for the
          new child nodes.  When ``None`` the children are auto-labeled by
          their path, ``"<node>.<action>"``.

        OUTPUT: the list of labels of the new child nodes

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame()
            sage: g.add_player('Alice'); g.add_player('Bob')
            'Alice'
            'Bob'
            sage: g.append_move('root', 'Alice', ['L', 'R'])
            ['root.L', 'root.R']
            sage: g.append_move('root.L', 'Bob', ['l', 'r'], children=['a', 'b'])
            ['a', 'b']
        """
        pygambit().require()
        parent = self._node(node)
        self._gambit_().append_move(parent, player, actions)
        return self._label_children(parent, node, actions, children)

    def append_chance_move(self, node, actions, probs=None, children=None):
        r"""
        Append a chance (nature) move at the terminal node ``node``.

        This makes ``node`` a chance node with one branch per action; if
        ``probs`` is given it sets the probabilities of those branches.

        INPUT:

        - ``node`` -- string; the label of the (terminal) node to move at

        - ``actions`` -- list of strings; the labels of the chance branches

        - ``probs`` -- (default: ``None``) list of branch probabilities; each
          may be a string such as ``'1/2'`` or a Sage rational

        - ``children`` -- (default: ``None``) list of strings; labels for the
          new child nodes (auto-labeled by path when ``None``)

        OUTPUT: the list of labels of the new child nodes

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame()
            sage: g.add_player('Alice')
            'Alice'
            sage: g.append_chance_move('root', ['H', 'T'], probs=['1/2', '1/2'])
            ['root.H', 'root.T']
            sage: g._node('root').infoset.is_chance
            True
        """
        pygambit().require()
        game = self._gambit_()
        parent = self._node(node)
        game.append_move(parent, game.players.chance, actions)
        if probs is not None:
            game.set_chance_probs(parent.infoset, probs)
        return self._label_children(parent, node, actions, children)

    def append_infoset(self, node, like, children=None):
        r"""
        Append a move at ``node`` in the same information set as ``like``.

        This makes ``node`` a decision node belonging to the information set of
        the node ``like``: the same player moves with the same actions, and the
        two nodes cannot be told apart (imperfect information).

        INPUT:

        - ``node`` -- string; the label of the (terminal) node to move at

        - ``like`` -- string; the label of a node whose information set to join

        - ``children`` -- (default: ``None``) list of strings; labels for the
          new child nodes (auto-labeled by path when ``None``)

        OUTPUT: the list of labels of the new child nodes

        EXAMPLES:

        Bob moves after Alice without observing her choice::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame()
            sage: g.add_player('Alice'); g.add_player('Bob')
            'Alice'
            'Bob'
            sage: g.append_move('root', 'Alice', ['L', 'R'])
            ['root.L', 'root.R']
            sage: g.append_move('root.L', 'Bob', ['l', 'r'])
            ['root.L.l', 'root.L.r']
            sage: g.append_infoset('root.R', 'root.L')
            ['root.R.l', 'root.R.r']
            sage: len(g.infosets)
            2
        """
        pygambit().require()
        parent = self._node(node)
        infoset = self._node(like).infoset
        self._gambit_().append_infoset(parent, infoset)
        actions = [a.label for a in infoset.actions]
        return self._label_children(parent, node, actions, children)

    def set_outcome(self, node, payoffs):
        r"""
        Set the payoffs awarded at the terminal node ``node``.

        INPUT:

        - ``node`` -- string; the label of a terminal node

        - ``payoffs`` -- a list of payoffs, one per player (in the order of
          :meth:`players`), or ``None`` to clear the outcome

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame()
            sage: g.add_player('Alice'); g.add_player('Bob')
            'Alice'
            'Bob'
            sage: g.append_move('root', 'Alice', ['L', 'R'])
            ['root.L', 'root.R']
            sage: g.set_outcome('root.L', [2, 5])
            sage: float(g._node('root.L').outcome['Alice'])
            2.0
        """
        pygambit().require()
        game = self._gambit_()
        if payoffs is not None:
            payoffs = game.add_outcome(list(payoffs))
        game.set_outcome(self._node(node), payoffs)

    def set_chance_probs(self, node, probs):
        r"""
        Set the branch probabilities of the chance move at ``node``.

        INPUT:

        - ``node`` -- string; the label of a node at which a chance move sits
          (see :meth:`append_chance_move`)

        - ``probs`` -- list of branch probabilities; each may be a string such
          as ``'1/3'`` or a Sage rational

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame()
            sage: g.add_player('Alice')
            'Alice'
            sage: g.append_chance_move('root', ['H', 'T'])
            ['root.H', 'root.T']
            sage: g.set_chance_probs('root', ['1/3', '2/3'])
            sage: float(list(g._node('root').infoset.actions)[0].prob)
            0.3333333333333333
        """
        pygambit().require()
        self._gambit_().set_chance_probs(self._node(node).infoset, probs)

    def insert_move(self, node, player, actions):
        r"""
        Insert a new move for ``player`` immediately above ``node``.

        A new decision node with ``actions`` many actions is inserted
        immediately above ``node``; ``node`` (keeping its label and subtree)
        becomes the first child of the new node, which is itself left
        unlabeled.

        INPUT:

        - ``node`` -- string; the label of the node to insert the move above

        - ``player`` -- string; the label of the player who moves

        - ``actions`` -- integer; the number of actions of the new move

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame()
            sage: g.add_player('Alice')
            'Alice'
            sage: g.append_move('root', 'Alice', ['L', 'R'])
            ['root.L', 'root.R']
            sage: g.insert_move('root', 'Alice', 3)
            sage: len(g._gambit_().root.children)   # a new 3-action move above 'root'
            3
            sage: g._node('root').parent == g._gambit_().root
            True
        """
        pygambit().require()
        self._gambit_().insert_move(self._node(node), player, int(actions))

    def delete_tree(self, node):
        r"""
        Delete the subtree below ``node``, making it a terminal node.  The
        node keeps its label (so it can still be referred to, e.g. to give it
        an outcome).

        INPUT:

        - ``node`` -- string; the label of the node to prune

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame()
            sage: g.add_player('Alice')
            'Alice'
            sage: g.append_move('root', 'Alice', ['L', 'R'])
            ['root.L', 'root.R']
            sage: g.delete_tree('root')
            sage: g._node('root').is_terminal
            True
        """
        pygambit().require()
        n = self._node(node)
        self._gambit_().delete_tree(n)
        n.label = node   # gambit clears the label on deletion; restore it

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
            sage: g = ExtensiveFormGame()
            sage: g.add_player('Alice'); g.add_player('Bob')
            'Alice'
            'Bob'
            sage: g.append_move('root', 'Alice', ['L', 'R'])
            ['root.L', 'root.R']
            sage: g.set_outcome('root.L', [2, 5]); g.set_outcome('root.R', [3, 1])
            sage: print(g.to_efg().splitlines()[0])
            EFG 2 R "Untitled extensive game" { "Alice" "Bob" }
        """
        pygambit().require()
        return self._gambit_().to_efg()

    def plot(self, **kwargs):
        r"""
        Plot the game tree.

        This is a Sage-native drawing of the tree (gambit's ``pygambit`` has no
        tree-drawing API): decision nodes are labeled with the moving player,
        chance nodes with ``'Chance'``, terminal nodes with their payoffs, and
        edges with the action labels.  Nodes belonging to the same information
        set share a color, so imperfect information is visible.

        INPUT:

        - ``**kwargs`` -- passed on to the underlying
          :meth:`~sage.graphs.generic_graph.GenericGraph.plot`

        OUTPUT: a :class:`~sage.plot.graphics.Graphics` object

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: from sage.plot.graphics import Graphics
            sage: g = ExtensiveFormGame()
            sage: g.add_player('Alice'); g.add_player('Bob')
            'Alice'
            'Bob'
            sage: g.append_move('root', 'Alice', ['L', 'R'])
            ['root.L', 'root.R']
            sage: g.append_move('root.L', 'Bob', ['l', 'r'])
            ['root.L.l', 'root.L.r']
            sage: g.append_infoset('root.R', 'root.L')
            ['root.R.l', 'root.R.r']
            sage: for leaf, payoff in zip(['root.L.l', 'root.L.r', 'root.R.l',
            ....:         'root.R.r'], [[3, 1], [0, 0], [0, 0], [1, 3]]):
            ....:     g.set_outcome(leaf, payoff)
            sage: isinstance(g.plot(), Graphics)
            True
        """
        pygambit().require()
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
                    graph.add_edge(index, child_index, action.label)
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
            sage: g = ExtensiveFormGame()
            sage: g.add_player('Alice'); g.add_player('Bob')
            'Alice'
            'Bob'
            sage: g.append_move('root', 'Alice', ['L', 'R'])
            ['root.L', 'root.R']
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
            sage: g = ExtensiveFormGame()
            sage: g.add_player('Alice'); g.add_player('Bob')
            'Alice'
            'Bob'
            sage: g.append_move('root', 'Alice', ['L', 'R'])
            ['root.L', 'root.R']
            sage: g.set_outcome('root.L', [2, 5]); g.set_outcome('root.R', [3, 1])
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
            sage: g = ExtensiveFormGame()
            sage: g.add_player('Alice'); g.add_player('Bob')
            'Alice'
            'Bob'
            sage: g.append_move('root', 'Alice', ['L', 'R'])
            ['root.L', 'root.R']
            sage: g.set_outcome('root.L', [2, 5]); g.set_outcome('root.R', [3, 1])
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

        - ``game`` -- (default: ``None``) the slug of a catalog game to load,
          e.g. ``'bagwell1995'``.  When ``None`` no game is loaded.  The
          catalog game must be an extensive form (tree) game.

        - ``info`` -- boolean (default: ``True``); when ``True`` return the
          table of available games (a :class:`pandas.DataFrame` with ``Game``
          slugs and ``Title`` columns).

        OUTPUT: the catalog table when ``info`` is ``True``, otherwise ``None``.

        EXAMPLES::

            sage: # optional - pygambit
            sage: from sage.game_theory.extensive_form_game import ExtensiveFormGame
            sage: g = ExtensiveFormGame()
            sage: 'bagwell1995' in list(g.load_from_gambit_catalog()['Game'])
            True
            sage: g.load_from_gambit_catalog('bagwell1995', info=False)
            sage: g
            An extensive form game with 2 players

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
