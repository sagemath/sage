# sage.doctest: needs sage.graphs
"""
Matching games

This module implements a class for matching games (stable marriage problems)
[DI1989]_. At present the extended Gale-Shapley algorithm is implemented
which can be used to obtain stable matchings.

AUTHORS:

- James Campbell and Vince Knight 06-2014: Original version
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
from sage.structure.sage_object import SageObject
from sage.rings.integer_ring import ZZ
from copy import deepcopy
from sage.graphs.bipartite_graph import BipartiteGraph


class MatchingGame(SageObject):
    r"""
    A matching game.

    A matching game (also called a stable matching problem) models a situation
    in a population of `N` suitors and `N` reviewers. Suitors and reviewers
    rank their preferences and attempt to find a match.

    Formally, a matching game of size `N` is defined by two disjoint sets `S`
    and `R` of size `N`. Associated to each element of `S` and `R` is a
    preference list:

    .. MATH::

        f : S \to R^N
        \text{ and }
        g : R \to S^N.

    Here is an example of matching game on 4 players:

    .. MATH::

        S = \{J, K, L, M\}, \\
        R = \{A, B, C, D\}.

    With preference functions:

    .. MATH::

        f(s) = \begin{cases}
        (A, D, C, B) & \text{ if } s=J,\\
        (A, B, C, D) & \text{ if } s=K,\\
        (B, D, C, A) & \text{ if } s=L,\\
        (C, A, B, D) & \text{ if } s=M,\\
        \end{cases}

        g(s) = \begin{cases}
        (L, J, K, M) & \text{ if } s=A,\\
        (J, M, L, K) & \text{ if } s=B,\\
        (K, M, L, J) & \text{ if } s=C,\\
        (M, K, J, L) & \text{ if } s=D.\\
        \end{cases}

    INPUT:

    Two potential inputs are accepted (see below to see the effect of each):

    - ``reviewer/suitors_preferences`` -- dictionary containing the
      preferences of all players:

      * key - each reviewer/suitors
      * value - a tuple of suitors/reviewers

    OR:

    - ``integer`` -- integer simply representing the number of reviewers
      and suitors

    To implement the above game in Sage::

        sage: suitr_pref = {'J': ('A', 'D', 'C', 'B'),
        ....:               'K': ('A', 'B', 'C', 'D'),
        ....:               'L': ('B', 'D', 'C', 'A'),
        ....:               'M': ('C', 'A', 'B', 'D')}
        sage: reviewr_pref = {'A': ('L', 'J', 'K', 'M'),
        ....:                 'B': ('J', 'M', 'L', 'K'),
        ....:                 'C': ('K', 'M', 'L', 'J'),
        ....:                 'D': ('M', 'K', 'J', 'L')}
        sage: m = MatchingGame([suitr_pref, reviewr_pref])
        sage: m
        A matching game with 4 suitors and 4 reviewers
        sage: m.suitors()
        ('J', 'K', 'L', 'M')
        sage: m.reviewers()
        ('A', 'B', 'C', 'D')

    A matching `M` is any bijection between `S` and `R`. If `s \in S` and
    `r \in R` are matched by `M` we denote:

    .. MATH::

        M(s) = r.

    On any given matching game, one intends to find a matching that is stable.
    In other words, so that no one individual has an incentive to break their
    current match.

    Formally, a stable matching is a matching that has no blocking pairs.
    A blocking pair is any pair `(s, r)` such that `M(s) \neq r` but `s`
    prefers `r` to `M(r)` and `r` prefers `s` to `M^{-1}(r)`.

    To obtain the stable matching in Sage we use the ``solve`` method which
    uses the extended Gale-Shapley algorithm [DI1989]_::

        sage: m.solve()
        {'J': 'A', 'K': 'C', 'L': 'D', 'M': 'B'}

    Matchings have a natural representations as bipartite graphs::

        sage: plot(m)                                                                   # needs sage.plot
        Graphics object consisting of 13 graphics primitives

    The above plots the bipartite graph associated with the matching.
    This plot can be accessed directly::

        sage: graph = m.bipartite_graph()
        sage: graph
        Bipartite graph on 8 vertices

    It is possible to initiate a matching game without having to name each
    suitor and reviewer::

        sage: n = 8
        sage: big_game = MatchingGame(n)
        sage: big_game.suitors()
        (1, 2, 3, 4, 5, 6, 7, 8)
        sage: big_game.reviewers()
        (-1, -2, -3, -4, -5, -6, -7, -8)

    If we attempt to obtain the stable matching for the above game,
    without defining the preference function we obtain an error::

        sage: big_game.solve()
        Traceback (most recent call last):
        ...
        ValueError: suitor preferences are not complete

    To continue we have to populate the preference dictionary. Here
    is one example where the preferences are simply the corresponding
    element of the permutation group::

        sage: from itertools import permutations
        sage: suitr_preferences = list(permutations([-i-1 for i in range(n)]))
        sage: revr_preferences = list(permutations([i+1 for i in range(n)]))
        sage: for player in range(n):
        ....:     big_game.suitors()[player].pref = suitr_preferences[player]
        ....:     big_game.reviewers()[player].pref = revr_preferences[-player]
        sage: big_game.solve()
        {1: -1, 2: -8, 3: -6, 4: -7, 5: -5, 6: -4, 7: -3, 8: -2}

    Note that we can also combine the two ways of creating a game. For example
    here is an initial matching game::

        sage: suitrs = {'Romeo': ('Juliet', 'Rosaline'),
        ....:            'Mercutio': ('Juliet', 'Rosaline')}
        sage: revwrs = {'Juliet': ('Romeo', 'Mercutio'),
        ....:           'Rosaline': ('Mercutio', 'Romeo')}
        sage: g = MatchingGame(suitrs, revwrs)

    Let us assume that all of a sudden a new pair of suitors and reviewers is
    added but their names are not known::

        sage: g.add_reviewer()
        sage: g.add_suitor()
        sage: g.reviewers()
        (-3, 'Juliet', 'Rosaline')
        sage: g.suitors()
        (3, 'Mercutio', 'Romeo')

    Note that when adding a reviewer or a suitor all preferences are wiped::

        sage: [s.pref for s in g.suitors()]
        [[], [], []]
        sage: [r.pref for r in g.reviewers()]
        [[], [], []]

    If we now try to solve the game we will get an error as we have not
    specified the preferences which will need to be updated::

        sage: g.solve()
        Traceback (most recent call last):
        ...
        ValueError: suitor preferences are not complete

    Here we update the preferences so that the new reviewers and suitors
    do not affect things too much (they prefer each other and are the least
    preferred of the others)::

        sage: g.suitors()[1].pref = suitrs['Mercutio'] + (-3,)
        sage: g.suitors()[2].pref = suitrs['Romeo'] + (-3,)
        sage: g.suitors()[0].pref = (-3, 'Juliet', 'Rosaline')
        sage: g.reviewers()[2].pref = revwrs['Rosaline'] + (3,)
        sage: g.reviewers()[1].pref = revwrs['Juliet'] + (3,)
        sage: g.reviewers()[0].pref = (3, 'Romeo', 'Mercutio')

    Now the game can be solved::

        sage: D = g.solve()
        sage: D['Mercutio']
        'Rosaline'
        sage: D['Romeo']
        'Juliet'
        sage: D[3]
        -3

    Note that the above could be equivalently (and more simply) carried out
    by simply updated the original preference dictionaries::

        sage: for key in suitrs:
        ....:     suitrs[key] = suitrs[key] + (-3,)
        sage: for key in revwrs:
        ....:     revwrs[key] = revwrs[key] + (3,)
        sage: suitrs[3] = (-3, 'Juliet', 'Rosaline')
        sage: revwrs[-3] = (3, 'Romeo', 'Mercutio')
        sage: g = MatchingGame(suitrs, revwrs)
        sage: D = g.solve()
        sage: D['Mercutio']
        'Rosaline'
        sage: D['Romeo']
        'Juliet'
        sage: D[3]
        -3

    It can be shown that the Gale-Shapley algorithm will return the stable
    matching that is optimal from the point of view of the suitors and is in
    fact the worst possible matching from the point of view of the reviewers.
    To quickly obtain the matching that is optimal for the reviewers we
    use the ``solve`` method with the ``invert=True`` option::

        sage: left_dict = {'a': ('A', 'B', 'C'),
        ....:              'b': ('B', 'C', 'A'),
        ....:              'c': ('B', 'A', 'C')}
        sage: right_dict = {'A': ('b', 'c', 'a'),
        ....:               'B': ('a', 'c', 'b'),
        ....:               'C': ('a', 'b', 'c')}
        sage: quick_game = MatchingGame([left_dict, right_dict])
        sage: quick_game.solve()
        {'a': 'A', 'b': 'C', 'c': 'B'}
        sage: quick_game.solve(invert=True)
        {'A': 'c', 'B': 'a', 'C': 'b'}

    EXAMPLES:

    8 player letter game::

        sage: suitr_pref = {'J': ('A', 'D', 'C', 'B'),
        ....:               'K': ('A', 'B', 'C', 'D'),
        ....:               'L': ('B', 'D', 'C', 'A'),
        ....:               'M': ('C', 'A', 'B', 'D')}
        sage: reviewr_pref = {'A': ('L', 'J', 'K', 'M'),
        ....:                 'B': ('J', 'M', 'L', 'K'),
        ....:                 'C': ('K', 'M', 'L', 'J'),
        ....:                 'D': ('M', 'K', 'J', 'L')}
        sage: m = MatchingGame([suitr_pref, reviewr_pref])
        sage: m.suitors()
        ('J', 'K', 'L', 'M')
        sage: m.reviewers()
        ('A', 'B', 'C', 'D')

    Also works for numbers::

        sage: suit = {0: (3, 4),
        ....:         1: (3, 4)}
        sage: revr = {3: (0, 1),
        ....:         4: (1, 0)}
        sage: g = MatchingGame([suit, revr])

    Can create a game from an integer. This gives default set of preference
    functions::

        sage: g = MatchingGame(3)
        sage: g
        A matching game with 3 suitors and 3 reviewers

    We have an empty set of preferences for a default named set of
    preferences::

        sage: for s in g.suitors():
        ....:     s, s.pref
        (1, [])
        (2, [])
        (3, [])
        sage: for r in g.reviewers():
        ....:     r, r.pref
        (-1, [])
        (-2, [])
        (-3, [])

    Before trying to solve such a game the algorithm will check if it is
    complete or not::

        sage: g.solve()
        Traceback (most recent call last):
        ...
        ValueError: suitor preferences are not complete

    To be able to obtain the stable matching we must input the preferences::

        sage: for s in g.suitors():
        ....:   s.pref = (-1, -2, -3)
        sage: for r in g.reviewers():
        ....:   r.pref = (1, 2, 3)
        sage: g.solve()
        {1: -1, 2: -2, 3: -3}
    """
    def __init__(self, generator, revr=None):
        r"""
        Initialize a matching game and check the inputs.

        TESTS::

            sage: suit = {0: (3, 4), 1: (3, 4)}
            sage: revr = {3: (0, 1), 4: (1, 0)}
            sage: g = MatchingGame([suit, revr])
            sage: TestSuite(g).run()

            sage: g = MatchingGame(3)
            sage: TestSuite(g).run()

            sage: g2 = MatchingGame(QQ(3))
            sage: g == g2
            True

        The above shows that the input can be either two  dictionaries
        or an integer::

            sage: g = MatchingGame(suit, 3)
            Traceback (most recent call last):
            ...
            TypeError: generator must be an integer or a pair of 2 dictionaries

            sage: g = MatchingGame(matrix(2, [1, 2, 3, 4]))
            Traceback (most recent call last):
            ...
            TypeError: generator must be an integer or a pair of 2 dictionaries

            sage: g = MatchingGame('1,2,3', 'A,B,C')
            Traceback (most recent call last):
            ...
            TypeError: generator must be an integer or a pair of 2 dictionaries
        """
        self._suitors = []
        self._reviewers = []
        if revr is not None:
            generator = [generator, revr]

        if generator in ZZ:
            for i in range(generator):
                self.add_suitor()
                self.add_reviewer()
        elif isinstance(generator[0], dict) and isinstance(generator[1], dict):
            for i in generator[0]:
                self.add_suitor(i)
            for k in generator[1]:
                self.add_reviewer(k)

            for i in self._suitors:
                i.pref = generator[0][i._name]
            for k in self._reviewers:
                k.pref = generator[1][k._name]
        else:
            raise TypeError("generator must be an integer or a pair of 2 dictionaries")

    def _repr_(self):
        r"""
        Return a basic representation of the game stating how many
        players are in the game.

        EXAMPLES:

        Matching game with 2 reviewers and 2 suitors::

            sage: M = MatchingGame(2)
            sage: M
            A matching game with 2 suitors and 2 reviewers
        """
        txt = 'A matching game with {} suitors and {} reviewers'
        return txt.format(len(self._suitors), len(self._reviewers))

    def _latex_(self):
        r"""
        Create the LaTeX representation of the dictionaries for suitors
        and reviewers.

        EXAMPLES::

            sage: suit = {0: (3, 4), 1: (3, 4)}
            sage: revr = {3: (0, 1), 4: (1, 0)}
            sage: g = MatchingGame([suit, revr])
            sage: latex(g)
            \text{Suitors:}
            \begin{aligned}
            \\ 0 & \to (3, 4)
            \\ 1 & \to (3, 4)
            \end{aligned}
            \text{Reviewers:}
            \begin{aligned}
            \\ 3 & \to (0, 1)
            \\ 4 & \to (1, 0)
            \end{aligned}
        """
        output = "\\text{Suitors:}\n\\begin{aligned}"
        for suitor in self._suitors:
            output += "\n\\\\ %s & \\to %s" % (suitor, suitor.pref)
        output += "\n\\end{aligned}\n\\text{Reviewers:}\n\\begin{aligned}"
        for reviewer in self._reviewers:
            output += "\n\\\\ %s & \\to %s" % (reviewer, reviewer.pref)
        return output + "\n\\end{aligned}"

    def __eq__(self, other):
        """
        Check equality.

            sage: suit = {0: (3, 4), 1: (3, 4)}
            sage: revr = {3: (0, 1), 4: (1, 0)}
            sage: g = MatchingGame([suit, revr])
            sage: g2 = MatchingGame([suit, revr])
            sage: g == g2
            True

        Here the two sets of suitors have different preferences::

            sage: suit1 = {0: (3, 4), 1: (3, 4)}
            sage: revr1 = {3: (1, 0), 4: (1, 0)}
            sage: g1 = MatchingGame([suit1, revr1])
            sage: suit2 = {0: (4, 3), 1: (3, 4)}
            sage: revr2 = {3: (1, 0), 4: (1, 0)}
            sage: g2 = MatchingGame([suit2, revr2])
            sage: g == g2
            False

        Here the two sets of reviewers have different preferences::

            sage: suit1 = {0: (3, 4), 1: (3, 4)}
            sage: revr1 = {3: (0, 1), 4: (1, 0)}
            sage: g1 = MatchingGame([suit1, revr1])
            sage: suit2 = {0: (3, 4), 1: (3, 4)}
            sage: revr2 = {3: (1, 0), 4: (0, 1)}
            sage: g2 = MatchingGame([suit2, revr2])
            sage: g == g2
            False

        Note that if two games are created with players ordered differently
        they can still be equal::

            sage: g1 = MatchingGame(1)
            sage: g1.add_reviewer(-2)
            sage: g1.add_reviewer(-3)
            sage: g1.add_suitor(3)
            sage: g1.add_suitor(2)
            sage: g1.reviewers()
            (-1, -2, -3)
            sage: g1.suitors()
            (1, 2, 3)

            sage: g2 = MatchingGame(1)
            sage: g2.add_reviewer(-2)
            sage: g2.add_reviewer(-3)
            sage: g2.add_suitor(2)
            sage: g2.add_suitor(3)
            sage: g2.reviewers()
            (-1, -2, -3)
            sage: g2.suitors()
            (1, 2, 3)

            sage: g1 == g2
            True
        """
        return (isinstance(other, MatchingGame)
                and set(self._suitors) == set(other._suitors)
                and set(self._reviewers) == set(other._reviewers)
                and all(r1.pref == r2.pref for r1, r2 in
                        zip(set(self._reviewers), set(other._reviewers)))
                and all(s1.pref == s2.pref for s1, s2 in
                        zip(set(self._suitors), set(other._suitors))))

    __hash__ = None
   # not hashable because this is mutable.

    def plot(self):
        r"""
        Create the plot representing the stable matching for the game.
        Note that the game must be solved for this to work.

        EXAMPLES:

        An error is returned if the game is not solved::

            sage: suit = {0: (3, 4),
            ....:         1: (3, 4)}
            sage: revr = {3: (0, 1),
            ....:         4: (1, 0)}
            sage: g = MatchingGame([suit, revr])
            sage: plot(g)                                                               # needs sage.plot
            Traceback (most recent call last):
            ...
            ValueError: game has not been solved yet

            sage: g.solve()
            {0: 3, 1: 4}
            sage: plot(g)                                                               # needs sage.plot
            Graphics object consisting of 7 graphics primitives
        """
        pl = self.bipartite_graph()
        return pl.plot()

    def bipartite_graph(self):
        r"""
        Construct a ``BipartiteGraph`` Object of the game.
        This method is similar to the plot method.
        Note that the game must be solved for this to work.

        EXAMPLES:

        An error is returned if the game is not solved::

            sage: suit = {0: (3, 4),
            ....:         1: (3, 4)}
            sage: revr = {3: (0, 1),
            ....:         4: (1, 0)}
            sage: g = MatchingGame([suit, revr])
            sage: g.bipartite_graph()
            Traceback (most recent call last):
            ...
            ValueError: game has not been solved yet

            sage: g.solve()
            {0: 3, 1: 4}
            sage: g.bipartite_graph()
            Bipartite graph on 4 vertices
        """
        self._is_solved()
        graph = BipartiteGraph(self._sol_dict)
        return graph

    def _is_solved(self):
        r"""
        Raise an error if the game has not been solved yet.

        EXAMPLES::

            sage: suit = {0: (3, 4),
            ....:         1: (3, 4)}
            sage: revr = {3: (0, 1),
            ....:         4: (1, 0)}
            sage: g = MatchingGame([suit, revr])
            sage: g._is_solved()
            Traceback (most recent call last):
            ...
            ValueError: game has not been solved yet
            sage: g.solve()
            {0: 3, 1: 4}
            sage: g._is_solved()
        """
        suitor_check = all(s.partner for s in self._suitors)
        reviewer_check = all(r.partner for r in self._reviewers)
        if not suitor_check or not reviewer_check:
            raise ValueError("game has not been solved yet")

    def _is_complete(self):
        r"""
        Raise an error if all players do not have acceptable preferences.

        EXAMPLES:

        Not enough reviewers::

            sage: suit = {0: (3, 4),
            ....:         1: (3, 4)}
            sage: revr = {3: (0, 1)}
            sage: g = MatchingGame([suit, revr])
            sage: g._is_complete()
            Traceback (most recent call last):
            ...
            ValueError: must have the same number of reviewers as suitors

        Not enough suitors::

            sage: suit = {0: (3, 4)}
            sage: revr = {1: (0, 2),
            ....:         3: (0, 1)}
            sage: g = MatchingGame([suit, revr])
            sage: g._is_complete()
            Traceback (most recent call last):
            ...
            ValueError: must have the same number of reviewers as suitors

        Suitors preferences are incomplete::

            sage: suit = {0: (3, 8),
            ....:         1: (0, 0)}
            sage: revr = {3: (0, 1),
            ....:         4: (1, 0)}
            sage: g = MatchingGame([suit, revr])
            sage: g._is_complete()
            Traceback (most recent call last):
            ...
            ValueError: suitor preferences are not complete

        Reviewer preferences are incomplete::

            sage: suit = {0: (3, 4),
            ....:         1: (3, 4)}
            sage: revr = {3: (0, 2, 1),
            ....:         4: (1, 0)}
            sage: g = MatchingGame([suit, revr])
            sage: g._is_complete()
            Traceback (most recent call last):
            ...
            ValueError: reviewer preferences are not complete

        Suitor preferences have repetitions::

            sage: suit = {0: (3,  4),
            ....:         1: (3, 4)}
            sage: revr = {3: (0, 0, 1),
            ....:         4: (1, 0)}
            sage: g = MatchingGame([suit, revr])
            sage: g._is_complete()
            Traceback (most recent call last):
            ...
            ValueError: reviewer preferences contain repetitions

        Reviewer preferences have repetitions::

            sage: suit = {0: (3,  4, 3),
            ....:         1: (3, 4)}
            sage: revr = {3: (0, 1),
            ....:         4: (1, 0)}
            sage: g = MatchingGame([suit, revr])
            sage: g._is_complete()
            Traceback (most recent call last):
            ...
            ValueError: suitor preferences contain repetitions
        """
        if len(self._suitors) != len(self._reviewers):
            raise ValueError("must have the same number of reviewers as suitors")

        for suitor in self._suitors:
            if set(suitor.pref) != set(self._reviewers):
                raise ValueError("suitor preferences are not complete")

        for reviewer in self._reviewers:
            if set(reviewer.pref) != set(self._suitors):
                raise ValueError("reviewer preferences are not complete")

        for reviewer in self._reviewers:
            if len(set(reviewer.pref)) < len(reviewer.pref):
                raise ValueError("reviewer preferences contain repetitions")

        for suitor in self._suitors:
            if len(set(suitor.pref)) < len(suitor.pref):
                raise ValueError("suitor preferences contain repetitions")

    def add_suitor(self, name=None):
        r"""
        Add a suitor to the game.

        INPUT:

        - ``name`` -- can be a string or a number; if left blank will
          automatically generate an integer

        EXAMPLES:

        Creating a two player game::

            sage: g = MatchingGame(2)
            sage: g.suitors()
            (1, 2)

        Adding a suitor without specifying a name::

            sage: g.add_suitor()
            sage: g.suitors()
            (1, 2, 3)

        Adding a suitor while specifying a name::

            sage: g.add_suitor('D')
            sage: g.suitors()
            (1, 2, 3, 'D')

        Note that now our game is no longer complete::

            sage: g._is_complete()
            Traceback (most recent call last):
            ...
            ValueError: must have the same number of reviewers as suitors

        Note that an error is raised if one tries to add a suitor
        with a name that already exists::

            sage: g.add_suitor('D')
            Traceback (most recent call last):
            ...
            ValueError: a suitor with name "D" already exists

        If we add a suitor without passing a name then the name
        of the suitor will not use one that is already chosen::

            sage: suit = {0: (-1,  -2),
            ....:         2: (-2, -1)}
            sage: revr = {-1: (0, 1),
            ....:         -2: (1, 0)}
            sage: g = MatchingGame([suit, revr])
            sage: g.suitors()
            (0, 2)

            sage: g.add_suitor()
            sage: g.suitors()
            (0, 2, 3)
        """
        if name is None:
            name = len(self._suitors) + 1
            while name in self._suitors:
                name += 1
        if any(s._name == name for s in self._suitors):
            raise ValueError('a suitor with name "{}" already exists'.format(name))

        new_suitor = Player(name)
        self._suitors.append(new_suitor)
        for r in self._reviewers:
            r.pref = []

    def add_reviewer(self, name=None):
        r"""
        Add a reviewer to the game.

        INPUT:

        - ``name`` -- can be a string or number; if left blank will
          automatically generate an integer

        EXAMPLES:

        Creating a two player game::

            sage: g = MatchingGame(2)
            sage: g.reviewers()
            (-1, -2)

        Adding a suitor without specifying a name::

            sage: g.add_reviewer()
            sage: g.reviewers()
            (-1, -2, -3)

        Adding a suitor while specifying a name::

            sage: g.add_reviewer(10)
            sage: g.reviewers()
            (-1, -2, -3, 10)

        Note that now our game is no longer complete::

            sage: g._is_complete()
            Traceback (most recent call last):
            ...
            ValueError: must have the same number of reviewers as suitors

        Note that an error is raised if one tries to add a reviewer
        with a name that already exists::

            sage: g.add_reviewer(10)
            Traceback (most recent call last):
            ...
            ValueError: a reviewer with name "10" already exists

        If we add a reviewer without passing a name then the name
        of the reviewer will not use one that is already chosen::

            sage: suit = {0: (-1,  -3),
            ....:         1: (-3, -1)}
            sage: revr = {-1: (0, 1),
            ....:         -3: (1, 0)}
            sage: g = MatchingGame([suit, revr])
            sage: g.reviewers()
            (-1, -3)

            sage: g.add_reviewer()
            sage: g.reviewers()
            (-1, -3, -4)
        """
        if name is None:
            name = -len(self._reviewers) - 1
            while name in self._reviewers:
                name -= 1
        if any(r._name == name for r in self._reviewers):
            raise ValueError('a reviewer with name "{}" already exists'.format(name))

        new_reviewer = Player(name)
        self._reviewers.append(new_reviewer)
        for s in self._suitors:
            s.pref = []

    def suitors(self):
        """
        Return the suitors of ``self``.

        EXAMPLES::

            sage: g = MatchingGame(2)
            sage: g.suitors()
            (1, 2)
        """
        return tuple(sorted(self._suitors, key=lambda s:str(s._name)))

    def reviewers(self):
        """
        Return the reviewers of ``self``.

        EXAMPLES::

            sage: g = MatchingGame(2)
            sage: g.reviewers()
            (-1, -2)
        """
        return tuple(sorted(self._reviewers, key=lambda r:str(r._name)))

    def solve(self, invert=False):
        r"""
        Compute a stable matching for the game using the Gale-Shapley
        algorithm.

        EXAMPLES::

            sage: suitr_pref = {'J': ('A', 'D', 'C', 'B'),
            ....:               'K': ('A', 'B', 'C', 'D'),
            ....:               'L': ('B', 'C', 'D', 'A'),
            ....:               'M': ('C', 'A', 'B', 'D')}
            sage: reviewr_pref = {'A': ('L', 'J', 'K', 'M'),
            ....:                 'B': ('J', 'M', 'L', 'K'),
            ....:                 'C': ('M', 'K', 'L', 'J'),
            ....:                 'D': ('M', 'K', 'J', 'L')}
            sage: m = MatchingGame([suitr_pref, reviewr_pref])
            sage: m.solve()
            {'J': 'A', 'K': 'D', 'L': 'B', 'M': 'C'}

            sage: suitr_pref = {'J': ('A', 'D', 'C', 'B'),
            ....:               'K': ('A', 'B', 'C', 'D'),
            ....:               'L': ('B', 'C', 'D', 'A'),
            ....:               'M': ('C', 'A', 'B', 'D')}
            sage: reviewr_pref = {'A': ('L', 'J', 'K', 'M'),
            ....:                 'B': ('J', 'M', 'L', 'K'),
            ....:                 'C': ('M', 'K', 'L', 'J'),
            ....:                 'D': ('M', 'K', 'J', 'L')}
            sage: m = MatchingGame([suitr_pref, reviewr_pref])
            sage: m.solve(invert=True)
            {'A': 'L', 'B': 'J', 'C': 'M', 'D': 'K'}

            sage: suitr_pref = {1: (-1,)}
            sage: reviewr_pref = {-1: (1,)}
            sage: m = MatchingGame([suitr_pref, reviewr_pref])
            sage: m.solve()
            {1: -1}

            sage: suitr_pref = {}
            sage: reviewr_pref = {}
            sage: m = MatchingGame([suitr_pref, reviewr_pref])
            sage: m.solve()
            {}

        TESTS:

        This also works for players who are both a suitor and reviewer::

            sage: suit = {0: (3,4,2), 1: (3,4,2), 2: (2,3,4)}
            sage: revr = {2: (2,0,1), 3: (0,1,2), 4: (1,0,2)}
            sage: g = MatchingGame(suit, revr)
            sage: g.solve()
            {0: 3, 1: 4, 2: 2}
        """
        self._is_complete()

        for s in self._suitors:
            s.partner = None
        for r in self._reviewers:
            r.partner = None

        if invert:
            reviewers = deepcopy(self._suitors)
            suitors = deepcopy(self._reviewers)
        else:
            suitors = deepcopy(self._suitors)
            reviewers = deepcopy(self._reviewers)

        while any(s.partner is None for s in suitors):
            s = None
            for x in suitors:
                if x.partner is None:
                    s = x
                    break
            r = next((x for x in reviewers if x == s.pref[0]), None)
            if r.partner is None:
                r.partner = s
                s.partner = r
            elif r.pref.index(s._name) < r.pref.index(r.partner._name):
                r.partner.partner = None
                r.partner = s
                s.partner = r
            else:
                s.pref = s.pref[1:]

        if invert:
            suitors, reviewers = reviewers, suitors

        for i, j in zip(self._suitors, suitors):
            i.partner = j.partner
        for i, j in zip(self._reviewers, reviewers):
            i.partner = j.partner

        self._sol_dict = {}
        for s in self._suitors:
            self._sol_dict[s] = [s.partner]
        for r in self._reviewers:
            self._sol_dict[r] = [r.partner]

        if invert:
            return {key: self._sol_dict[key][0] for key in self._reviewers}
        return {key: self._sol_dict[key][0] for key in self._suitors}


class Player:
    r"""
    A class to act as a data holder for the players used of the
    matching games.

    These instances are used when initiating players and to keep track of
    whether or not partners have a preference.
    """
    def __init__(self, name):
        r"""
        TESTS::

            sage: from sage.game_theory.matching_game import Player
            sage: p = Player(10)
            sage: p
            10
            sage: p.pref
            []
            sage: p.partner is None
            True
        """
        self._name = name
        self.pref = []
        self.partner = None

    def __hash__(self):
        r"""
        TESTS::

            sage: from sage.game_theory.matching_game import Player
            sage: p = Player(10)
            sage: d = {p : (1, 2, 3)}
            sage: d
            {10: (1, 2, 3)}
        """
        return hash(self._name)

    def __repr__(self):
        r"""
        TESTS::

            sage: from sage.game_theory.matching_game import Player
            sage: p = Player(10)
            sage: p
            10

            sage: p = Player('Karl')
            sage: p
            'Karl'
        """
        return repr(self._name)

    def __eq__(self, other):
        r"""

        Tests equality of two players. This only checks the name of the player
        and not their preferences.

        TESTS::

            sage: from sage.game_theory.matching_game import Player
            sage: p = Player(10)
            sage: q = Player('Karl')
            sage: p == q
            False

            sage: from sage.game_theory.matching_game import Player
            sage: p = Player(10)
            sage: q = Player(10)
            sage: p == q
            True

            sage: from sage.game_theory.matching_game import Player
            sage: p = Player(10)
            sage: q = Player(10)
            sage: p.pref = (1, 2)
            sage: p.pref = (2, 1)
            sage: p == q
            True
        """
        if isinstance(other, Player):
            return self._name == other._name
        return self._name == other

    def __lt__(self, other):
        """
        Test less than inequality of two players. Allows for players to be
        sorted on their names.

        TESTS::

            sage: from sage.game_theory.matching_game import Player
            sage: p = Player('A')
            sage: q = Player('B')
            sage: p < q
            True
            sage: q < p
            False

            sage: p = Player(0)
            sage: q = Player(1)
            sage: p < q
            True
            sage: q < p
            False
        """
        if isinstance(other, Player):
            return self._name < other._name
        return self._name < other

    def __gt__(self, other):
        """
        Test greater than inequality of two players. Allows for players to be
        sorted on their names.

        TESTS::

            sage: from sage.game_theory.matching_game import Player
            sage: p = Player('A')
            sage: q = Player('B')
            sage: p > q
            False
            sage: q > p
            True

            sage: p = Player(0)
            sage: q = Player(1)
            sage: p > q
            False
            sage: q > p
            True
        """
        if isinstance(other, Player):
            return self._name > other._name
        return self._name > other

    def __ge__(self, other):
        """
        Test greater than or equal inequality of two players. Allows for
        players to be sorted on their names.

        TESTS::

            sage: from sage.game_theory.matching_game import Player
            sage: p = Player('A')
            sage: q = Player('B')
            sage: p >= q
            False
            sage: q >= p
            True

            sage: p = Player(0)
            sage: q = Player(1)
            sage: p >= q
            False
            sage: q >= p
            True

            sage: p = Player(0)
            sage: q = Player(0)
            sage: p >= q
            True

            sage: p = Player('C')
            sage: q = Player('C')
            sage: p >= q
            True
        """
        if isinstance(other, Player):
            return self._name >= other._name
        return self._name >= other

    def __le__(self, other):
        """
        Test less than or equal inequality of two players. Allows for
        players to be sorted on their names.

        TESTS::

            sage: from sage.game_theory.matching_game import Player
            sage: p = Player('A')
            sage: q = Player('B')
            sage: p <= q
            True
            sage: q <= p
            False

            sage: p = Player(0)
            sage: q = Player(1)
            sage: p <= q
            True
            sage: q <= p
            False

            sage: p = Player(0)
            sage: q = Player(0)
            sage: p <= q
            True

            sage: p = Player('C')
            sage: q = Player('C')
            sage: p <= q
            True
        """
        if isinstance(other, Player):
            return self._name <= other._name
        return self._name <= other

    def __ne__(self, other):
        """
        Test inequality of two players. Allows for
        players to be sorted on their names.

        TESTS::

            sage: from sage.game_theory.matching_game import Player
            sage: p = Player('A')
            sage: q = Player('B')
            sage: p != q
            True

            sage: p = Player(0)
            sage: q = Player(1)
            sage: p != q
            True

            sage: p = Player(0)
            sage: q = Player(0)
            sage: p != q
            False

            sage: p = Player('C')
            sage: q = Player('C')
            sage: p != q
            False
        """
        if isinstance(other, Player):
            return self._name != other._name
        return self._name != other
