r"""
Suffix tries and Suffix trees
"""
# ****************************************************************************
#       Copyright (C) 2008 Franco Saliola <saliola@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from itertools import chain

from sage.structure.sage_object import SageObject
from sage.sets.set import Set
from sage.combinat.words.words import Words
from sage.combinat.words.word import Word
from sage.misc.lazy_import import lazy_import
from sage.rings.integer import Integer

lazy_import('sage.graphs.digraph', 'DiGraph')


################################################################################
# Suffix Tries
################################################################################


class SuffixTrie(SageObject):
    def __init__(self, word):
        r"""
        Construct the suffix trie of the word w.

        The suffix trie of a finite word w is a data structure representing
        the factors of w. It is a tree whose edges are labelled with
        letters of w, and whose leafs correspond to suffixes of w.

        This is a straightforward implementation of Algorithm 1 from
        [Ukko1995]_.  It constructs the suffix trie of w[:i] from that
        of w[:i-1].

        A suffix trie is modelled as a deterministic finite-state automaton
        together with the suffix_link map. The set of states corresponds to
        factors of the word (below we write x' for the state corresponding
        to x); these are always 0, 1, .... The state 0 is the initial
        state, and it corresponds to the empty word.  For the purposes of
        the algorithm, there is also an auxiliary state -1. The transition
        function t is defined as::

                t(-1,a) = 0 for all letters a; and
                t(x',a) = y' for all x',y' \in Q such that y = xa,

        and the suffix link function is defined as::

                suffix_link(0) = -1;
                suffix_link(x') = y', if x = ay for some letter a.

        REFERENCES:

        - [Ukko1995]_

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words("cao")("cacao")
            sage: t = SuffixTrie(w); t
            Suffix Trie of the word: cacao

        ::

            sage: e = Words("ab")()
            sage: t = SuffixTrie(e); t
            Suffix Trie of the word:
            sage: t.process_letter("a"); t
            Suffix Trie of the word: a
            sage: t.process_letter("b"); t
            Suffix Trie of the word: ab
            sage: t.process_letter("a"); t
            Suffix Trie of the word: aba

        TESTS::

            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words("cao")("cacao")
            sage: s = SuffixTrie(w)
            sage: loads(dumps(s))
            Suffix Trie of the word: cacao
        """
        # Create the suffix trie for the empty word.
        self._active_state = 0
        self._transition_function = {}
        self._suffix_link = [-1]
        self._alphabet = word.parent().alphabet()

        # Process each letter, in order.
        W = word.parent()
        for letter in word:
            self._process_letter(W([letter]))

    def _process_letter(self, letter):
        r"""
        Process a letter. That is, modify the current suffix trie producing
        the suffix trie for ``self.word() + letter``.

        .. NOTE::

           ``letter`` must occur within the alphabet of the word.

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: t = SuffixTrie(Word("ababba"))
            sage: t._process_letter(Words("ab")("b")); t
            Suffix Trie of the word: ababbab
        """
        r = self._active_state
        old_s = None
        # While r is not the auxiliary vertex, or
        # there is not transition from r along letter, ...
        while r != -1 and (r, letter) not in self._transition_function:
            # adjoin a new state s
            s = len(self._suffix_link)
            self._suffix_link.append(None)
            # create a transition from r to s along letter
            self._transition_function[(r, letter)] = s
            if r != self._active_state:
                # update the suffix link
                self._suffix_link[old_s] = s
            old_s = s
            r = self._suffix_link[r]
        # update the suffix link for the last visited state
        if r == -1:
            self._suffix_link[old_s] = 0
        else:
            self._suffix_link[old_s] = self._transition_function[(r, letter)]
        # update the active state
        self._active_state = \
            self._transition_function[(self._active_state, letter)]

    def process_letter(self, letter):
        r"""
        Modify ``self`` to produce the suffix trie for ``self.word() +
        letter``.

        .. NOTE::

           ``letter`` must occur within the alphabet of the word.

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words("ab")("ababba")
            sage: t = SuffixTrie(w); t
            Suffix Trie of the word: ababba
            sage: t.process_letter("a"); t
            Suffix Trie of the word: ababbaa

        TESTS::

            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words("cao")("cacao")
            sage: t = SuffixTrie(w); t
            Suffix Trie of the word: cacao
            sage: t.process_letter("d")
            Traceback (most recent call last):
            ...
            ValueError: d not in alphabet
        """
        # Make certain that letter is a word containing one letter.
        letter = Words(self._alphabet)([letter])
        self._process_letter(letter)

    #####
    # The following are not necessary for constructing the suffix trie (just
    # the __init__ and process_letter are needed). They just add additional
    # functionality to the class.
    #####

    def _repr_(self) -> str:
        """
        TESTS::

            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: SuffixTrie(Word("abcba"))._repr_()
            'Suffix Trie of the word: abcba'
        """
        return 'Suffix Trie of the %s' % repr(self.word())

    def node_to_word(self, state=0):
        r"""
        Return the word obtained by reading the edge labels from 0 to
        ``state``.

        INPUT:

        - ``state`` -- (default: 0) a state

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words("abc")("abcba")
            sage: t = SuffixTrie(w)
            sage: t.node_to_word(10)
            word: abcba
            sage: t.node_to_word(7)
            word: abcb
        """
        if state == 0:
            return Words(self._alphabet)()
        # We first invert the transition function
        tf_inv = {b: a for a, b in self._transition_function.items()}

        # Starting from the active state,
        # read labels along the unique path to the root.
        u, letter = tf_inv[state]
        w = letter
        s = u
        while s != 0:
            u, letter = tf_inv[s]
            w = letter * w
            s = u
        return w

    def word(self):
        r"""
        Return the word whose suffix tree this is.

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words("abc")("abcba")
            sage: t = SuffixTrie(w)
            sage: t.word()
            word: abcba
            sage: t.word() == w
            True
        """
        return self.node_to_word(self._active_state)

    def __eq__(self, other) -> bool:
        r"""
        If ``self`` and ``other`` have the same transition function,
        the same suffix link, and the same word, then they are equal.

        TESTS::

            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: SuffixTrie(Word("cacao")) == SuffixTrie(Word("ababc"))
            False
            sage: W = Words("cao")
            sage: s = SuffixTrie(W("cacao"))
            sage: t = SuffixTrie(W())
            sage: t.process_letter("c")
            sage: t.process_letter("a")
            sage: t.process_letter("c")
            sage: t.process_letter("a")
            sage: t.process_letter("o")
            sage: t == s
            True
        """
        if not isinstance(other, SuffixTrie):
            return False
        return self._transition_function == other._transition_function \
            and self._suffix_link == other._suffix_link \
            and self.word() == other.word()

    def transition_function(self, node, word):
        r"""
        Return the state reached by beginning at ``node`` and following the
        arrows in the transition graph labelled by the letters of ``word``.

        INPUT:

        - ``node`` -- a node
        - ``word`` -- a word

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words([0,1])([0,1,0,1,1])
            sage: t = SuffixTrie(w)
            sage: all(t.transition_function(u, letter) == v
            ....:     for (u, letter), v in t._transition_function.items())
            True
        """
        if node == -1:
            return self.transition_function(0, word[1:])
        if word.is_empty():
            return 0
        if word.length() == 1:
            return self._transition_function[(node, word)]
        return self.transition_function(
            self._transition_function[(node, word[0:1])], word[1:])

    def states(self):
        r"""
        Return the states of the automaton defined by the suffix trie.

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words([0,1])([0,1,1])
            sage: t = SuffixTrie(w)
            sage: t.states()
            [0, 1, 2, 3, 4]

        ::

            sage: u = Words("aco")("cacao")
            sage: s = SuffixTrie(u)
            sage: s.states()
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
        """
        return list(range(len(self._transition_function)))

    def suffix_link(self, state):
        r"""
        Evaluate the suffix link map of the suffix trie on ``state``.

        Note that the suffix link map is not defined on -1.

        INPUT:

        - ``state`` -- a state

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words("cao")("cacao")
            sage: t = SuffixTrie(w)
            sage: list(map(t.suffix_link, range(13)))
            [-1, 0, 3, 0, 5, 1, 7, 2, 9, 10, 11, 12, 0]
            sage: t.suffix_link(0)
            -1

        TESTS::

            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words("cao")("cacao")
            sage: t = SuffixTrie(w)
            sage: t.suffix_link([1])
            Traceback (most recent call last):
            ...
            TypeError: [1] is not an integer
            sage: t.suffix_link(-1)
            Traceback (most recent call last):
            ...
            TypeError: suffix link is not defined for -1
            sage: t.suffix_link(17)
            Traceback (most recent call last):
            ...
            TypeError: 17 is not a state
        """
        if not isinstance(state, (int, Integer)):
            raise TypeError("%s is not an integer" % state)
        if state == -1:
            raise TypeError("suffix link is not defined for -1")
        if state not in range(len(self._suffix_link)):
            raise TypeError("%s is not a state" % state)
        return self._suffix_link[state]

    def active_state(self):
        r"""
        Return the active state of the suffix trie.

        This is the state corresponding to the word as a suffix of
        itself.

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words("cao")("cacao")
            sage: t = SuffixTrie(w)
            sage: t.active_state()
            8

        ::

            sage: u = Words([0,1])([0,1,1,0,1,0,0,1])
            sage: s = SuffixTrie(u)
            sage: s.active_state()
            22
        """
        return self._active_state

    def final_states(self):
        r"""
        Return the set of final states of the suffix trie.

        These are the states corresponding to the suffixes of
        ``self.word()``. They are obtained be repeatedly following the
        suffix link from the active state until we reach 0.

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words("cao")("cacao")
            sage: t = SuffixTrie(w)
            sage: t.final_states() == Set([8, 9, 10, 11, 12, 0])
            True
        """
        s = self._active_state
        F = [s]
        while s != 0:
            s = self._suffix_link[s]
            F.append(s)
        return Set(F)

    def has_suffix(self, word) -> bool:
        r"""
        Return ``True`` if and only if ``word`` is a suffix of ``self.word()``.

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words("cao")("cacao")
            sage: t = SuffixTrie(w)
            sage: [t.has_suffix(w[i:]) for i in range(w.length()+1)]
            [True, True, True, True, True, True]
            sage: [t.has_suffix(w[:i]) for i in range(w.length()+1)]
            [True, False, False, False, False, True]
        """
        # Find the state corresponding to word, and
        # check to see if s is a final state.
        s = self.transition_function(0, word)
        q = self._active_state
        if q == s:
            return True

        while q != 0:
            q = self._suffix_link[q]
            if q == s:
                return True
        return False

    def to_digraph(self):
        r"""
        Return a ``DiGraph`` object of the transition graph of the suffix
        trie.

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words("cao")("cac")
            sage: t = SuffixTrie(w)
            sage: d = t.to_digraph(); d                                                 # needs sage.graphs
            Digraph on 6 vertices
            sage: d.adjacency_matrix()                                                  # needs sage.graphs sage.modules
            [0 1 0 1 0 0]
            [0 0 1 0 0 0]
            [0 0 0 0 1 0]
            [0 0 0 0 0 1]
            [0 0 0 0 0 0]
            [0 0 0 0 0 0]
        """
        dag = {}
        for (u, letter), v in self._transition_function.items():
            dag.setdefault(u, {})[v] = letter
        return DiGraph(dag)

    def plot(self, layout='tree', tree_root=0, tree_orientation='up',
             vertex_colors=None, edge_labels=True, *args, **kwds):
        r"""
        Return a Graphics object corresponding to the transition graph of
        the suffix trie.

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: SuffixTrie(Word("cacao")).plot()                                      # needs sage.plot
            Graphics object consisting of 38 graphics primitives

        TESTS::

            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: type(SuffixTrie(Word("cacao")).plot())                                # needs sage.plot
            <class 'sage.plot.graphics.Graphics'>
        """
        tree = self.to_digraph()
        for u, v, label in tree.edge_iterator():
            tree.set_edge_label(u, v, label.string_rep())
        if vertex_colors is None:
            suffix_nodes = self.final_states()
            non_suffix_nodes = list(set(self.states()) - set(suffix_nodes))
            vertex_colors = {'#fec7b8': suffix_nodes, '#ffffff': non_suffix_nodes}
        return tree.plot(layout=layout, tree_root=tree_root,
                         tree_orientation=tree_orientation,
                         vertex_colors=vertex_colors, edge_labels=edge_labels,
                         *args, **kwds)

    def show(self, *args, **kwds):
        r"""
        Display the output of :meth:`plot`.

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import SuffixTrie
            sage: w = Words("cao")("cac")
            sage: t = SuffixTrie(w)
            sage: t.show()                                                              # needs sage.plot
        """
        self.plot(*args, **kwds).show()

################################################################################
# Suffix Trees
################################################################################


class ImplicitSuffixTree(SageObject):
    def __init__(self, word):
        r"""
        Construct the implicit suffix tree of a word w.

        The suffix tree of a word w is a compactification of the suffix
        trie for w. The compactification removes all nodes that have
        exactly one incoming edge and exactly one outgoing edge. It
        consists of two components: a tree and a word. Thus, instead of
        labelling the edges by factors of w, we can labelled them by
        indices of the occurrence of the factors in w.

        The following is a straightforward implementation of Ukkonen's
        on-line algorithm for constructing the
        implicit suffix tree [Ukko1995]_.  It constructs the suffix tree for
        w[:i] from that of w[:i-1].

        GENERAL IDEA. The suffix tree of w[:i+1] can be obtained from that
        of w[:i] by visiting each node corresponding to a suffix of w[:i]
        and modifying the tree by applying one of two rules (either append
        a new node to the tree, or split an edge into two). The "active
        state" is the node where the algorithm begins and the "suffix link"
        carries us to the next node that needs to be dealt with.

        TREE. The tree is modelled as an automaton, which is stored as a
        dictionary of dictionaries: it is keyed by the nodes of the tree,
        and the corresponding dictionary is keyed by pairs `(i,j)` of
        integers representing the word w[i-1:j]. This makes it faster to
        look up a particular transition beginning at a specific node.

        STATES/NODES. The states will always be -1, 0, 1, ..., n. The state
        -1 is special and is only used for the purposes of the algorithm.
        All transitions map -1 to 0, so this information is not explicitly
        stored in the transition function.

        EXPLICIT/IMPLICIT NODES. By definition, some of the nodes will not
        be states, but merely locations along an edge; these are called
        implicit nodes. A node r (implicit or explicit) is referenced as a
        pair (s,(k,p)) where s is an ancestor of r and w[k-1:p] is the word
        read by transitioning from s to r in the suffix trie. A reference
        pair is canonical if s is the closest ancestor of r.

        SUFFIX LINK. The algorithm makes use of a map from (some) nodes to
        other nodes, called the suffix link. This is stored as a
        dictionary.

        ACTIVE STATE. We store as ._active_state the active state of the
        tree, the state where the algorithm will begin when processing the
        next letter.

        RUNNING TIME. The running time and storage space of the algorithm
        is linear in the length of the word w (whereas for a suffix tree it
        is quadratic).

        REFERENCES:

        - [Ukko1995]_

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: w = Words("aco")("cacao")
            sage: t = ImplicitSuffixTree(w); t
            Implicit Suffix Tree of the word: cacao
            sage: ababb = Words([0,1])([0,1,0,1,1])
            sage: s = ImplicitSuffixTree(ababb); s
            Implicit Suffix Tree of the word: 01011

        TESTS::

            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: w = Words("cao")("cacao")
            sage: s = ImplicitSuffixTree(w)
            sage: loads(dumps(s))
            Implicit Suffix Tree of the word: cacao
        """
        # For constructing the suffix tree.
        self._transition_function = {0: {}}
        self._suffix_link = {0: -1}
        self._active_state = (0, (1, 1))
        self._letters = []
        for letter in word:
            self._letters.append(letter)
            self._process_letter(letter)
        # _word is not needed for constructing the suffix tree,
        # but it is useful for the other methods.
        self._word = word

    def _process_letter(self, letter):
        r"""
        This is the main part of Ukkonen's algorithm.

        This corresponds to the algorithm "update" in [Ukko1995]_.

        .. NOTE::

           This function is a helper and does not update ``self._data`` and
           ``self._word``.

        REFERENCES:

        - [Ukko1995]_

        TESTS::

            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: w = Words("aco")("caca")
            sage: t = ImplicitSuffixTree(w); t
            Implicit Suffix Tree of the word: caca
            sage: new_letter = "o"
            sage: t._letters.append("o")
            sage: t._process_letter("o")
            sage: t._word = Words("aco")("cacao")
            sage: t
            Implicit Suffix Tree of the word: cacao

        ::

            sage: W = Words([0,1])
            sage: s = ImplicitSuffixTree(W([0,1,0,1])); s
            Implicit Suffix Tree of the word: 0101
            sage: s._letters.append(1)
            sage: s._process_letter(1)
            sage: s._word = W([0,1,0,1,1])
            sage: s
            Implicit Suffix Tree of the word: 01011
        """
        s, (k, i) = self._active_state
        old_r = 0
        end_state, r = self._test_and_split(s, (k, i-1), letter)
        while not end_state:
            # adjoin a new state rr and create a transition from r to rr
            rr = len(self._transition_function)
            self._transition_function[rr] = {}
            self._transition_function[r][(i, None)] = rr
            # update the suffix link, if necessary
            if old_r != 0:
                self._suffix_link[old_r] = r
            old_r = r
            # follow the suffix link to the next state
            s, k = self._canonize(self._suffix_link[s], (k, i-1))
            end_state, r = self._test_and_split(s, (k, i-1), letter)
        # update the suffix link, if necessary
        if old_r != 0:
            self._suffix_link[old_r] = s
        # set the active state
        s, k = self._canonize(s, (k, i))
        self._active_state = (s, (k, i+1))

    def _test_and_split(self, s, k_p, letter):
        r"""
        Helper function for :meth:`_process_letter`. Test to see
        whether an edge needs to be split.

        OUTPUT:

        ``(True, state)``, where ``state`` is the next state to
        process (either a newly created state or the original ``s``).

        TESTS::

            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: w = Words("aco")("caca")
            sage: t = ImplicitSuffixTree(w)
            sage: t._letters.append(w.parent().alphabet().rank("o"))
            sage: t._test_and_split(0, (4,5), w.parent().alphabet().rank("o"))
            (False, 3)
        """
        k, p = k_p
        if k <= p:
            # find the transition from s that begins with k-th letter
            (kk, pp), ss = self._find_transition(s, self._letters[k - 1])
            if letter == self._letters[kk + p - k]:
                return (True, s)
            else:
                # replace transition above by transitions
                del self._transition_function[s][(kk, pp)]
                r = len(self._transition_function)
                self._transition_function[r] = {}
                self._transition_function[s][(kk, kk+p-k)] = r
                self._transition_function[r][(kk+p-k+1, pp)] = ss
                return (False, r)
        else:
            transition = self._find_transition(s, letter)
            if transition is None:
                return (False, s)
            else:
                return (True, s)

    def _canonize(self, s, k_p):
        r"""
        Given an implicit or explicit reference pair for a node, return
        the canonical reference pair.

        Recall that a node r is referenced as (s, (k,p)), where s is an
        ancestor or r and w[k-1:p] is the word obtained by reading the edge
        labels along the path from s to r. A reference pair is canonical if
        s is the closest ancestor of r.

        TESTS::

            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: t = ImplicitSuffixTree(Word("cacao"))
            sage: t._canonize(0,(3,5))
            (3, 5)
            sage: t._canonize(0,(2,5))
            (5, 3)
        """
        k, p = k_p
        if p < k:
            return (s, k)
        else:
            (kk, pp), ss = self._find_transition(s, self._letters[k - 1])
            while pp is not None and pp - kk <= p - k:
                k = k + pp - kk + 1
                s = ss
                if k <= p:
                    (kk, pp), ss = self._find_transition(s, self._letters[k-1])
            return (s, k)

    def _find_transition(self, state, letter):
        r"""
        Return the transition from state that begins with letter.

        This returns ``None`` if no such transition exists.

        The transitions are stored as a dictionary of dictionaries: keyed
        by the nodes, with the corresponding dictionary keyed by pairs
        `(i,j)` of integers representing the word w[i-1:j].

        ._transition_function = {..., node: {(i,j): target_node, ...} }

        TESTS::

            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: t = ImplicitSuffixTree(Word("cacao"))
            sage: t._find_transition(-1, "c")
            ((0, 0), 0)
            sage: t._find_transition(0, "a")
            ((2, 2), 5)
            sage: t._find_transition(0, "c")
            ((1, 2), 3)
            sage: t._find_transition(5, "c")
            ((3, None), 2)
            sage: t._find_transition(5, "a")

        ::

            sage: t = ImplicitSuffixTree(Word([0,1,0,1,1]))
            sage: t._find_transition(3, 1)
            ((5, None), 4)
        """
        if state == -1:
            return ((0, 0), 0)

        if state in self._transition_function:
            for (k, p), s in self._transition_function[state].items():
                if self._letters[k - 1] == letter:
                    return ((k, p), s)
        return None

    #####
    # The following are not necessary for constructing the implicit suffix
    # tree; they add additional functionality to the class.
    #####

    #####
    # Visualization
    #####

    def _repr_(self) -> str:
        r"""
        TESTS::

            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: ImplicitSuffixTree(Word("abcba"))._repr_()
            'Implicit Suffix Tree of the word: abcba'
        """
        return 'Implicit Suffix Tree of the %s' % repr(self.word())

    def word(self):
        r"""
        Return the word whose implicit suffix tree this is.

        TESTS::

            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: ImplicitSuffixTree(Word([0,1,0,1,0])).word() == Word([0,1,0,1,0])
            True
        """
        return self._word

    def transition_function_dictionary(self) -> dict:
        r"""
        Return the transition function as a dictionary of dictionaries.

        The format is consistent with the input format for ``DiGraph``.

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: W = Words("aco")
            sage: t = ImplicitSuffixTree(W("cac"))
            sage: t.transition_function_dictionary()
            {0: {1: (0, None), 2: (1, None)}}

        ::

            sage: W = Words([0,1])
            sage: t = ImplicitSuffixTree(W([0,1,0]))
            sage: t.transition_function_dictionary()
            {0: {1: (0, None), 2: (1, None)}}
        """
        d = {}
        for u, v, (i, j) in self.edge_iterator():
            d.setdefault(u, {})[v] = (i, j)
        return d

    def to_digraph(self, word_labels=False):
        r"""
        Return a ``DiGraph`` object of the transition graph of the suffix tree.

        INPUT:

        - ``word_labels`` -- boolean (default: ``False``); if ``False``, labels
          the edges by pairs `(i, j)`. If ``True``, labels the edges by ``word[i:j]``

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: W = Words([0,1,2])
            sage: t = ImplicitSuffixTree(W([0,1,0,1,2]))
            sage: t.to_digraph()                                                        # needs sage.graphs
            Digraph on 8 vertices
        """
        if not self._letters:
            d = {0: {}}
            return DiGraph(d)
        d = self.transition_function_dictionary()
        for u in d:
            for v, (i, j) in d[u].items():
                if word_labels:
                    d[u][v] = self._word[i:j]
                elif j is None:
                    d[u][v] = (i, len(self._letters))
        return DiGraph(d)

    def plot(self, word_labels=False, layout='tree', tree_root=0,
             tree_orientation='up', vertex_colors=None, edge_labels=True,
             *args, **kwds):
        r"""
        Return a Graphics object corresponding to the transition graph of
        the suffix tree.

        INPUT:

        - ``word_labels`` -- boolean (default: ``False``); if ``False``, labels
          the edges by pairs `(i, j)`; if ``True``, labels the edges by ``word[i:j]``
        - ``layout`` -- (default: ``'tree'``)
        - ``tree_root`` -- (default: 0)
        - ``tree_orientation`` -- (default: ``'up'``)
        - ``vertex_colors`` -- (default: ``None``)
        - ``edge_labels`` -- (default: ``True``)

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: ImplicitSuffixTree(Word('cacao')).plot(word_labels=True)              # needs sage.graphs sage.plot
            Graphics object consisting of 23 graphics primitives
            sage: ImplicitSuffixTree(Word('cacao')).plot(word_labels=False)             # needs sage.graphs sage.plot
            Graphics object consisting of 23 graphics primitives

        TESTS::

            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: type(ImplicitSuffixTree(Word('cacao')).plot(word_labels=True))        # needs sage.graphs sage.plot
            <class 'sage.plot.graphics.Graphics'>
            sage: type(ImplicitSuffixTree(Word('cacao')).plot(word_labels=False))       # needs sage.graphs sage.plot
            <class 'sage.plot.graphics.Graphics'>
        """
        tree = self.to_digraph(word_labels=word_labels)
        if word_labels:
            for u, v, label in tree.edge_iterator():
                tree.set_edge_label(u, v, label.string_rep())
        if vertex_colors is None:
            vertex_colors = {'#fec7b8': tree.vertices(sort=True)}
        return tree.plot(layout=layout, tree_root=tree_root,
                         tree_orientation=tree_orientation,
                         vertex_colors=vertex_colors, edge_labels=edge_labels,
                         *args, **kwds)

    def show(self, word_labels=None, *args, **kwds):
        r"""
        Display the output of :meth:`plot`.

        INPUT:

        - ``word_labels`` -- (default: ``None``) if ``False``, labels the edges
          by pairs `(i, j)`; if ``True``, labels the edges by ``word[i:j]``

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: w = Words("cao")("cacao")
            sage: t = ImplicitSuffixTree(w)
            sage: t.show(word_labels=True)                                              # needs sage.plot
            sage: t.show(word_labels=False)                                             # needs sage.plot
        """
        self.plot(word_labels=word_labels, *args, **kwds).show()

    #####
    # Various methods
    #####

    def __eq__(self, other) -> bool:
        r"""
        If ``self`` and ``other`` have the same transition function
        and the same word, then they are equal.

        TESTS::

            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: w = Words([0,1,2])([0,1,0,1,2])
            sage: u = Words([0,1,2])(iter([0,1,0,1,2]))[:5]
            sage: ImplicitSuffixTree(w) == ImplicitSuffixTree(u)
            True
        """
        if not isinstance(other, ImplicitSuffixTree):
            return False
        return self._transition_function == other._transition_function \
            and self._letters == other._letters

    def transition_function(self, word, node=0):
        r"""
        Return the node obtained by starting from ``node`` and following the
        edges labelled by the letters of ``word``.

        OUTPUT:

        ``("explicit", end_node)`` if we end at ``end_node``, or
        ``("implicit", edge, d)`` if we end `d` spots along an edge.

        INPUT:

        - ``word`` -- a word
        - ``node`` -- (default: 0) starting node

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: W = Words([0,1,2])
            sage: t = ImplicitSuffixTree(W([0,1,0,1,2]))
            sage: t.transition_function(W([0,1,0]))
            ('implicit', (3, 1), 1)
            sage: t.transition_function(W([0,1,2]))
            ('explicit', 4)
            sage: t.transition_function(W([0,1,2]), 5)
            ('explicit', 2)
            sage: t.transition_function(W([0,1]), 5)
            ('implicit', (5, 2), 2)
        """
        if word.is_empty():
            return "explicit", node
        (k, p), s = self._find_transition(node, word[0])
        if p is None:
            # test that word is a prefix of self._letters[k-1:]
            if word == self._word[k-1:(k-1)+word.length()]:
                if word.length() == len(self._letters) - k + 1:
                    return "explicit", s
                else:
                    edge = (node, s)
                    return "implicit", edge, word.length()
        else:
            # find longest common prefix
            m = min(p-k+1, word.length())
            i = 0
            while i < m and self._word[k-1+i] == word[i]:
                i += 1
            if i == p-k+1:
                return self.transition_function(word[p-k+1:], s)
            else:
                edge = (node, s)
                return "implicit", edge, i
            return "explicit", node

    def states(self) -> list:
        r"""
        Return the states (explicit nodes) of the suffix tree.

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: W = Words([0,1,2])
            sage: t = ImplicitSuffixTree(W([0,1,0,1,2]))
            sage: t.states()
            [0, 1, 2, 3, 4, 5, 6, 7]
        """
        return list(range(len(self._transition_function)))

    def suffix_link(self, state):
        r"""
        Evaluate the suffix link map of the implicit suffix tree on ``state``.

        Note that the suffix link is not defined for all states.

        The suffix link of a state `x'` that corresponds to the suffix `x` is
        defined to be -1 is `x'` is the root (0) and `y'` otherwise, where `y'`
        is the state corresponding to the suffix ``x[1:]``.

        INPUT:

        - ``state`` -- a state

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: W = Words([0,1,2])
            sage: t = ImplicitSuffixTree(W([0,1,0,1,2]))
            sage: t.suffix_link(3)
            5
            sage: t.suffix_link(5)
            0
            sage: t.suffix_link(0)
            -1
            sage: t.suffix_link(-1)
            Traceback (most recent call last):
                ...
            TypeError: there is no suffix link from -1
        """
        if state in self._suffix_link:
            return self._suffix_link[state]
        raise TypeError("there is no suffix link from %s" % state)

    def active_state(self):
        r"""
        Return the active state of the suffix tree.

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: W = Words([0,1,2])
            sage: t = ImplicitSuffixTree(W([0,1,0,1,2]))
            sage: t.active_state()
            (0, (6, 6))
        """
        return self._active_state

    def process_letter(self, letter):
        r"""
        Modify the current implicit suffix tree producing the implicit
        suffix tree for ``self.word() + letter``.

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: w = Words("aco")("cacao")
            sage: t = ImplicitSuffixTree(w[:-1]); t
            Implicit Suffix Tree of the word: caca
            sage: t.process_letter(w[-1]); t
            Implicit Suffix Tree of the word: cacao

        ::

            sage: W = Words([0,1])
            sage: s = ImplicitSuffixTree(W([0,1,0,1])); s
            Implicit Suffix Tree of the word: 0101
            sage: s.process_letter(W([1])[0]); s
            Implicit Suffix Tree of the word: 01011
        """
        self._word = self._word * self._word._parent([letter])
        self._letters.append(letter)
        self._process_letter(letter)

    def to_explicit_suffix_tree(self):
        r"""
        Convert ``self`` to an explicit suffix tree.

        It is obtained by processing an end of string letter as if it
        were a regular letter, except that no new leaf nodes are
        created (thus, the only thing that happens is that some
        implicit nodes become explicit).

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: w = Words("aco")("cacao")
            sage: t = ImplicitSuffixTree(w)
            sage: t.to_explicit_suffix_tree()

        ::

            sage: W = Words([0,1])
            sage: s = ImplicitSuffixTree(W([0,1,0,1,1]))
            sage: s.to_explicit_suffix_tree()
        """
        # append a new unique symbol to the word and process the new letter
        end_of_string = object()
        self._letters.append(end_of_string)
        s, (k, i) = self._active_state
        end_state, r = self._test_and_split(s, (k, i-1), end_of_string)
        while not end_state:
            s, k = self._canonize(self._suffix_link[s], (k, i-1))
            end_state, r = self._test_and_split(s, (k, i-1), end_of_string)
        # remove the end of string symbol from the word
        self._letters.pop()

    def edge_iterator(self):
        r"""
        Return an iterator over the edges of the suffix tree.

        The edge from `u` to `v` labelled by `(i,j)` is yielded as
        the tuple `(u,v,(i,j))`.

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: sorted( ImplicitSuffixTree(Word("aaaaa")).edge_iterator() )
            [(0, 1, (0, None))]
            sage: sorted( ImplicitSuffixTree(Word([0,1,0,1])).edge_iterator() )
            [(0, 1, (0, None)), (0, 2, (1, None))]
            sage: sorted( ImplicitSuffixTree(Word()).edge_iterator() )
            []
        """
        queue = [0]
        while queue:
            v = queue.pop()
            for (i, j), u in self._transition_function[v].items():
                yield (v, u, (i - 1, j))
                queue.append(u)

    def number_of_factors(self, n=None):
        r"""
        Count the number of distinct factors of ``self.word()``.

        INPUT:

        - ``n`` -- integer or ``None``

        OUTPUT:

        If ``n`` is an integer, returns the number of distinct factors
        of length ``n``. If ``n`` is ``None``, returns the total number of
        distinct factors.

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: t = ImplicitSuffixTree(Word([1,2,1,3,1,2,1]))
            sage: t.number_of_factors()
            22
            sage: t.number_of_factors(1)
            3
            sage: t.number_of_factors(9)
            0
            sage: t.number_of_factors(0)
            1

        ::

            sage: t = ImplicitSuffixTree(Word("cacao"))
            sage: t.number_of_factors()
            13
            sage: list(map(t.number_of_factors, range(10)))
            [1, 3, 3, 3, 2, 1, 0, 0, 0, 0]

        ::

            sage: t = ImplicitSuffixTree(Word("c"*1000))
            sage: t.number_of_factors()
            1001
            sage: t.number_of_factors(17)
            1
            sage: t.number_of_factors(0)
            1

        ::

            sage: ImplicitSuffixTree(Word()).number_of_factors()
            1

        ::

            sage: blueberry = ImplicitSuffixTree(Word("blueberry"))
            sage: blueberry.number_of_factors()
            43
            sage: list(map(blueberry.number_of_factors, range(10)))
            [1, 6, 8, 7, 6, 5, 4, 3, 2, 1]
        """
        if n is None:
            length_word = self.word().length()
            num_factors = 1  # empty word
            for u, v, (i, j) in self.edge_iterator():
                if j is None:
                    num_factors += length_word - i
                else:
                    num_factors += j - i
        elif isinstance(n, (int, Integer)):
            num_factors = 0
            queue = [(0, 0)]
            while queue:
                v, l = queue.pop()
                if l == n:
                    num_factors += 1
                if l < n:
                    if self._transition_function[v]:
                        for (i, j), u in self._transition_function[v].items():
                            if j is None:
                                j = self.word().length()
                            if j - i >= n - l:
                                num_factors += 1
                            else:
                                queue.append((u, l + j - i + 1))
        else:
            raise TypeError("not an integer or None: %s" % n)
        return num_factors

    def factor_iterator(self, n=None):
        r"""
        Generate distinct factors of ``self``.

        INPUT:

        - ``n`` -- integer or ``None``

        OUTPUT:

        If ``n`` is an integer, returns an iterator over all distinct
        factors of length ``n``. If ``n`` is ``None``, returns an iterator
        generating all distinct factors.

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree
            sage: sorted( ImplicitSuffixTree(Word("cacao")).factor_iterator() )
            [word: , word: a, word: ac, word: aca, word: acao, word: ao, word: c, word: ca, word: cac, word: caca, word: cacao, word: cao, word: o]
            sage: sorted( ImplicitSuffixTree(Word("cacao")).factor_iterator(1) )
            [word: a, word: c, word: o]
            sage: sorted( ImplicitSuffixTree(Word("cacao")).factor_iterator(2) )
            [word: ac, word: ao, word: ca]
            sage: sorted( ImplicitSuffixTree(Word([0,0,0])).factor_iterator() )
            [word: , word: 0, word: 00, word: 000]
            sage: sorted( ImplicitSuffixTree(Word([0,0,0])).factor_iterator(2) )
            [word: 00]
            sage: sorted( ImplicitSuffixTree(Word([0,0,0])).factor_iterator(0) )
            [word: ]
            sage: sorted( ImplicitSuffixTree(Word()).factor_iterator() )
            [word: ]
            sage: sorted( ImplicitSuffixTree(Word()).factor_iterator(2) )
            []
        """
        # Every factor is a prefix of a suffix, so we do a depth
        # first search of the implicit suffix tree of the word.
        w = self.word()
        wlen = self.word().length()
        if n is None:
            queue = [(0, 0, -1, 0)]
            yield w[0:0]
            while queue:
                v, i, j, l = queue.pop()
                for k in range(i, j+1):
                    yield w[j-l:k]
                for (i, j), u in self._transition_function[v].items():
                    if j is None:
                        j = wlen
                    queue.append((u, i, j, l+j-i+1))
        elif isinstance(n, (int, Integer)):
            queue = [(0, 0, -1, 0)]
            while queue:
                v, i, j, l = queue.pop()
                if l == n:
                    yield w[j-l:j]
                if l < n:
                    for (i, j), u in self._transition_function[v].items():
                        if j is None:
                            j = wlen
                        if j - i >= n - l:
                            yield w[i-l-1:i-l+n-1]
                        else:
                            queue.append((u, i, j, l+j-i+1))
        else:
            raise TypeError("not an integer or None: %s" % n)

    def LZ_decomposition(self):
        r"""
        Return a list of index of the beginning of the block of the Lempel-Ziv
        decomposition of ``self.word``

        The *Lempel-Ziv decomposition* is the factorisation `u_1...u_k` of a
        word `w=x_1...x_n` such that `u_i` is the longest prefix of `u_i...u_k`
        that has an occurrence starting before `u_i` or a letter if this prefix
        is empty.

        OUTPUT:

        Return a list ``iB`` of index such that the blocks of the decomposition
        are ``self.word()[iB[k]:iB[k+1]]``

        EXAMPLES::

            sage: w = Word('abababb')
            sage: T = w.suffix_tree()
            sage: T.LZ_decomposition()
            [0, 1, 2, 6, 7]
            sage: w = Word('abaababacabba')
            sage: T = w.suffix_tree()
            sage: T.LZ_decomposition()
            [0, 1, 2, 3, 6, 8, 9, 11, 13]
            sage: w = Word([0, 0, 0, 1, 1, 0, 1])
            sage: T = w.suffix_tree()
            sage: T.LZ_decomposition()
            [0, 1, 3, 4, 5, 7]
            sage: w = Word('0000100101')
            sage: T = w.suffix_tree()
            sage: T.LZ_decomposition()
            [0, 1, 4, 5, 9, 10]
        """
        iB = [0]
        i = 0
        w = self.word()
        while i < len(w):
            l = 0
            (x, y), successor = self._find_transition(0, w[i])
            x = x-1
            while x < i+l:
                if y is None:
                    l = len(w)-i
                else:
                    l += y-x
                if i+l >= len(w):
                    l = len(w)-i
                    break
                (x, y), successor = self._find_transition(successor, w[i+l])
                x = x-1
            i += max(1, l)
            iB.append(i)
        return iB

    def _count_and_skip(self, node, i, j):
        r"""
        Use count and skip trick to follow the path starting at ``node`` and
        reading ``self.word()[i:j]``. We assume that reading
        ``self.word()[i:j]`` is possible from ``node``

        INPUT:

        - ``node`` -- explicit node of ``self``
        - ``i`` -- beginning of factor ``T.word()[i:j]``
        - ``j`` -- end of factor ``T.word()[i:j]``

        OUTPUT:

        The node obtained by starting at ``node`` and following the edges
        labeled by the letters of ``T.word()[i:j]``.
        Return ``("explicit", end_node)`` if w ends at the node ``end_node``,
        and ``("implicit", edge, d)`` if it ends after reading ``d`` letters along
        the edge ``edge``.

        EXAMPLES::

            sage: T = Word('00110111011').suffix_tree()
            sage: T._count_and_skip(5, 2, 5)
            ('implicit', (9, 10), 2)
            sage: T._count_and_skip(0, 1, 4)
            ('explicit', 7)
            sage: T._count_and_skip(0, 8, 10)
            ('implicit', (2, 7), 1)
            sage: T = Word('cacao').suffix_tree()
            sage: T._count_and_skip(3, 2, 5)
            ('explicit', 1)
        """
        trans = self._find_transition(node, self._letters[i])
        while trans[0][1] is not None and trans[0][1] - trans[0][0] + 1 <= j - i:
            node = trans[1]
            i += trans[0][1] - trans[0][0] + 1
            if i == j:
                return ('explicit', node)
            else:
                trans = self._find_transition(node, self._letters[i])
        if trans[0][1] is None and len(self.word()) - trans[0][0] + 1 <= j - i:
            return ('explicit', trans[1])
        else:
            return ('implicit', (node, trans[1]), j - i)

    def suffix_walk(self, edge, l):
        r"""
        Return the state of "w" if the input state is "aw".

        If the input state ``(edge, l)`` is path labeled "aw" with "a" a letter, the output is
        the state which is path labeled "w".

        INPUT:

        - ``edge`` -- the edge containing the state
        - ``l`` -- the string-depth of the state on edge (``l``>0)

        OUTPUT:

        Return ``("explicit", end_node)`` if the state of w is an explicit
        state and ``("implicit", edge, d)`` is obtained by reading ``d``
        letters on ``edge``.

        EXAMPLES::

            sage: T = Word('00110111011').suffix_tree()
            sage: T.suffix_walk((0, 5), 1)
            ('explicit', 0)
            sage: T.suffix_walk((7, 3), 1)
            ('implicit', (9, 4), 1)
        """
        start, end = edge
        # Select the transition that corresponds to edge
        ij = next(ij for ij, target in self._transition_function[start].items()
                  if target == end)

        # self.word()[i-1:j] is the word on the edges
        i = ij[0] - 1
        parent = self.suffix_link(start)
        return self._count_and_skip(parent, i, i + l)

    def leftmost_covering_set(self):
        r"""
        Compute the leftmost covering set of square pairs in ``self.word()``.
        Return a square as a pair ``(i,l)`` designating factor
        ``self.word()[i:i+l]``.

        A  leftmost covering set is a set such that the leftmost occurrence
        `(j,l)` of a square in ``self.word()`` is covered by a pair
        `(i,l)` in the set for all types of squares. We say that `(j,l)` is
        covered by `(i,l)` if `(i,l)` (i+1,l), \ldots, (j,l)` are all
        squares.

        The set is returned in the form of a list ``P`` such that ``P[i]``
        contains all the lengths of squares starting at ``i`` in the set.
        The lists ``P[i]`` are sorted in decreasing order.

        The algorithm used is described in [DS2004]_.

        EXAMPLES::

            sage: w = Word('abaabaabbaaabaaba')
            sage: T = w.suffix_tree()
            sage: T.leftmost_covering_set()
            [[6], [6], [2], [], [], [], [], [2], [], [], [6, 2], [], [], [], [], [], []]
            sage: w = Word('abaca')
            sage: T = w.suffix_tree()
            sage: T.leftmost_covering_set()
            [[], [], [], [], []]
            sage: T = Word('aaaaa').suffix_tree()
            sage: T.leftmost_covering_set()
            [[4, 2], [], [], [], []]
        """

        def condition1_square_pairs(i):
            r"""
            Compute the squares that have their center (the last letter of the
            first  occurrence of ``w`` in ``ww``) in the `i`-th block of the
            LZ-decomposition and that start in the `i`-th block and end in the
            `(i+1)`-th.
            """
            for k in range(1, B[i+1]-B[i]+1):
                q = B[i+1]-k
                k1 = w.longest_forward_extension(B[i+1], q) if B[i+1] < len(w) else 0
                k2 = w.longest_backward_extension(B[i+1]-1, q-1) if q > 0 else 0
                start = max(q-k2, q-k+1)
                if k1+k2 >= k and k1 > 0 and start >= B[i]:
                    yield (start, 2*k)

        def condition2_square_pairs(i):
            r"""
            Compute the squares that have their center (the last letter of the
            first  occurrence of ``w``  in ``ww``)  in the `i`-th  block of the
            LZ-decomposition and that starts in the `(i-1)`-th block or before.
            Their end is either in the `i`-th or the `(i+1)`-th block.
            """
            if i+2 < len(B):
                end = B[i+2] - B[i] + 1
            else:
                end = B[i+1] - B[i] + 1
            for k in range(2, end):
                q = B[i] + k
                k1 = w.longest_forward_extension(B[i], q) if q < len(w) else 0
                k2 = w.longest_backward_extension(B[i]-1, q-1) if B[i] > 0 else 0
                start = max(B[i]-k2, B[i]-k+1)
                if k1+k2 >= k and k1 > 0 and start+k <= B[i+1] and k2 > 0:
                    yield (start, 2*k)

        w = self.word()
        B = self.LZ_decomposition()
        P = [[] for _ in w]
        for i in range(len(B) - 1):
            for i, l in chain(condition2_square_pairs(i), condition1_square_pairs(i)):
                P[i].append(l)
        for l in P:
            l.reverse()
        return P

    #####
    # Miscellaneous methods
    #####

    def uncompactify(self):
        r"""
        Return the tree obtained from ``self`` by splitting edges so that they
        are labelled by exactly one letter.

        The resulting tree is isomorphic to the suffix trie.

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree, SuffixTrie
            sage: abbab = Words("ab")("abbab")
            sage: s = SuffixTrie(abbab)
            sage: t = ImplicitSuffixTree(abbab)
            sage: t.uncompactify().is_isomorphic(s.to_digraph())                        # needs sage.graphs
            True
        """
        tree = self.to_digraph(word_labels=True)
        newtree = DiGraph()
        newtree.add_vertices(range(tree.order()))
        new_node = tree.order() + 1
        for u, v, label in tree.edge_iterator():
            if len(label) == 1:
                newtree.add_edge(u, v)
            else:
                newtree.add_edge(u, new_node, label[0])
                for w in label[1:-1]:
                    newtree.add_edge(new_node, new_node+1, w)
                    new_node += 1
                newtree.add_edge(new_node, v, label[-1])
                new_node += 1
        return newtree

    def trie_type_dict(self):
        r"""
        Return a dictionary in a format compatible with that of the suffix
        trie transition function.

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import ImplicitSuffixTree, SuffixTrie
            sage: W = Words("ab")
            sage: t = ImplicitSuffixTree(W("aba"))
            sage: d = t.trie_type_dict()
            sage: len(d)
            5
            sage: d                     # random
            {(4, word: b): 5, (0, word: a): 4, (0, word: b): 3, (5, word: a): 1, (3, word: a): 2}
        """
        d = {}
        new_node = len(self._transition_function)
        for u, dd in self._transition_function.items():
            for sl, v in dd.items():
                w = self._word[sl[0]-1:sl[1]]
                if w.length() == 1:
                    d[u, w] = v
                else:
                    d[u, w[0:1]] = new_node
                    for i in range(1, w.length()-1):
                        d[new_node, w[i:i+1]] = new_node + 1
                        new_node += 1
                    d[new_node, w[-1:]] = v
                    new_node += 1
        return d

################################################################################
# Decorated Suffix Tree
################################################################################


class DecoratedSuffixTree(ImplicitSuffixTree):
    r"""
    The decorated suffix tree of a word.

    A *decorated suffix tree* of a word `w` is the suffix tree of `w`
    marked with the end point of all squares in the `w`.

    The symbol ``$`` is appended to ``w`` to ensure that each final
    state is a leaf of the suffix tree.

    INPUT:

    - ``w`` -- a finite word

    EXAMPLES::

        sage: from sage.combinat.words.suffix_trees import DecoratedSuffixTree
        sage: w = Word('0011001')
        sage: DecoratedSuffixTree(w)
        Decorated suffix tree of : 0011001$
        sage: w = Word('0011001', '01')
        sage: DecoratedSuffixTree(w)
        Decorated suffix tree of : 0011001$

    ALGORITHM:

    When using ``'pair'`` as output, the squares are retrieved in linear
    time. The algorithm is an implementation of the one proposed in
    [DS2004]_.
    """
    def __init__(self, w):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import DecoratedSuffixTree
            sage: w = Word('0011001')
            sage: DST = DecoratedSuffixTree(w)

        We skip the ``_test_and_split`` test because it is not a test meant
        for the ``TestSuite``::

            sage: TestSuite(DST).run(skip='_test_and_split')

        Test that we do not allow ``'$'`` to appear in the word::

            sage: w = Word('0011001$')
            sage: DecoratedSuffixTree(w)
            Traceback (most recent call last):
            ...
            ValueError: the symbol '$' is reserved for this class
        """
        if "$" in w:
            raise ValueError("the symbol '$' is reserved for this class")
        end_symbol = '$'
        w = Word(str(w) + end_symbol)
        ImplicitSuffixTree.__init__(self, w)
        self.labeling = self._complete_labeling()

    def __repr__(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import DecoratedSuffixTree
            sage: w = Word('0011001')
            sage: t = DecoratedSuffixTree(w)
            sage: t.__repr__()
            'Decorated suffix tree of : 0011001$'
        """
        w = self.word()
        if len(w) > 40:
            w = str(w[:40]) + '...'
        return f"Decorated suffix tree of : {w}"

    def _partial_labeling(self):
        r"""
        Make a depth-first search in the suffix tree and mark some squares of a
        leftmost covering set of the tree.

        This is used by :meth:`_complete_labeling`.

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import DecoratedSuffixTree
            sage: w = Word('abaababbabba')
            sage: T = DecoratedSuffixTree(w)
            sage: T._partial_labeling()
            {(3, 4): [1], (5, 1): [3], (5, 6): [1], (11, 17): [1], (13, 8): [1], (15, 10): [2]}
        """
        def node_processing(node, parent, head):
            r"""
            Marks points along the edge ``(parent, node)`` if the string depth
            of parent is smaller than the length of the square at the head of
            ``P(node)``.
            Make it for all such square pairs and remove them from ``P(node)``.

            INPUT:

            - ``node`` -- a node of ``self``
            - ``parent`` -- the parent of a node in ``self``
            - ``head`` -- tuple indicating the head of the list ``P(node)``

            OUTPUT: ``(i, pos)``, the new head of ``P(node)``
            """
            i, pos = head
            pano = (parent, node)
            while pos < len(P[i]) and P[i][pos] > string_depth[parent]:
                label = P[i][pos] - string_depth[parent]
                if pano in labeling:
                    labeling[pano].append(label)
                else:
                    labeling[pano] = [label]
                pos += 1
            return (i, pos)

        def treat_node(current_node, parent):
            r"""
            Proceed to a depth-first search in ``self``, counting the
            string_depth of each node and processing each node for marking.

            To initiate the depth first search call ``self.treat_node(0,None)``

            INPUT:

            - ``current_node`` -- a node
            - ``parent`` -- parent of ``current_node`` in ``self``

            OUTPUT:

            The resulting list P(current_node) with current_node have been
            processed by ``node_processing``. The output is a pair ``(i,
            pos)`` such that ``P[i][pos:]`` is the list of current_node.
            """

            # Call recursively on children of current_node
            if current_node in D:
                node_list = (n, 0)
                for child in D[current_node]:
                    i, j = D[current_node][child]
                    if j is None:
                        j = n
                    string_depth[child] = string_depth[current_node] + j - i
                    child_list = treat_node(child, current_node)
                    if child_list[0] < node_list[0]:
                        node_list = child_list
            else:  # The node is a child
                node_list = (n - string_depth[current_node], 0)
            # Make treatment on current node head
            return node_processing(current_node, parent, node_list)

        P = self.leftmost_covering_set()
        D = self.transition_function_dictionary()
        string_depth = {0: 0}
        n = len(self.word())
        labeling = {}
        treat_node(0, None)
        return labeling

    def _complete_labeling(self):
        r"""
        Return a dictionary of edges of ``self``, with marked points for the end
        of each distinct squares that can be found in ``self.word()``.

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import DecoratedSuffixTree
            sage: w = Word('aabbaaba')
            sage: DecoratedSuffixTree(w)._complete_labeling()
            {(2, 7): [1], (5, 4): [1]}
        """

        def walk_chain(u, v, l, start):
            r"""
            Execute a chain of suffix walk until a walk is unsuccessful or it
            got to a point already registered in ``QP``. Registers all visited
            point in ``Q``.

            INPUT:

            - ``(u, v)`` -- edge on which the point is registered
            - ``l`` -- depth of the registered point on (u,v)
            - ``start`` -- beginning of the squares registered by the label
              ``(u, v), l``
            """
            # Mark the point in labeling
            if (u, v) in labeling:
                labeling[(u, v)].append(l)
            else:
                labeling[(u, v)] = [l]
            # Make the walk
            final_state = self.suffix_walk((u, v), l)
            successful = False
            if final_state[0] == 'explicit':
                parent = final_state[1]
                transition = self._find_transition(parent, self._letters[start])
                if transition is not None:
                    child = transition[1]
                    successful = True
                    depth = 1
            else:
                parent = final_state[1][0]
                child = final_state[1][1]
                depth = final_state[2]
                next_letter = self._letters[D[parent][child][0]+depth]
                if next_letter == self._letters[start]:
                    successful = True
                    depth += 1
            # If needed start a new walk
            if successful:
                if (parent, child) in prelabeling:
                    if depth not in prelabeling[(parent, child)]:
                        walk_chain(parent, child, depth, start+1)
                else:
                    walk_chain(parent, child, depth, start+1)

        def treat_node(current_node, i, j):
            r"""
            Execute a depth-first search on ``self`` and start a suffix walk
            for labeled points on each edges of T.

            The function is recursive, call
            ``treat_node(0,0,0)`` to initiate the search.

            INPUT:

            - ``current_node`` -- the node to treat
            - ``(i, j)`` -- pair of index such that the path from 0 to
              ``current_node`` reads ``self.word()[i:j]``
            """

            if current_node in D:
                for child in D[current_node]:
                    edge = (current_node, child)
                    edge_label = D[edge[0]][edge[1]]
                    treat_node(child, edge_label[0]-(j-i), edge_label[1])
                    if (current_node, child) in prelabeling:
                        for l in prelabeling[edge]:
                            square_start = edge_label[0] - (j - i)
                            walk_chain(current_node, child, l, square_start)

        prelabeling = self._partial_labeling()
        labeling = {}
        D = self.transition_function_dictionary()
        treat_node(0, 0, 0)
        return labeling

    def square_vocabulary(self, output='pair'):
        r"""
        Return the list of distinct squares of ``self.word``.

        Two types of outputs are available `pair` and `word`. The algorithm
        is only truly linear if `output` is set to `pair`. A pair is a tuple
        `(i, l)` that indicates the factor ``self.word()[i:i+l]``.
        The option ``'word'`` return word objects.

        INPUT:

        - ``output`` -- (default: ``'pair'``) either ``'pair'`` or ``'word'``

        EXAMPLES::

            sage: from sage.combinat.words.suffix_trees import DecoratedSuffixTree
            sage: w = Word('aabb')
            sage: sorted(DecoratedSuffixTree(w).square_vocabulary())
            [(0, 0), (0, 2), (2, 2)]
            sage: w = Word('00110011010')
            sage: sorted(DecoratedSuffixTree(w).square_vocabulary(output='word'))
            [word: , word: 00, word: 00110011, word: 01100110, word: 1010, word: 11]
        """
        def treat_node(current_node, i, j):
            if current_node in D:
                for child in D[current_node]:
                    edge = (current_node, child)
                    edge_label = (D[edge[0]][edge[1]])
                    treat_node(child, edge_label[0]-(j-i), edge_label[1])
                    if (current_node, child) in Q:
                        for l in Q[(current_node, child)]:
                            square_start = edge_label[0]-(j-i)
                            pair = (square_start, edge_label[0]+l-square_start)
                            squares.append(pair)

        if output not in ["pair", "word"]:
            msg = f"output should be 'pair' or 'word'; got {output}"
            raise ValueError(msg)
        D = self.transition_function_dictionary()
        Q = self.labeling
        squares = [(0, 0)]
        treat_node(0, 0, 0)
        if output == "pair":
            return squares
        return [self.word()[i:i + l] for i, l in squares]
