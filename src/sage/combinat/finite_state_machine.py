# -*- coding: utf-8 -*-
"""
Finite State Machines, Automata, Transducers

This module adds support for finite state machines, automata and
transducers. See class :class:`FiniteStateMachine` and the examples
below for details creating one.

Examples
========


A simple finite state machine
-----------------------------

We can easily create a finite state machine by

::

    sage: fsm = FiniteStateMachine()
    sage: fsm
    Finite state machine with 0 states

By default this is the empty finite state machine, so not very
interesting. Let's create some states and transitions::

    sage: from sage.combinat.finite_state_machine import FSMState, FSMTransition
    sage: day = FSMState('day')
    sage: night = FSMState('night')
    sage: sunrise = FSMTransition(night, day)
    sage: sunset = FSMTransition(day, night)

And now let's add those states and transitions to our finite state machine::

    sage: fsm.add_transition(sunrise)
    Transition from 'night' to 'day': -|-
    sage: fsm.add_transition(sunset)
    Transition from 'day' to 'night': -|-

Note that the states are added automatically, since they are present
in the transitions. We could add the states manually by

::

    sage: fsm.add_state(day)
    'day'
    sage: fsm.add_state(night)
    'night'

Anyhow, we got the following finite state machine::

    sage: fsm
    Finite state machine with 2 states

We can also visualize it as a graph by

::

    sage: fsm.graph()
    Digraph on 2 vertices

Alternatively, we could have created the finite state machine above
simply by

::

    sage: FiniteStateMachine([('night', 'day'), ('day', 'night')])
    Finite state machine with 2 states

or by

::

    sage: fsm = FiniteStateMachine()
    sage: day = fsm.add_state('day')
    sage: night = fsm.add_state('night')
    sage: sunrise = fsm.add_transition(night, day)
    sage: sunset = fsm.add_transition(day, night)
    sage: fsm
    Finite state machine with 2 states

A simple Automaton (recognizing NAFs)
---------------------------------------

We want to build an automaton which recognizes non-adjacent forms
(NAFs), i.e., sequences which have no adjacent non-zeros.
We use `0`, `1`, and `-1` as digits::

    sage: NAF = Automaton(
    ....:     {'A': [('A', 0), ('B', 1), ('B', -1)], 'B': [('A', 0)]})
    sage: NAF.state('A').is_initial = True
    sage: NAF.state('A').is_final = True
    sage: NAF.state('B').is_final = True
    sage: NAF
    Automaton with 2 states

Of course, we could have specified the initial and final states
directly in the definition of ``NAF`` by ``initial_states=['A']`` and
``final_states=['A', 'B']``.

So let's test the automaton with some input::

    sage: NAF([0])[0]
    True
    sage: NAF([0, 1])[0]
    True
    sage: NAF([1, -1])[0]
    False
    sage: NAF([0, -1, 0, 1])[0]
    True
    sage: NAF([0, -1, -1, -1, 0])[0]
    False
    sage: NAF([-1, 0, 0, 1, 1])[0]
    False

Alternatively, we could call that by

::

    sage: NAF.process([-1, 0, 0, 1, 1])[0]
    False

A simple transducer (binary inverter)
-------------------------------------

Let's build a simple transducer, which rewrites a binary word by
iverting each bit::

    sage: inverter = Transducer({'A': [('A', 0, 1), ('A', 1, 0)]},
    ....:     initial_states=['A'], final_states=['A'])

We can look at the states and transitions::

    sage: inverter.states()
    ['A']
    sage: for t in inverter.transitions():
    ....:     print t
    Transition from 'A' to 'A': 0|1
    Transition from 'A' to 'A': 1|0

Now we apply a word to it and see what the transducer does::

    sage: inverter([0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1])
    (True, 'A', [1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0])

``True`` means, that we landed in a final state, that state is labeled
``'A'``, and we also got an output.


A transducer which performs division by `3` in binary
-----------------------------------------------------

Now we build a transducer, which divides a binary number by 3.
The labels of the states are the remainder of the division.
The transition function is

::

    sage: def f(state_from, read):
    ....:     if state_from + read <= 1:
    ....:         state_to = 2*state_from + read
    ....:         write = 0
    ....:     else:
    ....:         state_to = 2*state_from + read - 3
    ....:         write = 1
    ....:     return (state_to, write)

We get the transducer with

::

    sage: D = Transducer(f, initial_states=[0], final_states=[0],
    ....:                input_alphabet=[0, 1])

Now we want to divide 13 by 3::

    sage: D([1, 1, 0, 1])
    (False, 1, [0, 1, 0, 0])

So we have 13 : 3 = 4 and the reminder is 1. ``False`` means 13 is not
divisible by 3.


Using the hook-functions
------------------------

Let's use the previous example "divison by `3`" to demonstrate the
optional state and transition parameters ``hook``.

First, we define, what those functions should do. In our case, this is
just saying in which state we are and which transition we take

::

    sage: def state_hook(state, process):
    ....:     print "We are now in State %s." % (state.label(),)
    sage: from sage.combinat.finite_state_machine import FSMWordSymbol
    sage: def transition_hook(transition, process):
    ....:     print ("Currently we go from %s to %s, "
    ....:            "reading %s and writing %s." % (
    ....:                transition.from_state, transition.to_state,
    ....:                FSMWordSymbol(transition.word_in),
    ....:                FSMWordSymbol(transition.word_out)))

Now, let's add these hook-functions to the existing transducer::

    sage: for s in D.iter_states():
    ....:     s.hook = state_hook
    sage: for t in D.iter_transitions():
    ....:     t.hook = transition_hook

Rerunning the process again now gives the following output::

    sage: D.process([1, 1, 0, 1])
    We are now in State 0.
    Currently we go from 0 to 1, reading 1 and writing 0.
    We are now in State 1.
    Currently we go from 1 to 0, reading 1 and writing 1.
    We are now in State 0.
    Currently we go from 0 to 0, reading 0 and writing 0.
    We are now in State 0.
    Currently we go from 0 to 1, reading 1 and writing 0.
    We are now in State 1.
    (False, 1, [0, 1, 0, 0])

The example above just explains the basic idea of using
hook-functions. In the following, we will use those hooks more seriously.


Detecting sequences with same number of `0` and `1`
---------------------------------------------------

Suppose we have a binary input and want to accept all sequences with
the same number of `0` and `1`. This cannot be done with a finite
automaton. Anyhow, we can make usage of the hook functions to extend
our finite automaton by a counter::

    sage: from sage.combinat.finite_state_machine import FSMState, FSMTransition
    sage: C = Automaton()
    sage: def update_counter(state, process):
    ....:     l = process.read_letter()
    ....:     process.fsm.counter += 1 if l == 1 else -1
    ....:     if process.fsm.counter > 0:
    ....:         next_state = 'positive'
    ....:     elif process.fsm.counter < 0:
    ....:         next_state = 'negative'
    ....:     else:
    ....:         next_state = 'zero'
    ....:     return FSMTransition(state, process.fsm.state(next_state),
    ....:                          l, process.fsm.counter)
    sage: C.add_state(FSMState('zero', hook=update_counter,
    ....:             is_initial=True, is_final=True))
    'zero'
    sage: C.add_state(FSMState('positive', hook=update_counter))
    'positive'
    sage: C.add_state(FSMState('negative', hook=update_counter))
    'negative'

Now, let's input some sequence::

    sage: C.counter = 0; C([1, 1, 1, 1, 0, 0])
    (False, 'positive', [1, 2, 3, 4, 3, 2])

The result is False, since there are four `1` but only two `0`. We
land in the state ``positive`` and we can also see the values of the
counter in each step.

Let's try some other examples::

    sage: C.counter = 0; C([1, 1, 0, 0])
    (True, 'zero', [1, 2, 1, 0])
    sage: C.counter = 0; C([0, 1, 0, 0])
    (False, 'negative', [-1, 0, -1, -2])


AUTHORS:

- Daniel Krenn (2012-03-27): initial version
- Clemens Heuberger (2012-04-05): initial version
- Sara Kropf (2012-04-17): initial version
- Clemens Heuberger (2013-08-21): release candidate for Sage patch
- Daniel Krenn (2013-08-21): release candidate for Sage patch
- Sara Kropf (2013-08-21): release candidate for Sage patch
- Clemens Heuberger (2013-09-02): documentation improved
- Daniel Krenn (2013-09-13): comments from trac worked in
- Clemens Heuberger (2013-11-03): output (labels) of determinisation,
    product, composition, etc. changed (for consistency),
    representation of state changed, documentation improved
- Daniel Krenn (2013-11-04): whitespaces in documentation corrected
- Clemens Heuberger (2013-11-04): full_group_by added
- Daniel Krenn (2013-11-04): next release candidate for Sage patch
- Sara Kropf (2013-11-08): fix for adjacency matrix
- Clemens Heuberger (2013-11-11): fix for prepone_output
- Daniel Krenn (2013-11-11): comments from trac 15078 included:
    docstring of FiniteStateMachine rewritten, Automaton and Transducer
    inherited from FiniteStateMachine
- Daniel Krenn (2013-11-25): documentation improved according to
    comments from trac 15078

ACKNOWLEDGEMENT:

- Daniel Krenn, Clemens Heuberger and Sara Kropf are supported by the
  Austrian Science Fund (FWF): P 24644-N26.

"""

#*****************************************************************************
#  Copyright (C) 2012, 2013 Daniel Krenn <math+sage@danielkrenn.at>
#                2012, 2013 Clemens Heuberger <clemens.heuberger@aau.at>
#                2012, 2013 Sara Kropf <sara.kropf@aau.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject
from sage.graphs.digraph import DiGraph
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.calculus.var import var
from sage.misc.latex import latex
from sage.functions.trig import cos, sin, atan2
from sage.symbolic.constants import pi

from copy import copy
from copy import deepcopy

import itertools
from collections import defaultdict


def full_group_by(l, key=lambda x: x):
    """
    Group iterable ``l`` by values of ``key``.

    INPUT:

    - iterable ``l``
    - key function ``key``

    OUTPUT:

    A list of pairs ``(k, elements)`` such that ``key(e)=k`` for all
    ``e`` in ``elements``.

    This is similar to ``itertools.groupby`` except that lists are
    returned instead of iterables and no prior sorting is required.

    We do not require

    - that the keys are sortable (in contrast to the
      approach via ``sorted`` and ``itertools.groupby``) and
    - that the keys are hashable (in contrast to the
      implementation proposed in `<http://stackoverflow.com/a/15250161>`_).

    However, it is required

    - that distinct keys have distinct ``str``-representations.

    The implementation is inspired by
    `<http://stackoverflow.com/a/15250161>`_, but non-hashable keys are
    allowed.

    EXAMPLES::

        sage: from sage.combinat.finite_state_machine import full_group_by
        sage: t = [2/x, 1/x, 2/x]
        sage: r = full_group_by([0,1,2], key=lambda i:t[i])
        sage: sorted(r, key=lambda p:p[1])
        [(2/x, [0, 2]), (1/x, [1])]
        sage: from itertools import groupby
        sage: for k, elements in groupby(sorted([0,1,2],
        ....:                            key=lambda i:t[i]),
        ....:                            key=lambda i:t[i]):
        ....:     print k, list(elements)
        2/x [0]
        1/x [1]
        2/x [2]

    Note that the behavior is different from ``itertools.groupby``
    because neither `1/x<2/x` nor `2/x<1/x` does hold.

    Here, the result ``r`` has been sorted in order to guarantee a
    consistent order for the doctest suite.
    """
    elements = defaultdict(list)
    original_keys = {}
    for item in l:
        k = key(item)
        s = str(k)
        if s in original_keys:
            if original_keys[s]!=k:
                raise ValueError("Two distinct elements with representation "
                                 "%s " % s)
        else:
            original_keys[s]=k
        elements[s].append(item)
    return [(original_keys[s], values ) for (s, values) in elements.items()]

#*****************************************************************************

FSMEmptyWordSymbol = '-'

def FSMLetterSymbol(letter):
    """
    Returns a string associated to the input letter.

    INPUT:

    - ``letter`` -- the input letter or ``None`` (representing the
      empty word).

    OUTPUT:

    If ``letter`` is ``None`` the symbol for the empty word
    ``FSMEmptyWordSymbol`` is returned, otherwise the string
    associated to the letter.

    EXAMPLES::

        sage: from sage.combinat.finite_state_machine import FSMLetterSymbol
        sage: FSMLetterSymbol(0)
        '0'
        sage: FSMLetterSymbol(None)
        '-'
    """
    return FSMEmptyWordSymbol if letter is None else repr(letter)


def FSMWordSymbol(word):
    """
    Returns a string of ``word``. It may returns the symbol of the
    empty word ``FSMEmptyWordSymbol``.

    INPUT:

    - ``word`` -- the input word.

    OUTPUT:

    A string of ``word``.

    EXAMPLES::

        sage: from sage.combinat.finite_state_machine import FSMWordSymbol
        sage: FSMWordSymbol([0, 1, 1])
        '0,1,1'
    """
    if not isinstance(word, list):
        return FSMLetterSymbol(word)
    if len(word) == 0:
        return FSMEmptyWordSymbol
    s = ''
    for letter in word:
        s += (',' if len(s) > 0 else '') + FSMLetterSymbol(letter)
    return s


#*****************************************************************************


def is_FSMState(S):
    """
    Tests whether or not ``S`` inherits from :class:`FSMState`.

    TESTS::

        sage: from sage.combinat.finite_state_machine import is_FSMState, FSMState
        sage: is_FSMState(FSMState('A'))
        True
    """
    return isinstance(S, FSMState)


class FSMState(SageObject):
    """
    Class for a state of a finite state machine.

    INPUT:

    - ``label`` -- the label of the state.

    - ``word_out`` -- (default: ``None``) a word that is written when
      the state is reached.

    - ``is_initial`` -- (default: ``False``)

    - ``is_final`` -- (default: ``False``)

    - ``hook`` -- (default: ``None``) A function which is called when
      the state is reached during processing input.

    OUTPUT:

    Returns a state of a finite state machine.

    EXAMPLES::

        sage: from sage.combinat.finite_state_machine import FSMState
        sage: A = FSMState('state 1', word_out=0, is_initial=True)
        sage: A
        'state 1'
        sage: A.label()
        'state 1'
        sage: B = FSMState('state 2')
        sage: A == B
        False

    """
    def __init__(self, label, word_out=None,
                 is_initial=False, is_final=False,
                 hook=None):
        """
        See :class:`FSMState` for more information.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: FSMState('final', is_final=True)
            'final'
        """
        if label is None or label == "":
            raise ValueError, "You have to specify a label for the state."
        self._label_ = label

        if isinstance(word_out, list):
            self.word_out = word_out
        elif word_out is not None:
            self.word_out = [word_out]
        else:
            self.word_out = []

        self.is_initial = is_initial
        self.is_final = is_final
        if hook is not None:
            if hasattr(hook, '__call__'):
                self.hook = hook
            else:
                raise TypeError, 'Wrong argument for hook.'


    def label(self):
        """
        Returns the label of the state.

        INPUT:

        Nothing.

        OUTPUT:

        The label of the state.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('state')
            sage: A.label()
            'state'
        """
        return self._label_


    def __copy__(self):
        """
        Returns a (shallow) copy of the state.

        INPUT:

        Nothing.

        OUTPUT:

        A new state.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A')
            sage: copy(A)
            'A'
        """
        new = FSMState(self.label(), self.word_out,
                       self.is_initial, self.is_final)
        if hasattr(self, 'hook'):
            new.hook = self.hook
        return new


    copy = __copy__


    def __deepcopy__(self, memo):
        """
        Returns a deep copy of the state.

        INPUT:

        - ``memo`` -- a dictionary storing already processed elements.

        OUTPUT:

        A new state.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A')
            sage: deepcopy(A)
            'A'
        """
        try:
            label = self._deepcopy_relabel_
        except AttributeError:
            label = deepcopy(self.label(), memo)
        new = FSMState(label, deepcopy(self.word_out, memo),
                       self.is_initial, self.is_final)
        if hasattr(self, 'hook'):
            new.hook = deepcopy(self.hook, memo)
        return new


    def deepcopy(self, memo=None):
        """
        Returns a deep copy of the state.

        INPUT:

        - ``memo`` -- (default: ``None``) a dictionary storing already
          processed elements.

        OUTPUT:

        A new state.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A')
            sage: deepcopy(A)
            'A'
        """
        return deepcopy(self, memo)


    def relabeled(self, label, memo=None):
        """
        Returns a deep copy of the state with a new label.

        INPUT:

        - ``label`` -- the label of new state.

        - ``memo`` -- (default: ``None``) a dictionary storing already
          processed elements.

        OUTPUT:

        A new state.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A')
            sage: A.relabeled('B')
            'B'

        """
        self._deepcopy_relabel_ = label
        new = deepcopy(self, memo)
        del self._deepcopy_relabel_
        return new


    def __hash__(self):
        """
        Returns a hash value for the object.

        INPUT:

        Nothing.

        OUTPUT:

        The hash of this state.

        TESTS::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A')
            sage: hash(A) #random
            -269909568
        """
        return hash(self.label())


    def _repr_(self):
        """
        Returns the string "label".

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        TESTS:

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: FSMState('A')._repr_()
            "'A'"
        """
        return repr(self.label())


    def __eq__(left, right):
        """
        Returns True if two states are the same, i.e., if they have
        the same labels.

        Note that the hooks and whether the states are initial or
        final are not checked.

        INPUT:

        - ``left`` -- a state.

        - ``right`` -- a state.

        OUTPUT:

        True or False.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A')
            sage: B = FSMState('A', is_initial=True)
            sage: A == B
            True
        """
        if not is_FSMState(right):
            return False
        return left.label() == right.label()


    def __ne__(left, right):
        """
        Tests for inequality, complement of __eq__.

        INPUT:

        - ``left`` -- a state.

        - ``right`` -- a state.

        OUTPUT:

        True or False.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A', is_initial=True)
            sage: B = FSMState('A', is_final=True)
            sage: A != B
            False
        """
        return (not (left == right))


    def __nonzero__(self):
        """
        Returns True.

        INPUT:

        Nothing.

        OUTPUT:

        True or False.

        TESTS::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: FSMState('A').__nonzero__()
            True
        """
        return True  # A state cannot be zero (see __init__)


#*****************************************************************************


def is_FSMTransition(T):
    """
    Tests whether or not ``T`` inherits from :class:`FSMTransition`.

    TESTS::

        sage: from sage.combinat.finite_state_machine import is_FSMTransition, FSMTransition
        sage: is_FSMTransition(FSMTransition('A', 'B'))
        True
    """
    return isinstance(T, FSMTransition)


class FSMTransition(SageObject):
    """
    Class for a transition of a finite state machine.

    INPUT:

    - ``from_state`` -- state from which transition starts.

    - ``to_state`` -- state in which transition ends.

    - ``word_in`` -- the input word of the transitions (when the
      finite state machine is used as automaton)

    - ``word_out`` -- the output word of the transitions (when the
      finite state machine is used as transducer)

    OUTPUT:

    A transition of a finite state machine.

    EXAMPLES::

        sage: from sage.combinat.finite_state_machine import FSMState, FSMTransition
        sage: A = FSMState('A')
        sage: B = FSMState('B')
        sage: S = FSMTransition(A, B, 0, 1)
        sage: T = FSMTransition('A', 'B', 0, 1)
        sage: T == S
        True
        sage: U = FSMTransition('A', 'B', 0)
        sage: U == T
        False

    """
    def __init__(self, from_state, to_state,
                 word_in=None, word_out=None,
                 hook=None):
        """
        See :class:`FSMTransition` for more information.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: FSMTransition('A', 'B', 0, 1)
            Transition from 'A' to 'B': 0|1
        """
        if is_FSMState(from_state):
            self.from_state = from_state
        else:
            self.from_state = FSMState(from_state)
        if is_FSMState(to_state):
            self.to_state = to_state
        else:
            self.to_state = FSMState(to_state)

        if isinstance(word_in, list):
            self.word_in = word_in
        elif word_in is not None:
            self.word_in = [word_in]
        else:
            self.word_in = []

        if isinstance(word_out, list):
            self.word_out = word_out
        elif word_out is not None:
            self.word_out = [word_out]
        else:
            self.word_out = []

        if hook is not None:
            if hasattr(hook, '__call__'):
                self.hook = hook
            else:
                raise TypeError, 'Wrong argument for hook.'


    def __copy__(self):
        """
        Returns a (shallow) copy of the transition.

        INPUT:

        Nothing.

        OUTPUT:

        A new transition.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: t = FSMTransition('A', 'B', 0)
            sage: copy(t)
            Transition from 'A' to 'B': 0|-
        """
        new = FSMTransition(self.from_state, self.to_state,
                            self.word_in, self.word_out)
        if hasattr(self, 'hook'):
            new.hook = self.hook
        return new


    copy = __copy__

    def __deepcopy__(self, memo):
        """
        Returns a deep copy of the transition.

        INPUT:

        - ``memo`` -- a dictionary storing already processed elements.

        OUTPUT:

        A new transition.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: t = FSMTransition('A', 'B', 0)
            sage: deepcopy(t)
            Transition from 'A' to 'B': 0|-
        """
        new = FSMTransition(deepcopy(self.from_state, memo),
                            deepcopy(self.to_state, memo),
                            deepcopy(self.word_in, memo),
                            deepcopy(self.word_out, memo))
        if hasattr(self, 'hook'):
            new.hook = deepcopy(self.hook, memo)
        return new


    def deepcopy(self, memo=None):
        """
        Returns a deep copy of the transition.

        INPUT:

        - ``memo`` -- (default: ``None``) a dictionary storing already
          processed elements.

        OUTPUT:

        A new transition.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: t = FSMTransition('A', 'B', 0)
            sage: deepcopy(t)
            Transition from 'A' to 'B': 0|-
        """
        return deepcopy(self, memo)


    def __hash__(self):
        """
        Since transitions are mutable, they should not be hashable, so
        we return a type error.

        INPUT:

        Nothing.

        OUTPUT:

        The hash of this transition.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: hash(FSMTransition('A', 'B'))
            Traceback (most recent call last):
            ...
            TypeError: Transitions are mutable, and thus not hashable.

        """
        raise TypeError, "Transitions are mutable, and thus not hashable."


    def _repr_(self):
        """
        Represents a transitions as from state to state and input, output.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: FSMTransition('A', 'B', 0, 0)._repr_()
            "Transition from 'A' to 'B': 0|0"

        """
        return "Transition from %s to %s: %s" % (repr(self.from_state),
                                                 repr(self.to_state),
                                                 self._in_out_label_())


    def _in_out_label_(self):
        """
        Returns the input and output of a transition as
        "word_in|word_out".

        INPUT:

        Nothing.

        OUTPUT:

        A string of the input and output labels.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: FSMTransition('A', 'B', 0, 1)._in_out_label_()
            '0|1'
        """
        return "%s|%s" % (FSMWordSymbol(self.word_in),
                          FSMWordSymbol(self.word_out))


    def __eq__(left, right):
        """
        Returns True if the two transitions are the same, i.e., if the
        both go from the same states to the same states and read and
        write the same words.

        Note that the hooks are not checked.

        INPUT:

        - ``left`` -- a transition.

        - ``right`` -- a transition.

        OUTPUT:

        True or False.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState, FSMTransition
            sage: A = FSMState('A', is_initial=True)
            sage: t1 = FSMTransition('A', 'B', 0, 1)
            sage: t2 = FSMTransition(A, 'B', 0, 1)
            sage: t1 == t2
            True
        """
        if not is_FSMTransition(right):
            raise TypeError, 'Only instances of FSMTransition ' \
                'can be compared.'
        return left.from_state == right.from_state \
            and left.to_state == right.to_state \
            and left.word_in == right.word_in \
            and left.word_out == right.word_out


    def __ne__(left, right):
        """

        INPUT:

        - ``left`` -- a transition.

        - ``right`` -- a transition.

        OUTPUT:

        True or False.
        Tests for inequality, complement of __eq__.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState, FSMTransition
            sage: A = FSMState('A', is_initial=True)
            sage: t1 = FSMTransition('A', 'B', 0, 1)
            sage: t2 = FSMTransition(A, 'B', 0, 1)
            sage: t1 != t2
            False
        """
        return (not (left == right))


    def __nonzero__(self):
        """
        Returns True.

        INPUT:

        Nothing.

        OUTPUT:

        True or False.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: FSMTransition('A', 'B', 0).__nonzero__()
            True
        """
        return True  # A transition cannot be zero (see __init__)


#*****************************************************************************


def is_FiniteStateMachine(FSM):
    """
    Tests whether or not ``FSM`` inherits from :class:`FiniteStateMachine`.

    TESTS::

        sage: from sage.combinat.finite_state_machine import is_FiniteStateMachine
        sage: is_FiniteStateMachine(FiniteStateMachine())
        True
        sage: is_FiniteStateMachine(Automaton())
        True
        sage: is_FiniteStateMachine(Transducer())
        True
    """
    return isinstance(FSM, FiniteStateMachine)


class FiniteStateMachine(SageObject):
    """
    Class for a finite state machine.

    A finite state machine is a finite set of states connected by
    transitions.

    INPUT:

    - ``data`` -- can be any of the following:

      #. a dictionary of dictionaries (of transitions),

      #. a dictionary of lists (of states or transitions),

      #. a list (of transitions),

      #. a function (transition function),

      #. an other instance of a finite state machine.

    - ``initial_states`` and ``final_states`` -- the initial and
      final states of this machine

    - ``input_alphabet`` and ``output_alphabet`` -- the input and
      output alphabets of this machine

    - ``determine_alphabets`` -- If True, then the function
      ``determine_alphabets()`` is called after ``data`` was read and
      processed, if ``False``, then not. If it is ``None``, then it is
      decided during the construction of the finite state machine
      whether ``determine_alphabets()`` should be called.

    - ``store_states_dict`` -- If ``True``, then additionally the states
      are stored in an interal dictionary for speed up.

    OUTPUT:

    A finite state machine.

    The object creation of :class:`Automaton` and :class:`Transducer`
    is the same as the one described here (i.e. just replace the word
    ``FiniteStateMachine`` by ``Automaton`` or ``Transducer``).

    Each transition of an automaton has an input label. Automata can,
    for example, be determinised (see
    :meth:`Automaton.determinisation`) and minimized (see
    :meth:`Automaton.minimization`). Each transition of a transducer
    has an input and an output label. Transducers can, for example, be
    simplified (see :meth:`Transducer.simplification`).

    EXAMPLES::

        sage: from sage.combinat.finite_state_machine import FSMState, FSMTransition

    See documentation for more examples.

    We illustrate the different input formats:

    #.  The input-data can be a dictionary of dictionaries, where

        - the keys of the outer dictionary are state-labels (from-states of
          transitions),
        - the keys of the inner dictionaries are state-labels (to-states of
          transitions),
        - the values of the inner dictionaries specify the transition
          more precisely.

        The easiest is to use a tuple consisting of an input and an
        output word::

            sage: FiniteStateMachine({'a':{'b':(0, 1), 'c':(1, 1)}})
            Finite state machine with 3 states

        Instead of the tuple anything iterable (e.g. a list) can be
        used as well.

        If you want to use the arguments of ``FSMTransition``
        directly, you can use a dictionary::

            sage: FiniteStateMachine({'a':{'b':{'word_in':0, 'word_out':1},
            ....:                          'c':{'word_in':1, 'word_out':1}}})
            Finite state machine with 3 states

        In the case you already have instances of ``FSMTransition``, it is
        possible to use them directly::

            sage: FiniteStateMachine({'a':{'b':FSMTransition('a', 'b', 0, 1),
            ....:                          'c':FSMTransition('a', 'c', 1, 1)}})
            Finite state machine with 3 states

    #.  The input-data can be a dictionary of lists, where the keys
        are states or label of states.

        The list-elements can be states::

            sage: a = FSMState('a')
            sage: b = FSMState('b')
            sage: c = FSMState('c')
            sage: FiniteStateMachine({a:[b, c]})
            Finite state machine with 3 states

        Or the list-elements can simply be labels of states::

            sage: FiniteStateMachine({'a':['b', 'c']})
            Finite state machine with 3 states

        The list-elements can also be transitions::

            sage: FiniteStateMachine({'a':[FSMTransition('a', 'b', 0, 1),
            ....:                          FSMTransition('a', 'c', 1, 1)]})
            Finite state machine with 3 states

        Or they can be tuples of a label, an input word and an output
        word specifying a transition::

            sage: FiniteStateMachine({'a':[('b', 0, 1), ('c', 1, 1)]})
            Finite state machine with 3 states

    #.  The input-data can be a list, where its elements specify
        transitions::

            sage: FiniteStateMachine([FSMTransition('a', 'b', 0, 1),
            ....:                     FSMTransition('a', 'c', 1, 1)])
            Finite state machine with 3 states

        It is possible to skip ``FSMTransition`` in the example above::

            sage: FiniteStateMachine([('a', 'b', 0, 1), ('a', 'c', 1, 1)])
            Finite state machine with 3 states

        The parameters of the transition are given in tuples. Anyhow,
        anything iterable (e.g. a list) is possible.

        You can also name the parameters of the transition. For this
        purpose you take a dictionary::

            sage: FiniteStateMachine([{'from_state':'a', 'to_state':'b',
            ....:                      'word_in':0, 'word_out':1},
            ....:                     {'from_state':'a', 'to_state':'c',
            ....:                      'word_in':1, 'word_out':1}])
            Finite state machine with 3 states

        Other arguments, which :class:`FSMTransition` accepts, can be
        added, too.

    #.  The input-data can also be function acting as transition
        function:

        This function has two input arguments:

        #. a label of a state (from which the transition starts),

        #. a letter of the (input-)alphabet (as input-label of the transition).

        It returns a tuple with the following entries:

        #. a label of a state (to which state the transition goes),

        #. a letter of or a word over the (output-)alphabet (as
           output-label of the transition).

        It may also output a list of such tuples if several
        transitions from the from-state and the input letter exist
        (this means that the finite state machine is
        non-deterministic).

        If the transition does not exist, the function should raise a
        ``LookupError`` or return an empty list.

        When constructing a finite state machine in this way, some
        inital states and an input alphabet have to be specified.

        ::

            sage: def f(state_from, read):
            ....:     if int(state_from) + read <= 2:
            ....:         state_to = 2*int(state_from)+read
            ....:         write = 0
            ....:     else:
            ....:         state_to = 2*int(state_from) + read - 5
            ....:         write = 1
            ....:     return (str(state_to), write)
            sage: F = FiniteStateMachine(f, input_alphabet=[0, 1],
            ....:                        initial_states=['0'],
            ....:                        final_states=['0'])
            sage: F([1, 0, 1])
            (True, '0', [0, 0, 1])

    #.  The input-data can be an other instance of a finite state machine::

            sage: FiniteStateMachine(FiniteStateMachine([]))
            Traceback (most recent call last):
            ...
            NotImplementedError


        TESTS::

            sage: a = FSMState('S_a', 'a')
            sage: b = FSMState('S_b', 'b')
            sage: c = FSMState('S_c', 'c')
            sage: d = FSMState('S_d', 'd')
            sage: FiniteStateMachine({a:[b, c], b:[b, c, d],
            ....:                     c:[a, b], d:[a, c]})
            Finite state machine with 4 states

        We have several constructions which lead to the same finite
        state machine::

            sage: A = FSMState('A')
            sage: B = FSMState('B')
            sage: C = FSMState('C')
            sage: FSM1 = FiniteStateMachine(
            ....:  {A:{B:{'word_in':0, 'word_out':1},
            ....:   C:{'word_in':1, 'word_out':1}}})
            sage: FSM2 = FiniteStateMachine({A:{B:(0, 1), C:(1, 1)}})
            sage: FSM3 = FiniteStateMachine(
            ....:  {A:{B:FSMTransition(A, B, 0, 1),
            ....:      C:FSMTransition(A, C, 1, 1)}})
            sage: FSM4 = FiniteStateMachine({A:[(B, 0, 1), (C, 1, 1)]})
            sage: FSM5 = FiniteStateMachine(
            ....:  {A:[FSMTransition(A, B, 0, 1), FSMTransition(A, C, 1, 1)]})
            sage: FSM6 = FiniteStateMachine(
            ....:  [{'from_state':A, 'to_state':B, 'word_in':0, 'word_out':1},
            ....:   {'from_state':A, 'to_state':C, 'word_in':1, 'word_out':1}])
            sage: FSM7 = FiniteStateMachine([(A, B, 0, 1), (A, C, 1, 1)])
            sage: FSM8 = FiniteStateMachine(
            ....:  [FSMTransition(A, B, 0, 1), FSMTransition(A, C, 1, 1)])

            sage: FSM1 == FSM2 == FSM3 == FSM4 == FSM5 == FSM6 == FSM7 == FSM8
            True

        It is possible to skip ``FSMTransition`` in the example above.

        Some more tests for different input-data::

            sage: FiniteStateMachine({'a':{'a':[0, 0], 'b':[1, 1]},
            ....:                     'b':{'b':[1, 0]}})
            Finite state machine with 2 states

            sage: a = FSMState('S_a', 'a')
            sage: b = FSMState('S_b', 'b')
            sage: c = FSMState('S_c', 'c')
            sage: d = FSMState('S_d', 'd')
            sage: t1 = FSMTransition(a, b)
            sage: t2 = FSMTransition(b, c)
            sage: t3 = FSMTransition(b, d)
            sage: t4 = FSMTransition(c, d)
            sage: FiniteStateMachine([t1, t2, t3, t4])
            Finite state machine with 4 states
    """

    #*************************************************************************
    # init
    #*************************************************************************


    def __init__(self,
                 data=None,
                 initial_states=None, final_states=None,
                 input_alphabet=None, output_alphabet=None,
                 determine_alphabets=None,
                 store_states_dict=True):
        """
        See :class:`FiniteStateMachine` for more information.

        TEST::

            sage: FiniteStateMachine()
            Finite state machine with 0 states
        """
        self._states_ = []  # List of states in the finite state
                            # machine.  Each state stores a list of
                            # outgoing transitions.
        if store_states_dict:
            self._states_dict_ = {}

        if initial_states is not None:
            if not hasattr(initial_states, '__iter__'):
                raise TypeError, 'Initial states must be iterable ' \
                    '(e.g. a list of states).'
            for s in initial_states:
                state = self.add_state(s)
                state.is_initial = True

        if final_states is not None:
            if not hasattr(final_states, '__iter__'):
                raise TypeError, 'Final states must be iterable ' \
                    '(e.g. a list of states).'
            for s in final_states:
                state = self.add_state(s)
                state.is_final = True

        self.input_alphabet = input_alphabet
        self.output_alphabet = output_alphabet

        if data is None:
            pass
        elif is_FiniteStateMachine(data):
            raise NotImplementedError
        elif hasattr(data, 'iteritems'):
            # data is a dict (or something similar),
            # format: key = from_state, value = iterator of transitions
            for (sf, iter_transitions) in data.iteritems():
                self.add_state(sf)
                if hasattr(iter_transitions, 'iteritems'):
                    for (st, transition) in iter_transitions.iteritems():
                        self.add_state(st)
                        if is_FSMTransition(transition):
                            self.add_transition(transition)
                        elif hasattr(transition, 'iteritems'):
                            self.add_transition(sf, st, **transition)
                        elif hasattr(transition, '__iter__'):
                            self.add_transition(sf, st, *transition)
                        else:
                            self.add_transition(sf, st, transition)
                elif hasattr(iter_transitions, '__iter__'):
                    for transition in iter_transitions:
                        if hasattr(transition, '__iter__'):
                            L = [sf]
                            L.extend(transition)
                        elif is_FSMTransition(transition):
                            L = transition
                        else:
                            L = [sf, transition]
                        self.add_transition(L)
                else:
                    raise TypeError, 'Wrong input data for transition.'
            if determine_alphabets is None and input_alphabet is None \
                    and output_alphabet is None:
                determine_alphabets = True
        elif hasattr(data, '__iter__'):
            # data is a something that is iterable,
            # items are transitions
            for transition in data:
                if is_FSMTransition(transition):
                    self.add_transition(transition)
                elif hasattr(transition, 'iteritems'):
                    self.add_transition(transition)
                elif hasattr(transition, '__iter__'):
                    self.add_transition(transition)
                else:
                    raise TypeError, 'Wrong input data for transition.'
            if determine_alphabets is None and input_alphabet is None \
                    and output_alphabet is None:
                determine_alphabets = True
        elif hasattr(data, '__call__'):
            self.add_from_transition_function(data)
        else:
            raise TypeError, 'Cannot decide what to do with data.'

        if determine_alphabets:
            self.determine_alphabets()


    #*************************************************************************
    # copy and hash
    #*************************************************************************


    def __copy__(self):
        """
        Returns a (shallow) copy of the finite state machine.

        INPUT:

        Nothing.

        OUTPUT:

        A new finite state machine.

        TESTS::

            sage: copy(FiniteStateMachine())
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


    copy = __copy__

    def empty_copy(self, memo=None):
        """
        Returns an empty deep copy of the finite state machine, i.e.,
        input_alphabet, output_alphabet are preserved, but states and
        transitions are not.

        INPUT:

        - ``memo`` -- a dictionary storing already processed elements.

        OUTPUT:

        A new finite state machine.

        EXAMPLES::

            sage: F = FiniteStateMachine([('A', 'A', 0, 2), ('A', 'A', 1, 3)],
            ....:                        input_alphabet=[0,1],
            ....:                        output_alphabet=[2,3])
            sage: FE = F.empty_copy(); FE
            Finite state machine with 0 states
            sage: FE.input_alphabet
            [0, 1]
            sage: FE.output_alphabet
            [2, 3]
        """
        new = self.__class__()
        new.input_alphabet = deepcopy(self.input_alphabet, memo)
        new.output_alphabet = deepcopy(self.output_alphabet, memo)
        return new

    def __deepcopy__(self, memo):
        """
        Returns a deep copy of the finite state machine.

        INPUT:

        - ``memo`` -- a dictionary storing already processed elements.

        OUTPUT:

        A new finite state machine.

        EXAMPLES::

            sage: F = FiniteStateMachine([('A', 'A', 0, 1), ('A', 'A', 1, 0)])
            sage: deepcopy(F)
            Finite state machine with 1 states
        """
        relabel = hasattr(self, '_deepcopy_relabel_')
        new = self.empty_copy(memo=memo)
        relabel_iter = itertools.count(0)
        for state in self.iter_states():
            if relabel:
                state._deepcopy_relabel_ = relabel_iter.next()
            s = deepcopy(state, memo)
            if relabel:
                del state._deepcopy_relabel_
            new.add_state(s)
        for transition in self.iter_transitions():
            new.add_transition(deepcopy(transition, memo))
        return new


    def deepcopy(self, memo=None):
        """
        Returns a deep copy of the finite state machine.

        INPUT:

        - ``memo`` -- (default: ``None``) a dictionary storing already
          processed elements.

        OUTPUT:

        A new finite state machine.

        EXAMPLES::

            sage: F = FiniteStateMachine([('A', 'A', 0, 1), ('A', 'A', 1, 0)])
            sage: deepcopy(F)
            Finite state machine with 1 states
        """
        return deepcopy(self, memo)


    def relabeled(self, memo=None):
        """
        Returns a deep copy of the finite state machine, but the
        states are relabeled by integers starting with 0.

        INPUT:

        - ``memo`` -- (default: ``None``) a dictionary storing already
          processed elements.

        OUTPUT:

        A new finite state machine.

        EXAMPLES::

            sage: FSM1 = FiniteStateMachine([('A', 'B'), ('B', 'C'), ('C', 'A')])
            sage: FSM1.states()
            ['A', 'B', 'C']
            sage: FSM2 = FSM1.relabeled()
            sage: FSM2.states()
            [0, 1, 2]
        """
        self._deepcopy_relabel_ = True
        new = deepcopy(self, memo)
        del self._deepcopy_relabel_
        return new


    def __hash__(self):
        """
        Since finite state machines are mutable, they should not be
        hashable, so we return a type error.

        INPUT:

        Nothing.

        OUTPUT:

        The hash of this finite state machine.

        EXAMPLES::

            sage: hash(FiniteStateMachine())
            Traceback (most recent call last):
            ...
            TypeError: Finite state machines are mutable, and thus not hashable.
        """
        if getattr(self, "_immutable", False):
            return hash((tuple(self.states()), tuple(self.transitions())))
        raise TypeError, "Finite state machines are mutable, " \
            "and thus not hashable."


    #*************************************************************************
    # operators
    #*************************************************************************


    def __add__(self, other):
        """
        Returns the disjoint union of the finite state machines self and other.

        INPUT:

        - ``other`` -- a finite state machine.

        OUTPUT:

        A new finite state machine.

        TESTS::

            sage: FiniteStateMachine() + FiniteStateMachine([('A', 'B')])
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if is_FiniteStateMachine(other):
            return self.disjoint_union(other)


    def __iadd__(self, other):
        """
        TESTS::

            sage: F = FiniteStateMachine()
            sage: F += FiniteStateMachine()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


    def __mul__(self, other):
        """
        TESTS::

            sage: FiniteStateMachine() * FiniteStateMachine([('A', 'B')])
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if is_FiniteStateMachine(other):
            return self.intersection(other)


    def __imul__(self, other):
        """
        TESTS::

            sage: F = FiniteStateMachine()
            sage: F *= FiniteStateMachine()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


    def __call__(self, *args, **kwargs):
        """
        Calls either method :meth:`.composition` or :meth:`.process`.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A', is_initial=True, is_final=True)
            sage: binary_inverter = Transducer({A:[(A, 0, 1), (A, 1, 0)]})
            sage: binary_inverter([0, 1, 0, 0, 1, 1])
            (True, 'A', [1, 0, 1, 1, 0, 0])

        ::

            sage: F = Transducer([('A', 'B', 1, 0), ('B', 'B', 1, 1),
            ....:                 ('B', 'B', 0, 0)],
            ....:                initial_states=['A'], final_states=['B'])
            sage: G = Transducer([(1, 1, 0, 0), (1, 2, 1, 0),
            ....:                 (2, 2, 0, 1), (2, 1, 1, 1)],
            ....:                initial_states=[1], final_states=[1])
            sage: H = G(F)
            sage: H.states()
            [('A', 1), ('B', 1), ('B', 2)]
        """
        if len(args) == 0:
            raise TypeError, "Called with too few arguments."
        if is_FiniteStateMachine(args[0]):
            return self.composition(*args, **kwargs)
        if hasattr(args[0], '__iter__'):
            return self.process(*args, **kwargs)
        raise TypeError, "Do not know what to do with that arguments."


    #*************************************************************************
    # tests
    #*************************************************************************


    def __nonzero__(self):
        """
        Returns True if the finite state machine consists of at least
        one state.

        INPUT:

        Nothing.

        OUTPUT:

        True or False.

        TESTS::

            sage: FiniteStateMachine().__nonzero__()
            False
        """
        return len(self._states_) > 0


    def __eq__(left, right):
        """
        Returns True if the two finite state machines are equal, i.e.,
        if they have the same states and the same transitions.

        INPUT:

        - ``left`` -- a finite state machine.

        - ``right`` -- a finite state machine.

        OUTPUT:

        True or False.

        EXAMPLES::

            sage: F = FiniteStateMachine([('A', 'B', 1)])
            sage: F == FiniteStateMachine()
            False
        """
        if not is_FiniteStateMachine(right):
            raise TypeError, 'Only instances of FiniteStateMachine ' \
                'can be compared.'
        if len(left._states_) != len(right._states_):
            return False
        for state in left.iter_states():
            if state not in right._states_:
                return False
            left_transitions = state.transitions
            right_transitions = right.state(state).transitions
            if len(left_transitions) != len(right_transitions):
                return False
            for t in left_transitions:
                if t not in right_transitions:
                    return False
        return True


    def __ne__(left, right):
        """
        Tests for inequality, complement of :meth:`.__eq__`.

        INPUT:

        - ``left`` -- a finite state machine.

        - ``right`` -- a finite state machine.

        OUTPUT:

        True or False.

        EXAMPLES::

            sage: E = FiniteStateMachine([('A', 'B', 0)])
            sage: F = Automaton([('A', 'B', 0)])
            sage: G = Transducer([('A', 'B', 0, 1)])
            sage: E == F
            True
            sage: E == G
            False
        """
        return (not (left == right))


    def __contains__(self, item):
        """
        Returns true, if the finite state machine contains the
        state or transition item. Note that only the labels of the
        states and the input and output words are tested.

        INPUT:

        - ``item`` -- a state or a transition.

        OUTPUT:

        True or False.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState, FSMTransition
            sage: F = FiniteStateMachine([('A', 'B', 0), ('B', 'A', 1)])
            sage: FSMState('A', is_initial=True) in F
            True
            sage: 'A' in F
            False
            sage: FSMTransition('A', 'B', 0) in F
            True
        """
        if is_FSMState(item):
            return self.has_state(item)
        if is_FSMTransition(item):
            return self.has_transition(item)
        return False


    #*************************************************************************
    # representations / LaTeX
    #*************************************************************************


    def _repr_(self):
        """
        Represents the finite state machine as "Finite state machine
        with n states" where n is the number of states.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: FiniteStateMachine()._repr_()
            'Finite state machine with 0 states'
        """
        return "Finite state machine with %s states" % len(self._states_)


    def _latex_(self):
        r"""
        Returns a LaTeX code for the graph of the finite state machine.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: F = FiniteStateMachine([('A', 'B', 1, 2)])
            sage: F._latex_()
            '\\begin{tikzpicture}[auto]\n\\node[state] (v0) at (3.000000,0.000000) {\\text{\\texttt{A}}}\n;\\node[state] (v1) at (-3.000000,0.000000) {\\text{\\texttt{B}}}\n;\\path[->] (v0) edge node {$ $} (v1);\n\\end{tikzpicture}'
        """
        result = "\\begin{tikzpicture}[auto]\n"
        j = 0;
        for vertex in self.states():
            if not hasattr(vertex, "coordinates"):
                vertex.coordinates = (3*cos(2*pi*j/len(self.states())),
                                      3*sin(2*pi*j/len(self.states())))
            options = ""
            if vertex in self.final_states():
                options += ",accepting"
            if hasattr(vertex, "format_label"):
                label = vertex.format_label()
            elif hasattr(self, "format_state_label"):
                label = self.format_state_label(vertex)
            else:
                label = latex(vertex.label())
            result += "\\node[state%s] (v%d) at (%f,%f) {%s}\n;" % (
                options, j, vertex.coordinates[0],
                vertex.coordinates[1], label)
            vertex._number_ = j
            j += 1
        adjacent = {}
        for source in self.states():
            for target in self.states():
                transitions = filter(lambda transition: \
                                         transition.to_state == target,
                                     source.transitions)
                adjacent[source, target] = transitions

        for ((source, target), transitions) in adjacent.iteritems():
            if len(transitions) > 0:
                labels = []
                for transition in transitions:
                    if hasattr(transition, "format_label"):
                        labels.append(transition.format_label())
                        continue
                    elif hasattr(self, "format_transition_label"):
                        format_transition_label = self.format_transition_label
                    else:
                        format_transition_label = latex
                    labels.append(self._latex_transition_label_(
                            transition, format_transition_label))
                label = ", ".join(labels)
                if source != target:
                    if len(adjacent[target, source]) > 0:
                        angle = atan2(
                            target.coordinates[1] - source.coordinates[1],
                            target.coordinates[0]-source.coordinates[0])*180/pi
                        angle_source = ".%.2f" % ((angle+5).n(),)
                        angle_target = ".%.2f" % ((angle+175).n(),)
                    else:
                        angle_source = ""
                        angle_target = ""
                    result += "\\path[->] (v%d%s) edge node {$%s$} (v%d%s);\n" % (
                        source._number_, angle_source, label,
                        target._number_, angle_target)
                else:
                    result += "\\path[->] (v%d) edge[loop above] node {$%s$} ();\n" % (
                        source._number_, label)

        result += "\\end{tikzpicture}"
        return result


    def _latex_transition_label_(self, transition, format_function=latex):
        r"""
        Returns the proper transition label.

        INPUT:

        - ``transition`` - a transition

        - ``format_function'' - a function formatting the labels

        OUTPUT:

        A string.

        TESTS::

            sage: F = FiniteStateMachine([('A', 'B', 0, 1)])
            sage: t = F.transitions()[0]
            sage: F._latex_transition_label_(t)
            ' '
        """
        return ' '


    #*************************************************************************
    # other
    #*************************************************************************


    def _matrix_(self, R=None):
        """
        Returns the adjacency matrix of the finite state machine.
        See :meth:`.adjacency_matrix` for more information.

        EXAMPLES::

            sage: B = FiniteStateMachine({0: {0: (0, 0), 'a': (1, 0)},
            ....:                         'a': {2: (0, 0), 3: (1, 0)},
            ....:                         2:{0:(1, 1), 4:(0, 0)},
            ....:                         3:{'a':(0, 1), 2:(1, 1)},
            ....:                         4:{4:(1, 1), 3:(0, 1)}},
            ....:                        initial_states=[0])
            sage: B._matrix_()
            [1 1 0 0 0]
            [0 0 1 1 0]
            [x 0 0 0 1]
            [0 x x 0 0]
            [0 0 0 x x]
        """
        return self.adjacency_matrix()


    def adjacency_matrix(self, input=None,
                         entry=(lambda transition:var('x')**transition.word_out[0])):
        """
        Returns the adjacency matrix of the underlying graph.

        INPUT:

        - ``input`` -- Only transitions with input label ``input`` are
          respected.

        - ``entry`` -- The function ``entry`` takes a transition and
          the return value is written in the matrix as the entry
          ``(transition.from_state, transition.to_state)``.

        OUTPUT:

        A matrix.

        If any label of a state is not an integer, the finite state
        machine is relabeled at the beginning.  If there are more than
        one transitions between two states, then the different return
        values of ``entry`` are added up.

        The default value of entry takes the variable ``x`` to the
        power of the output word of the transition.

        EXAMPLES::

            sage: B = FiniteStateMachine({0:{0:(0, 0), 'a':(1, 0)},
            ....:                         'a':{2:(0, 0), 3:(1, 0)},
            ....:                         2:{0:(1, 1), 4:(0, 0)},
            ....:                         3:{'a':(0, 1), 2:(1, 1)},
            ....:                         4:{4:(1, 1), 3:(0, 1)}},
            ....:                        initial_states=[0])
            sage: B.adjacency_matrix()
            [1 1 0 0 0]
            [0 0 1 1 0]
            [x 0 0 0 1]
            [0 x x 0 0]
            [0 0 0 x x]
            sage: B.adjacency_matrix(entry=(lambda transition: 1))
            [1 1 0 0 0]
            [0 0 1 1 0]
            [1 0 0 0 1]
            [0 1 1 0 0]
            [0 0 0 1 1]
            sage: B.adjacency_matrix(1, entry=(lambda transition:
            ....:     exp(I*transition.word_out[0]*var('t'))))
            [      0       1       0       0       0]
            [      0       0       0       1       0]
            [e^(I*t)       0       0       0       0]
            [      0       0 e^(I*t)       0       0]
            [      0       0       0       0 e^(I*t)]

        """
        relabeledFSM = self
        l = len(relabeledFSM.states())
        for state in self.states():
            if state.label() not in ZZ or state.label() >= l or state.label() < 0:
                relabeledFSM = self.relabeled()
                break
        dictionary = {}
        for transition in relabeledFSM.iter_transitions():
            if input is None or transition.word_in == [input]:
                if (transition.from_state.label(), transition.to_state.label()) in dictionary:
                    dictionary[(transition.from_state.label(), transition.to_state.label())] += entry(transition)
                else:
                    dictionary[(transition.from_state.label(), transition.to_state.label())] = entry(transition)
        return matrix(len(relabeledFSM.states()), dictionary)


    def determine_alphabets(self, reset=True):
        """
        Determines the input and output alphabet according to the
        transitions in self.

        INPUT:

        - ``reset`` -- If reset is True, then the existing input
          alphabet is erased, otherwise new letters are appended to
          the existing alphabet.

        OUTPUT:

        Nothing.

        After this operation the input alphabet and the output
        alphabet of self are a list of letters.

        EXAMPLES::

            sage: T = Transducer([(1, 1, 1, 0), (1, 2, 2, 1),
            ....:                 (2, 2, 1, 1), (2, 2, 0, 0)],
            ....:                determine_alphabets=False)
            sage: (T.input_alphabet, T.output_alphabet)
            (None, None)
            sage: T.determine_alphabets()
            sage: (T.input_alphabet, T.output_alphabet)
            ([0, 1, 2], [0, 1])
       """
        if reset:
            ain = set()
            aout = set()
        else:
            ain = set(self.input_alphabet)
            aout = set(self.output_alphabet)

        for t in self.iter_transitions():
            for letter in t.word_in:
                ain.add(letter)
            for letter in t.word_out:
                aout.add(letter)
        self.input_alphabet = list(ain)
        self.output_alphabet = list(aout)


    #*************************************************************************
    # get states and transitions
    #*************************************************************************


    def states(self):
        """
        Returns the states of the finite state machine.

        INPUT:

        Nothing.

        OUTPUT:

        The states of the finite state machine as list.

        EXAMPLES::

            sage: FSM = Automaton([('1', '2', 1), ('2', '2', 0)])
            sage: FSM.states()
            ['1', '2']

        """
        return copy(self._states_)


    def iter_states(self):
        """
        Returns an iterator of the states.

        INPUT:

        Nothing.

        OUTPUT:

        An iterator of the states of the finite state machine.

        EXAMPLES::

            sage: FSM = Automaton([('1', '2', 1), ('2', '2', 0)])
            sage: [s.label() for s in FSM.iter_states()]
            ['1', '2']
        """
        return iter(self._states_)


    def transitions(self, from_state=None):
        """
        Returns a list of all transitions.

        INPUT:

        - ``from_state`` -- (default: ``None``) If ``from_state`` is
          given, then a list of transitions starting there is given.

        OUTPUT:

        A list of all transitions.

        EXAMPLES::

            sage: FSM = Automaton([('1', '2', 1), ('2', '2', 0)])
            sage: FSM.transitions()
            [Transition from '1' to '2': 1|-,
             Transition from '2' to '2': 0|-]
        """
        return list(self.iter_transitions(from_state))


    def iter_transitions(self, from_state=None):
        """
        Returns an iterator of all transitions.

        INPUT:

        - ``from_state`` -- (default: ``None``) If ``from_state`` is
          given, then a list of transitions starting there is given.

        OUTPUT:

        An iterator of all transitions.

        EXAMPLES::

            sage: FSM = Automaton([('1', '2', 1), ('2', '2', 0)])
            sage: [(t.from_state.label(), t.to_state.label())
            ....:     for t in FSM.iter_transitions('1')]
            [('1', '2')]
            sage: [(t.from_state.label(), t.to_state.label())
            ....:     for t in FSM.iter_transitions('2')]
            [('2', '2')]
            sage: [(t.from_state.label(), t.to_state.label())
            ....:     for t in FSM.iter_transitions()]
            [('1', '2'), ('2', '2')]
        """
        if from_state is None:
            return self._iter_transitions_all_()
        else:
            return iter(self.state(from_state).transitions)


    def _iter_transitions_all_(self):
        """
        Returns an iterator over all transitions.

        INPUT:

        Nothing.

        OUTPUT:

        An iterator over all transitions.

        EXAMPLES::

            sage: FSM = Automaton([('1', '2', 1), ('2', '2', 0)])
            sage: [(t.from_state.label(), t.to_state.label())
            ....:     for t in FSM._iter_transitions_all_()]
            [('1', '2'), ('2', '2')]
        """
        for state in self.iter_states():
            for t in state.transitions:
                yield t


    def initial_states(self):
        """
        Returns a list of all initial states.

        INPUT:

        Nothing.

        OUTPUT:

        A list of all initial states.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A', is_initial=True)
            sage: B = FSMState('B')
            sage: F = FiniteStateMachine([(A, B, 1, 0)])
            sage: F.initial_states()
            ['A']
        """
        return list(self.iter_initial_states())


    def iter_initial_states(self):
        """
        Returns an iterator of the initial states.

        INPUT:

        Nothing.

        OUTPUT:

        An iterator over all initial states.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A', is_initial=True)
            sage: B = FSMState('B')
            sage: F = FiniteStateMachine([(A, B, 1, 0)])
            sage: [s.label() for s in F.iter_initial_states()]
            ['A']
        """
        return itertools.ifilter(lambda s:s.is_initial, self.iter_states())


    def final_states(self):
        """
        Returns a list of all final states.

        INPUT:

        Nothing.

        OUTPUT:

        A list of all final states.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A', is_final=True)
            sage: B = FSMState('B', is_initial=True)
            sage: C = FSMState('C', is_final=True)
            sage: F = FiniteStateMachine([(A, B), (A, C)])
            sage: F.final_states()
            ['A', 'C']
        """
        return list(self.iter_final_states())


    def iter_final_states(self):
        """
        Returns an iterator of the final states.

        INPUT:

        Nothing.

        OUTPUT:

        An iterator over all initial states.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A', is_final=True)
            sage: B = FSMState('B', is_initial=True)
            sage: C = FSMState('C', is_final=True)
            sage: F = FiniteStateMachine([(A, B), (A, C)])
            sage: [s.label() for s in F.iter_final_states()]
            ['A', 'C']
        """
        return itertools.ifilter(lambda s:s.is_final, self.iter_states())


    def state(self, state):
        """
        Returns the state of the finite state machine.

        INPUT:

        - ``state`` -- If ``state`` is not an instance of
          :class:`FSMState`, then it is assumed that it is the label
          of a state.

        OUTPUT:

        Returns the state of the finite state machine corresponding to
        ``state``.

        If no state is found, then a ``LookupError`` is thrown.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A')
            sage: FSM = FiniteStateMachine([(A, 'B'), ('C', A)])
            sage: FSM.state('A') == A
            True
            sage: FSM.state('xyz')
            Traceback (most recent call last):
            ...
            LookupError: No state with label xyz found.
        """
        def what(s, switch):
            if switch:
                return s.label()
            else:
                return s
        switch = is_FSMState(state)

        try:
            return self._states_dict_[what(state, switch)]
        except AttributeError:
            for s in self.iter_states():
                if what(s, not switch) == state:
                    return s
        except KeyError:
            pass
        raise LookupError, \
            "No state with label %s found." % (what(state, switch),)


    def transition(self, transition):
        """
        Returns the transition of the finite state machine.

        INPUT:

        - ``transition`` -- If ``transition`` is not an instance of
          :class:`FSMTransition`, then it is assumed that it is a
          tuple ``(from_state, to_state, word_in, word_out)``.

        OUTPUT:

        Returns the transition of the finite state machine
        corresponding to ``transition``.

        If no transition is found, then a ``LookupError`` is thrown.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: t = FSMTransition('A', 'B', 0)
            sage: F = FiniteStateMachine([t])
            sage: F.transition(('A', 'B', 0))
            Transition from 'A' to 'B': 0|-
            sage: id(t) == id(F.transition(('A', 'B', 0)))
            True
       """
        if not is_FSMTransition(transition):
            transition = FSMTransition(*transition)
        for s in self.iter_transitions(transition.from_state):
            if s == transition:
                return s
        raise LookupError, "No transition found."


    #*************************************************************************
    # properties (state and transitions)
    #*************************************************************************


    def has_state(self, state):
        """
        Returns whether ``state`` is one of the states of the finite
        state machine.

        INPUT:

        - ``state`` can be a :class:`FSMState` or a label of a state.

        OUTPUT:

        True or False.

        EXAMPLES::

            sage: FiniteStateMachine().has_state('A')
            False
        """
        try:
            self.state(state)
            return True
        except LookupError:
            return False


    def has_transition(self, transition):
        """
        Returns whether ``transition`` is one of the transitions of
        the finite state machine.

        INPUT:

        - ``transition`` has to be a :class:`FSMTransition`.

        OUTPUT:

        True or False.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: t = FSMTransition('A', 'A', 0, 1)
            sage: FiniteStateMachine().has_transition(t)
            False
            sage: FiniteStateMachine().has_transition(('A', 'A', 0, 1))
            Traceback (most recent call last):
            ...
            TypeError: Transition is not an instance of FSMTransition.
        """
        if is_FSMTransition(transition):
            return transition in self.iter_transitions()
        raise TypeError, "Transition is not an instance of FSMTransition."


    def has_initial_state(self, state):
        """
        Returns whether ``state`` is one of the initial states of the
        finite state machine.

        INPUT:

        - ``state`` can be a :class:`FSMState` or a label.

        OUTPUT:

        True or False.

        EXAMPLES::

            sage: F = FiniteStateMachine([('A', 'A')], initial_states=['A'])
            sage: F.has_initial_state('A')
            True
        """
        try:
            return self.state(state).is_initial
        except LookupError:
            return False


    def has_initial_states(self):
        """
        Returns whether the finite state machine has an initial state.

        INPUT:

        Nothing.

        OUTPUT:

        True or False.

        EXAMPLES::

            sage: FiniteStateMachine().has_initial_states()
            False
        """
        return len(self.initial_states()) > 0


    def has_final_state(self, state):
        """
        Returns whether ``state`` is one of the final states of the
        finite state machine.

        INPUT:

        - ``state`` can be a :class:`FSMState` or a label.

        OUTPUT:

        True or False.

        EXAMPLES::

            sage: FiniteStateMachine(final_states=['A']).has_final_state('A')
            True
        """
        try:
            return self.state(state).is_final
        except LookupError:
            return False


    def has_final_states(self):
        """
        Returns whether the finite state machine has a final state.

        INPUT:

        Nothing.

        OUTPUT:

        True or False.

        EXAMPLES::

            sage: FiniteStateMachine().has_final_states()
            False
        """
        return len(self.final_states()) > 0


    #*************************************************************************
    # properties
    #*************************************************************************


    def is_deterministic(self):
        """
        Returns whether the finite finite state machine is deterministic.

        INPUT:

        Nothing.

        OUTPUT:

        True or False.

        A finite state machine is considered to be deterministic if
        each transition has input label of length one and for each
        pair `(q,a)` where `q` is a state and `a` is an element of the
        input alphabet, there is at most one transition from `q` with
        input label `a`.

        TESTS::

            sage: fsm = FiniteStateMachine()
            sage: fsm.add_transition(('A', 'B', 0, []))
            Transition from 'A' to 'B': 0|-
            sage: fsm.is_deterministic()
            True
            sage: fsm.add_transition(('A', 'C', 0, []))
            Transition from 'A' to 'C': 0|-
            sage: fsm.is_deterministic()
            False
            sage: fsm.add_transition(('A', 'B', [0,1], []))
            Transition from 'A' to 'B': 0,1|-
            sage: fsm.is_deterministic()
            False
        """
        for state in self.states():
            for transition in state.transitions:
                if len(transition.word_in) != 1:
                    return False

            transition_classes_by_word_in = full_group_by(
                state.transitions,
                key=lambda t:t.word_in)

            for key,transition_class in transition_classes_by_word_in:
                if len(transition_class) > 1:
                    return False
        return True


    def is_connected(self):
        """
        TESTS::

            sage: FiniteStateMachine().is_connected()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


    #*************************************************************************
    # let the finite state machine work
    #*************************************************************************


    def process(self, *args, **kwargs):
        """
        Returns whether the finite state machine accepts the input, the state
        where the computation stops and which output is generated.

        INPUT:

        - ``input_tape`` -- The input tape can be a list with entries from
          the input alphabet.

        - ``initial_state`` -- (default: ``None``) The state in which
          to start. If this parameter is ``None`` and there is only
          one initial state in the machine, then this state is taken.

        OUTPUT:

        A triple, where

        - the first entry is true if the input string is accepted,

        - the second gives the reached state after processing the
          input tape, and

        - the third gives a list of the output labels used during
          processing (in the case the finite state machine runs as
          transducer).

        Note that in the case the finite state machine is not
        deterministic, one possible path is gone. This means, that in
        this case the output can be wrong. Use
        :meth:`.determinisation` to get a deterministic finite state
        machine and try again.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A', is_initial = True, is_final = True)
            sage: binary_inverter = Transducer({A:[(A, 0, 1), (A, 1, 0)]})
            sage: binary_inverter.process([0, 1, 0, 0, 1, 1])
            (True, 'A', [1, 0, 1, 1, 0, 0])

        Alternatively, we can invoke this function by::

            sage: binary_inverter([0, 1, 0, 0, 1, 1])
            (True, 'A', [1, 0, 1, 1, 0, 0])

        ::

            sage: NAF_ = FSMState('_', is_initial = True, is_final = True)
            sage: NAF1 = FSMState('1', is_final = True)
            sage: NAF = Automaton(
            ....:     {NAF_: [(NAF_, 0), (NAF1, 1)], NAF1: [(NAF_, 0)]})
            sage: [NAF.process(w)[0] for w in [[0], [0, 1], [1, 1], [0, 1, 0, 1],
            ....: [0, 1, 1, 1, 0], [1, 0, 0, 1, 1]]]
            [True, True, False, True, False, False]

        """
        it = self.iter_process(*args, **kwargs)
        for _ in it:
            pass
        return (it.accept_input, it.current_state, it.output_tape)


    def iter_process(self, input_tape=None, initial_state=None):
        """
        See `process` for more informations.

        EXAMPLES::

            sage: inverter = Transducer({'A': [('A', 0, 1), ('A', 1, 0)]},
            ....:     initial_states=['A'], final_states=['A'])
            sage: it = inverter.iter_process(input_tape=[0, 1, 1])
            sage: for _ in it:
            ....:     pass
            sage: it.output_tape
            [1, 0, 0]
        """
        return FSMProcessIterator(self, input_tape, initial_state)


    #*************************************************************************
    # change finite state machine (add/remove state/transitions)
    #*************************************************************************


    def add_state(self, state):
        """
        Adds a state to the finite state machine and returns the new
        state. If the state already exists, that existing state is
        returned.

        INPUT:

        - ``state`` is either an instance of FSMState or, otherwise, a
          label of a state.

        OUTPUT:

        The new or existing state.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: F = FiniteStateMachine()
            sage: A = FSMState('A', is_initial=True)
            sage: F.add_state(A)
            'A'
        """
        try:
            return self.state(state)
        except LookupError:
            pass
        # at this point we know that we have a new state
        if is_FSMState(state):
            s = state
        else:
            s = FSMState(state)
        s.transitions = list()
        self._states_.append(s)
        try:
            self._states_dict_[s.label()] = s
        except AttributeError:
            pass
        return s


    def add_states(self, states):
        """
        Adds several states. See add_state for more information.

        INPUT:

        - ``states`` -- a list of states or iterator over states.

        OUTPUT:

        Nothing.

        EXAMPLES::

            sage: F = FiniteStateMachine()
            sage: F.add_states(['A', 'B'])
            sage: F.states()
            ['A', 'B']
        """
        for state in states:
            self.add_state(state)


    def add_transition(self, *args, **kwargs):
        """
        Adds a transition to the finite state machine and returns the
        new transition. If the transition already exists, that
        existing transition is returned.

        INPUT:

        The following forms are all accepted:

        ::

            sage: from sage.combinat.finite_state_machine import FSMState, FSMTransition
            sage: A = FSMState('A')
            sage: B = FSMState('B')

            sage: FSM = FiniteStateMachine()
            sage: FSM.add_transition(FSMTransition(A, B, 0, 1))
            Transition from 'A' to 'B': 0|1

            sage: FSM = FiniteStateMachine()
            sage: FSM.add_transition(A, B, 0, 1)
            Transition from 'A' to 'B': 0|1

            sage: FSM = FiniteStateMachine()
            sage: FSM.add_transition(A, B, word_in=0, word_out=1)
            Transition from 'A' to 'B': 0|1

            sage: FSM = FiniteStateMachine()
            sage: FSM.add_transition('A', 'B', {'word_in': 0, 'word_out': 1})
            Transition from 'A' to 'B': {'word_in': 0, 'word_out': 1}|-

            sage: FSM = FiniteStateMachine()
            sage: FSM.add_transition(from_state=A, to_state=B,
            ....:                    word_in=0, word_out=1)
            Transition from 'A' to 'B': 0|1

            sage: FSM = FiniteStateMachine()
            sage: FSM.add_transition({'from_state': A, 'to_state': B,
            ....:                    'word_in': 0, 'word_out': 1})
            Transition from 'A' to 'B': 0|1

            sage: FSM = FiniteStateMachine()
            sage: FSM.add_transition((A, B, 0, 1))
            Transition from 'A' to 'B': 0|1

            sage: FSM = FiniteStateMachine()
            sage: FSM.add_transition([A, B, 0, 1])
            Transition from 'A' to 'B': 0|1

        If the states ``A`` and ``B`` are not instances of FSMState, then
        it is assumed that they are labels of states.

        OUTPUT:

        The new or existing transition.
        """
        if len(args) + len(kwargs) == 0:
            return
        if len(args) + len(kwargs) == 1:
            if len(args) == 1:
                d = args[0]
                if is_FSMTransition(d):
                    return self._add_fsm_transition_(d)
            else:
                d = kwargs.itervalues().next()
            if hasattr(d, 'iteritems'):
                args = []
                kwargs = d
            elif hasattr(d, '__iter__'):
                args = d
                kwargs = {}
            else:
                raise TypeError, "Cannot decide what to do with input."

        data = dict(zip(
                ('from_state', 'to_state', 'word_in', 'word_out', 'hook'),
                args))
        data.update(kwargs)

        data['from_state'] = self.add_state(data['from_state'])
        data['to_state'] = self.add_state(data['to_state'])

        return self._add_fsm_transition_(FSMTransition(**data))


    def _add_fsm_transition_(self, t):
        """
        Adds a transition.

        INPUT:

        - ``t`` -- an instance of ``FSMTransition``.

        OUTPUT:

        The new transition.

        TESTS::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: F = FiniteStateMachine()
            sage: F._add_fsm_transition_(FSMTransition('A', 'B'))
            Transition from 'A' to 'B': -|-
        """
        try:
            return self.transition(t)
        except LookupError:
            pass
        from_state = self.add_state(t.from_state)
        self.add_state(t.to_state)
        from_state.transitions.append(t)
        return t


    def add_from_transition_function(self, function, initial_states=None,
                                     explore_existing_states=True):
        """
        Constructs a finite state machine from a transition function.

        INPUT:

        - ``function`` may return a tuple (new_state, output_word) or a
          list of such tuples.

        - ``initial_states`` -- If no initial states are given, the
          already existing initial states of self are taken.

        - If ``explore_existing_states`` is True (default), then
          already existing states in self (e.g. already given final
          states) will also be processed if they are reachable from
          the initial states.

        OUTPUT:

        Nothing.

        EXAMPLES::

            sage: F = FiniteStateMachine(initial_states=['A'],
            ....:                        input_alphabet=[0, 1])
            sage: def f(state, input):
            ....:     return [('A', input), ('B', 1-input)]
            sage: F.add_from_transition_function(f)
            sage: F.transitions()
            [Transition from 'A' to 'A': 0|0,
            Transition from 'A' to 'B': 0|1,
            Transition from 'A' to 'A': 1|1,
            Transition from 'A' to 'B': 1|0,
            Transition from 'B' to 'A': 0|0,
            Transition from 'B' to 'B': 0|1,
            Transition from 'B' to 'A': 1|1,
            Transition from 'B' to 'B': 1|0]

        Initial states can also be given as a parameter::

            sage: F = FiniteStateMachine(input_alphabet=[0,1])
            sage: def f(state, input):
            ....:     return [('A', input), ('B', 1-input)]
            sage: F.add_from_transition_function(f,initial_states=['A'])
            sage: F.initial_states()
            ['A']

        Already existing states in the finite state machine (the final
        states in the example below) are also explored::

            sage: F = FiniteStateMachine(initial_states=[0],
            ....:                        final_states=[1],
            ....:                        input_alphabet=[0])
            sage: def transition_function(state, letter):
            ....:     return(1-state, [])
            sage: F.add_from_transition_function(transition_function)
            sage: F.transitions()
            [Transition from 0 to 1: 0|-,
             Transition from 1 to 0: 0|-]

        If ``explore_existing_states=False``, however, this behavior
        is turned off, i.e., already existing states are not
        explored::

            sage: F = FiniteStateMachine(initial_states=[0],
            ....:                        final_states=[1],
            ....:                        input_alphabet=[0])
            sage: def transition_function(state, letter):
            ....:     return(1-state, [])
            sage: F.add_from_transition_function(transition_function,
            ....:                                explore_existing_states=False)
            sage: F.transitions()
            [Transition from 0 to 1: 0|-]

        TEST::

            sage: F = FiniteStateMachine(initial_states=['A'])
            sage: def f(state, input):
            ....:     return [('A', input), ('B', 1-input)]
            sage: F.add_from_transition_function(f)
            Traceback (most recent call last):
            ...
            ValueError: No input alphabet is given.
            Try calling determine_alphabets().
        """
        if self.input_alphabet is None:
            raise ValueError, ("No input alphabet is given. "
                               "Try calling determine_alphabets().")

        if initial_states is None:
            not_done = self.initial_states()
        elif hasattr(initial_states, '__iter__'):
            not_done = []
            for s in initial_states:
                state = self.add_state(s)
                state.is_initial = True
                not_done.append(state)
        else:
            raise TypeError, 'Initial states must be iterable ' \
                '(e.g. a list of states).'
        if len(not_done) == 0:
            raise ValueError, "No state is initial."
        if explore_existing_states:
            ignore_done = self.states()
            for s in not_done:
                try:
                    ignore_done.remove(s)
                except ValueError:
                    pass
        else:
            ignore_done = []
        while len(not_done) > 0:
            s = not_done.pop(0)
            for letter in self.input_alphabet:
                try:
                    return_value = function(s.label(), letter)
                except LookupError:
                    continue
                if not hasattr(return_value, "pop"):
                    return_value = [return_value]
                try:
                    for (st_label, word) in return_value:
                        if not self.has_state(st_label):
                            not_done.append(self.add_state(st_label))
                        elif len(ignore_done) > 0:
                            u = self.state(st_label)
                            if u in ignore_done:
                                not_done.append(u)
                                ignore_done.remove(u)
                        self.add_transition(s, st_label,
                                            word_in=letter, word_out=word)
                except TypeError:
                    raise ValueError("The callback function for add_from_transition is expected to return a pair (new_state, output_label) or a list of such pairs. For the state %s and the input letter %s, it however returned %s, which is not acceptable." % (s.label(), letter, return_value))


    def add_transitions_from_function(self, function, labels_as_input=True):
        """
        Adds a transition if ``function(state, state)`` says that there is one.

        INPUT:

        - ``function`` -- a transition function. Given two states
          ``from_state`` and ``to_state`` (or their labels, if
          ``label_as_input`` is true), this function shall return a
          tuple ``(word_in, word_out)`` to add a transition from
          ``from_state`` to ``to_state`` with input and output labels
          ``word_in`` and ``word_out``, respectively. If no such
          addition is to be added, the transition function shall
          return ``None``.

        - ``label_as_input`` -- (default: ``True``)

        OUTPUT:

        Nothing.

        EXAMPLES::

            sage: F = FiniteStateMachine()
            sage: F.add_states(['A', 'B', 'C'])
            sage: def f(state1, state2):
            ....:     if state1 == 'C':
            ....:         return None
            ....:     return (0, 1)
            sage: F.add_transitions_from_function(f)
            sage: len(F.transitions())
            6
        """
        for s_from in self.iter_states():
            for s_to in self.iter_states():
                if labels_as_input:
                    t = function(s_from.label(), s_to.label())
                else:
                    t = function(s_from, s_to)
                if hasattr(t, '__getitem__'):
                    label_in = t[0]
                    try:
                        label_out = t[1]
                    except LookupError:
                        label_out = None
                    self.add_transition(s_from, s_to, label_in, label_out)


    def delete_transition(self, t):
        """
        Deletes a transition by removing it from the list of transitions of
        the state, where the transition starts.

        INPUT:

        - ``t`` -- a transition.

        OUTPUT:

        Nothing.

        EXAMPLES::

            sage: F = FiniteStateMachine([('A', 'B', 0), ('B', 'A', 1)])
            sage: F.delete_transition(('A', 'B', 0))
            sage: F.transitions()
            [Transition from 'B' to 'A': 1|-]
        """
        transition = self.transition(t)
        transition.from_state.transitions.remove(transition)


    def delete_state(self, s):
        """
        Deletes a state and all transitions coming or going to this state.

        INPUT:

        - ``s`` -- s has to be a label of a state or :class:`FSMState`.

        OUTPUT:

        Nothing.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: t1 = FSMTransition('A', 'B', 0)
            sage: t2 = FSMTransition('B', 'B', 1)
            sage: F = FiniteStateMachine([t1, t2])
            sage: F.delete_state('A')
            sage: F. transitions()
            [Transition from 'B' to 'B': 1|-]
        """
        state = self.state(s)
        for transition in self.transitions():
            if transition.to_state == state:
                self.delete_transition(transition)
        self._states_.remove(state)


    def remove_epsilon_transitions(self):
        """
        TESTS::

            sage: FiniteStateMachine().remove_epsilon_transitions()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


    def accessible_components(self):
        """
        Returns a new finite state machine with the accessible states
        of self and all transitions between those states.

        INPUT:

        Nothing.

        OUTPUT:

        A finite state machine with the accessible states of self and
        all transitions between those states.

        A state is accessible if there is a directed path from an
        initial state to the state. If self has no initial states then
        a copy of the finite state machine self is returned.

        EXAMPLES::

            sage: F = Automaton([(0, 0, 0), (0, 1, 1), (1, 1, 0), (1, 0, 1)],
            ....:               initial_states=[0])
            sage: F.accessible_components()
            Automaton with 2 states

        ::

            sage: F = Automaton([(0, 0, 1), (0, 0, 1), (1, 1, 0), (1, 0, 1)],
            ....:               initial_states=[0])
            sage: F.accessible_components()
            Automaton with 1 states
        """
        if len(self.initial_states()) == 0:
            return deepcopy(self)

        memo = {}
        def accessible(sf, read):
            trans = filter(lambda x: x.word_in[0] == read,
                           self.transitions(sf))
            return map(lambda x: (deepcopy(x.to_state, memo), x.word_out),
                       trans)

        new_initial_states=map(lambda x: deepcopy(x, memo),
                               self.initial_states())
        result = self.empty_copy()
        result.add_from_transition_function(accessible,
                                            initial_states=new_initial_states)
        for final_state in self.iter_final_states():
            try:
                new_final_state=result.state(final_state.label)
                new_final_state.is_final=True
            except LookupError:
                pass
        return result


    # *************************************************************************
    # creating new finite state machines
    # *************************************************************************


    def disjoint_union(self, other):
        """
        TESTS::

            sage: F = FiniteStateMachine([('A', 'A')])
            sage: FiniteStateMachine().disjoint_union(F)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


    def concatenation(self, other):
        """
        TESTS::

            sage: F = FiniteStateMachine([('A', 'A')])
            sage: FiniteStateMachine().concatenation(F)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


    def Kleene_closure(self):
        """
        TESTS::

            sage: FiniteStateMachine().Kleene_closure()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


    def intersection(self, other):
        """
        TESTS::

            sage: F = FiniteStateMachine([('A', 'A')])
            sage: FiniteStateMachine().intersection(F)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


    def product_FiniteStateMachine(self, other, function,
                                   new_input_alphabet=None,
                                   only_accessible_components=True):
        """
        Returns a new finite state machine whose states are
        pairs of states of the original finite state machines.

        INPUT:

        - ``other`` -- a finite state machine.

        - ``function`` has to accept two transitions from `A` to `B`
          and `C` to `D` and returns a pair ``(word_in, word_out)``
          which is the label of the transition `(A, C)` to `(B,
          D)`. If there is no transition from `(A, C)` to `(B, D)`,
          then ``function`` should raise a ``LookupError``.

        - ``new_input_alphabet`` (optional)-- the new input alphabet
          as a list.

        - ``only_accessible_components`` -- If true (default), then
          the result is piped through ``accessible_components``. If no
          ``new_input_alphabet`` is given, it is determined by
          ``determine_alphabets``.

        OUTPUT:

        A finite state machine whose states are pairs of states of the
        original finite state machines.

        The labels of the transitions are defined by ``function``.

        EXAMPLES::

            sage: F = Automaton([('A', 'B', 1), ('A', 'A', 0), ('B', 'A', 2)],
            ....:               initial_states=['A'], final_states=['B'],
            ....:               determine_alphabets=True)
            sage: G = Automaton([(1, 1, 1)], initial_states=[1], final_states=[1])
            sage: def addition(transition1, transition2):
            ....:     return (transition1.word_in[0] + transition2.word_in[0],
            ....:             None)
            sage: H = F.product_FiniteStateMachine(G, addition, [0, 1, 2, 3], only_accessible_components=False)
            sage: H.transitions()
            [Transition from ('A', 1) to ('B', 1): 2|-,
             Transition from ('A', 1) to ('A', 1): 1|-,
             Transition from ('B', 1) to ('A', 1): 3|-]
            sage: H1 = F.product_FiniteStateMachine(G, addition, [0, 1, 2, 3], only_accessible_components=False)
            sage: H1.states()[0].label()[0] is F.states()[0]
            True
            sage: H1.states()[0].label()[1] is G.states()[0]
            True

        ::

            sage: F = Automaton([(0,1,1/4), (0,0,3/4), (1,1,3/4), (1,0,1/4)],
            ....:                initial_states=[0] )
            sage: G = Automaton([(0,0,1), (1,1,3/4), (1,0,1/4)],
            ....:                initial_states=[0] )
            sage: H = F.product_FiniteStateMachine(
            ....:         G, lambda t1,t2: (t1.word_in[0]*t2.word_in[0], None))
            sage: H.states()
            [(0, 0), (1, 0)]

        ::

            sage: F = Automaton([(0,1,1/4), (0,0,3/4), (1,1,3/4), (1,0,1/4)],
            ....:                initial_states=[0] )
            sage: G = Automaton([(0,0,1), (1,1,3/4), (1,0,1/4)],
            ....:                initial_states=[0] )
            sage: H = F.product_FiniteStateMachine(G,
            ....:                                  lambda t1,t2: (t1.word_in[0]*t2.word_in[0], None),
            ....:                                  only_accessible_components=False)
            sage: H.states()
            [(0, 0), (1, 0), (0, 1), (1, 1)]
        """
        result = self.empty_copy()
        if new_input_alphabet is not None:
            result.input_alphabet = new_input_alphabet
        else:
            result.input_alphabet = None

        for transition1 in self.transitions():
            for transition2 in other.transitions():
                try:
                    word = function(transition1, transition2)
                except LookupError:
                    continue
                result.add_transition((transition1.from_state,
                                       transition2.from_state),
                                      (transition1.to_state,
                                       transition2.to_state),
                                      word[0], word[1])
        for state in result.states():
            if all(map(lambda s: s.is_initial, state.label())):
                state.is_initial = True
            if all(map(lambda s: s.is_final, state.label())):
                state.is_final = True

        if only_accessible_components:
            if new_input_alphabet is None:
                result.determine_alphabets()
            return result.accessible_components()
        else:
            return result


    def cartesian_product(self, other, only_accessible_components=True):
        """
        Returns a new finite state machine, which is the cartesian
        product of self and other.

        INPUT:

        - ``other`` -- a finite state machine.

        - ``only_accessible_components``

        OUTPUT:

        A new finite state machine, which is the cartesian product of
        self and other.

        The set of states of the new automaton is the cartesian
        product of the set of states of both given automata. There is
        a transition `((A, B), (C, D), a)` in the new automaton if
        there are transitions `(A, C, a)` and `(B, C, a)` in the old
        automata.

        EXAMPLES::

            sage: aut1 = Automaton([('1', '2', 1), ('2', '2', 1), ('2', '2', 0)],
            ....:                initial_states=['1'], final_states=['2'],
            ....:                determine_alphabets=True)
            sage: aut2 = Automaton([('A', 'A', 1), ('A', 'B', 1),
            ....:                 ('B', 'B', 0), ('B', 'A', 0)],
            ....:                initial_states=['A'], final_states=['B'],
            ....:                determine_alphabets=True)
            sage: res = aut1.cartesian_product(aut2)
            sage: res.transitions()
            [Transition from ('1', 'A') to ('2', 'A'): 1|-,
             Transition from ('1', 'A') to ('2', 'B'): 1|-,
             Transition from ('2', 'A') to ('2', 'A'): 1|-,
             Transition from ('2', 'A') to ('2', 'B'): 1|-,
             Transition from ('2', 'B') to ('2', 'B'): 0|-,
             Transition from ('2', 'B') to ('2', 'A'): 0|-]
        """
        def function(transition1, transition2):
            if transition1.word_in == transition2.word_in \
                    and transition1.word_out == transition2.word_out:
                return (transition1.word_in, transition1.word_out)
            else:
                raise LookupError

        return self.product_FiniteStateMachine(
            other, function,
            only_accessible_components = only_accessible_components)


    def composition(self, other, algorithm=None,
                    only_accessible_components=True):
        """
        Returns a new transducer which is the composition of self and
        other.

        INPUT:

        - ``other`` -- a transducer

        - ``algorithm`` -- can be one of the following

          - ``direct`` -- The composition is calculated directly.

            There can be arbitrarily many initial and final states,
            but the input and output labels must have length 1.


            WARNING: The output of other is fed into self.

          - ``explorative`` -- An explorative algorithm is used.

            At least the following restrictions apply, but are not
            checked:
            - both self and other have exactly one initial state
            - all input labels of transitions have length exactly 1

            The input alphabet of self has to be specified.

            This is a very limited implementation of composition.
            WARNING: The output of ``other`` is fed into ``self``.

          If algorithm is ``None``, then the algorithm is chosen
          automatically (at the moment always ``direct``).

        OUTPUT:

        A new transducer.

        The labels of the new finite state machine are pairs
        of states of the original finite state machines.

        EXAMPLES::

            sage: F = Transducer([('A', 'B', 1, 0), ('B', 'A', 0, 1)],
            ....:                initial_states=['A', 'B'], final_states=['B'],
            ....:                determine_alphabets=True)
            sage: G = Transducer([(1, 1, 1, 0), (1, 2, 0, 1),
            ....:                 (2, 2, 1, 1), (2, 2, 0, 0)],
            ....:                initial_states=[1], final_states=[2],
            ....:                determine_alphabets=True)
            sage: Hd = F.composition(G, algorithm='direct')
            sage: Hd.initial_states()
            [(1, 'B'), (1, 'A')]
            sage: Hd.transitions()
            [Transition from (1, 'B') to (1, 'A'): 1|1,
             Transition from (1, 'A') to (2, 'B'): 0|0,
             Transition from (2, 'B') to (2, 'A'): 0|1,
             Transition from (2, 'A') to (2, 'B'): 1|0]

        ::

            sage: F = Transducer([('A', 'B', 1, [1, 0]), ('B', 'B', 1, 1),
            ....:                 ('B', 'B', 0, 0)],
            ....:                initial_states=['A'], final_states=['B'])
            sage: G = Transducer([(1, 1, 0, 0), (1, 2, 1, 0),
            ....:                 (2, 2, 0, 1), (2, 1, 1, 1)],
            ....:                initial_states=[1], final_states=[1])
            sage: He = G.composition(F, algorithm='explorative')
            sage: He.transitions()
            [Transition from ('A', 1) to ('B', 2): 1|0,1,
             Transition from ('B', 2) to ('B', 2): 0|1,
             Transition from ('B', 2) to ('B', 1): 1|1,
             Transition from ('B', 1) to ('B', 1): 0|0,
             Transition from ('B', 1) to ('B', 2): 1|0]

        Be aware that after composition, different transitions may
        share the same output label (same python object)::

            sage: F = Transducer([ ('A','B',0,0), ('B','A',0,0)],
            ....:                initial_states=['A'],
            ....:                final_states=['A'])
            sage: F.transitions()[0].word_out is F.transitions()[1].word_out
            False
            sage: G = Transducer([('C','C',0,1)],)
            ....:                initial_states=['C'],
            ....:                final_states=['C'])
            sage: H = G.composition(F)
            sage: H.transitions()[0].word_out is H.transitions()[1].word_out
            True

        TESTS:

        Due to the limitations of the two algorithms the following
        (examples from above, but different algorithm used) does not
        give a full answer or does not work

        In the following, ``algorithm='explorative'`` is inadequate,
        as ``F`` has more than one initial state::

            sage: F = Transducer([('A', 'B', 1, 0), ('B', 'A', 0, 1)],
            ....:                initial_states=['A', 'B'], final_states=['B'],
            ....:                determine_alphabets=True)
            sage: G = Transducer([(1, 1, 1, 0), (1, 2, 0, 1),
            ....:                 (2, 2, 1, 1), (2, 2, 0, 0)],
            ....:                initial_states=[1], final_states=[2],
            ....:                determine_alphabets=True)
            sage: He = F.composition(G, algorithm='explorative')
            sage: He.initial_states()
            [(1, 'A')]
            sage: He.transitions()
            [Transition from (1, 'A') to (2, 'B'): 0|0,
             Transition from (2, 'B') to (2, 'A'): 0|1,
             Transition from (2, 'A') to (2, 'B'): 1|0]

        In the following example, ``algorithm='direct'`` is inappropriate
        as there are edges with output labels of length greater than 1::

            sage: F = Transducer([('A', 'B', 1, [1, 0]), ('B', 'B', 1, 1),
            ....:                 ('B', 'B', 0, 0)],
            ....:                initial_states=['A'], final_states=['B'])
            sage: G = Transducer([(1, 1, 0, 0), (1, 2, 1, 0),
            ....:                 (2, 2, 0, 1), (2, 1, 1, 1)],
            ....:                initial_states=[1], final_states=[1])
            sage: Hd = G.composition(F, algorithm='direct')
        """
        if algorithm == None:
            algorithm = 'direct'
        if algorithm == 'direct':
            return self._composition_direct_(other, only_accessible_components)
        elif algorithm == 'explorative':
            return self._composition_explorative_(other)
        else:
            raise ValueError, "Unknown algorithm %s." % (algorithm,)


    def _composition_direct_(self, other, only_accessible_components=True):
        """
        See :meth:`.composition` for details.

        TESTS::

            sage: F = Transducer([('A', 'B', 1, 0), ('B', 'A', 0, 1)],
            ....:                initial_states=['A', 'B'], final_states=['B'],
            ....:                determine_alphabets=True)
            sage: G = Transducer([(1, 1, 1, 0), (1, 2, 0, 1),
            ....:                 (2, 2, 1, 1), (2, 2, 0, 0)],
            ....:                initial_states=[1], final_states=[2],
            ....:                determine_alphabets=True)
            sage: Hd = F._composition_direct_(G)
            sage: Hd.initial_states()
            [(1, 'B'), (1, 'A')]
            sage: Hd.transitions()
            [Transition from (1, 'B') to (1, 'A'): 1|1,
             Transition from (1, 'A') to (2, 'B'): 0|0,
             Transition from (2, 'B') to (2, 'A'): 0|1,
             Transition from (2, 'A') to (2, 'B'): 1|0]

        """
        def function(transition1, transition2):
            if transition1.word_out == transition2.word_in:
                return (transition1.word_in, transition2.word_out)
            else:
                raise LookupError

        return other.product_FiniteStateMachine(
            self, function,
            only_accessible_components=only_accessible_components)


    def _composition_explorative_(self, other):
        """
        See :meth:`.composition` for details.

        TESTS::

            sage: F = Transducer([('A', 'B', 1, [1, 0]), ('B', 'B', 1, 1),
            ....:                 ('B', 'B', 0, 0)],
            ....:                initial_states=['A'], final_states=['B'])
            sage: G = Transducer([(1, 1, 0, 0), (1, 2, 1, 0),
            ....:                 (2, 2, 0, 1), (2, 1, 1, 1)],
            ....:                initial_states=[1], final_states=[1])
            sage: He = G._composition_explorative_(F)
            sage: He.transitions()
            [Transition from ('A', 1) to ('B', 2): 1|0,1,
             Transition from ('B', 2) to ('B', 2): 0|1,
             Transition from ('B', 2) to ('B', 1): 1|1,
             Transition from ('B', 1) to ('B', 1): 0|0,
             Transition from ('B', 1) to ('B', 2): 1|0]

        TODO:

        The explorative algorithm should be re-implemented using the
        process iterators of both finite state machines.
        """
        def composition_transition(state, input):
            (state1, state2) = state
            transition1 = None
            for transition in other.iter_transitions(state1):
                if transition.word_in == [input]:
                    transition1 = transition
                    break
            if transition1 is None:
                raise LookupError
            new_state1 = transition1.to_state
            new_state2 = state2
            output = []
            for o in transition1.word_out:
                transition2 = None
                for transition in self.iter_transitions(new_state2):
                    if transition.word_in == [o]:
                        transition2 = transition
                        break
                if transition2 is None:
                    raise LookupError
                new_state2 = transition2.to_state
                output += transition2.word_out
            return ((new_state1, new_state2), output)

        F = other.empty_copy()
        new_initial_states = [(other.initial_states()[0], self.initial_states()[0])]
        F.add_from_transition_function(composition_transition,
                                       initial_states=new_initial_states)

        for state in F.states():
            if all(map(lambda s: s.is_final, state.label())):
                state.is_final = True

        return F


    def input_projection(self):
        """
        Returns an automaton where the output of each transition of
        self is deleted.

        INPUT:

        Nothing

        OUTPUT:

        An automaton.

        EXAMPLES::

            sage: F = FiniteStateMachine([('A', 'B', 0, 1), ('A', 'A', 1, 1),
            ....:                         ('B', 'B', 1, 0)])
            sage: G = F.input_projection()
            sage: G.transitions()
            [Transition from 'A' to 'B': 0|-,
             Transition from 'A' to 'A': 1|-,
             Transition from 'B' to 'B': 1|-]
        """
        return self.projection(what='input')


    def output_projection(self):
        """
        Returns a automaton where the input of each transition of self
        is deleted and the new input is the original output.

        INPUT:

        Nothing

        OUTPUT:

        An automaton.

        EXAMPLES::

            sage: F = FiniteStateMachine([('A', 'B', 0, 1), ('A', 'A', 1, 1),
            ....:                         ('B', 'B', 1, 0)])
            sage: G = F.output_projection()
            sage: G.transitions()
            [Transition from 'A' to 'B': 1|-,
             Transition from 'A' to 'A': 1|-,
             Transition from 'B' to 'B': 0|-]
        """
        return self.projection(what='output')


    def projection(self, what='input'):
        """
        Returns an Automaton which transition labels are the projection
        of the transition labels of the input.

        INPUT:

        - ``what`` -- (default: ``input``) either ``input`` or ``output``.

        OUTPUT:

        An automaton.

        EXAMPLES::

            sage: F = FiniteStateMachine([('A', 'B', 0, 1), ('A', 'A', 1, 1),
            ....:                         ('B', 'B', 1, 0)])
            sage: G = F.projection(what='output')
            sage: G.transitions()
            [Transition from 'A' to 'B': 1|-,
             Transition from 'A' to 'A': 1|-,
             Transition from 'B' to 'B': 0|-]
        """
        new = Automaton()

        if what == 'input':
            new.input_alphabet = copy(self.input_alphabet)
        elif what == 'output':
            new.input_alphabet = copy(self.output_alphabet)
        else:
            raise NotImplementedError

        state_mapping = {}
        for state in self.iter_states():
            state_mapping[state] = new.add_state(deepcopy(state))
        for transition in self.iter_transitions():
            if what == 'input':
                new_word_in = transition.word_in
            elif what == 'output':
                new_word_in = transition.word_out
            else:
                raise NotImplementedError
            new.add_transition((state_mapping[transition.from_state],
                                state_mapping[transition.to_state],
                                new_word_in, None))
        return new


    def transposition(self):
        """
        Returns a new finite state machine, where all transitions of the
        input finite state machine are reversed.

        INPUT:

        Nothing.

        OUTPUT:

        A new finite state machine.

        EXAMPLES::

            sage: aut = Automaton([('A', 'A', 0), ('A', 'A', 1), ('A', 'B', 0)],
            ....:                 initial_states=['A'], final_states=['B'])
            sage: aut.transposition().transitions('B')
            [Transition from 'B' to 'A': 0|-]

        ::

            sage: aut = Automaton([('1', '1', 1), ('1', '2', 0), ('2', '2', 0)],
            ....:                 initial_states=['1'], final_states=['1', '2'])
            sage: aut.transposition().initial_states()
            ['1', '2']
        """
        transposition = self.empty_copy()

        for state in self.states():
            transposition.add_state(deepcopy(state))

        for transition in self.transitions():
            transposition.add_transition(
                transition.to_state.label(), transition.from_state.label(),
                transition.word_in, transition.word_out)

        for initial in self.initial_states():
            state = transposition.state(initial.label())
            if not initial.is_final:
                state.is_final = True
                state.is_initial = False

        for final in self.final_states():
            state = transposition.state(final.label())
            if not final.is_initial:
                state.is_final = False
                state.is_initial = True

        return transposition


    def split_transitions(self):
        """
        Returns a new transducer, where all transitions in self with input
        labels consisting of more than one letter
        are replaced by a path of the corresponding length.

        INPUT:

        Nothing.

        OUTPUT:

        A new transducer.

        EXAMPLES::

            sage: A = Transducer([('A', 'B', [1, 2, 3], 0)],
            ....:                initial_states=['A'], final_states=['B'])
            sage: A.split_transitions().states()
            [('A', ()), ('B', ()),
             ('A', (1,)), ('A', (1, 2))]
        """
        new = self.empty_copy()
        for state in self.states():
            new.add_state(FSMState((state, ()), is_initial=state.is_initial,
                                   is_final=state.is_final))
        for transition in self.transitions():
            for j in range(len(transition.word_in)-1):
                new.add_transition((
                        (transition.from_state, tuple(transition.word_in[:j])),
                        (transition.from_state, tuple(transition.word_in[:j+1])),
                        transition.word_in[j],
                        []))
            new.add_transition((
                    (transition.from_state, tuple(transition.word_in[:-1])),
                    (transition.to_state, ()),
                    transition.word_in[-1:],
                    transition.word_out))
        return new


    # *************************************************************************
    # simplifications
    # *************************************************************************


    def prepone_output(self):
        """
        Apply the following to each state `s` (except initial and
        final states) of the finite state machine as often as
        possible:

        If the letter a is prefix of the output label of all
        transitions from `s`, then remove it from all these labels and
        append it to all output labels of all transitions leading to
        `s`.

        We assume that the states have no output labels.

        INPUT:

        Nothing.

        OUTPUT:

        Nothing.

        EXAMPLES::

            sage: A = Transducer([('A', 'B', 1, 1), ('B', 'B', 0, 0), ('B', 'C', 1, 0)],
            ....:                initial_states=['A'], final_states=['C'])
            sage: A.prepone_output()
            sage: A.transitions()
            [Transition from 'A' to 'B': 1|1,0,
             Transition from 'B' to 'B': 0|0,
             Transition from 'B' to 'C': 1|-]

        ::

            sage: B = Transducer([('A', 'B', 0, 1), ('B', 'C', 1, [1, 1]), ('B', 'C', 0, 1)],
            ....:                initial_states=['A'], final_states=['C'])
            sage: B.prepone_output()
            sage: B.transitions()
            [Transition from 'A' to 'B': 0|1,1,
             Transition from 'B' to 'C': 1|1,
             Transition from 'B' to 'C': 0|-]

        If initial states are not labeled as such, unexpected results may be obtained::

            sage: C = Transducer([(0,1,0,0)])
            sage: C.prepone_output()
            prepone_output: All transitions leaving state 0 have an
            output label with prefix 0.  However, there is no inbound
            transition and it is not an initial state. This routine
            (possibly called by simplification) therefore erased this
            prefix from all outbound transitions.
            sage: C.transitions()
            [Transition from 0 to 1: 0|-]

        """
        def find_common_output(state):
            if len(filter(lambda transition: len(transition.word_out) == 0, self.transitions(state))) > 0:
                return ()
            first_letters = set(map(lambda transition: transition.word_out[0], self.transitions(state)))
            if len(first_letters) == 1:
                return (first_letters.pop(),)
            return ()

        changed = 1
        iteration = 0
        while changed > 0:
            changed = 0
            iteration += 1
            for state in self.states():
                if state.is_initial or state.is_final:
                    continue
                assert len(state.word_out) == 0, \
                    "prepone_output assumes that all states have empty output word, but state %s has output word %s" % \
                    (state, state.word_out)
                common_output = find_common_output(state)
                if len(common_output) > 0:
                    changed += 1
                    for transition in self.transitions(state):
                        assert transition.word_out[0] == common_output[0]
                        transition.word_out = transition.word_out[1:]
                    found_inbound_transition = False
                    for transition in self.transitions():
                        if transition.to_state == state:
                            transition.word_out = transition.word_out + [common_output[0]]
                            found_inbound_transition = True
                    if not found_inbound_transition:
                        print "prepone_output: All transitions leaving state %s have an output label with prefix %s. "\
                            "However, there is no inbound transition and it is not an initial state. "\
                            "This routine (possibly called by simplification) therefore erased this prefix from all "\
                            "outbound transitions." % (state, common_output[0])



    def equivalence_classes(self):
        """
        Returns a list of equivalence classes of states.

        INPUT:

        Nothing.

        OUTPUT:

        A list of equivalence classes of states.

        Two states `a` and `b` are equivalent, if and only if for each
        input label word_in the following holds:

        For paths `p_a` from `a` to `a'` with input label ``word_in``
        and output label ``word_out_a`` and `p_b` from `b` to `b'`
        with input label ``word_in`` and output label ``word_out_b``,
        we have ``word_out_a=word_out_b``, `a'` and `b'` have the same
        output label and are both final or both non-final.

        The function :meth:`.equivalence_classes` returns a list of
        the equivalence classes to this equivalence relation.

        This is one step of Moore's minimization algorithm.

        .. SEEALSO::

            :meth:`.minimization`

        EXAMPLES::

            sage: fsm = FiniteStateMachine([("A", "B", 0, 1), ("A", "B", 1, 0),
            ....:                           ("B", "C", 0, 0), ("B", "C", 1, 1),
            ....:                           ("C", "D", 0, 1), ("C", "D", 1, 0),
            ....:                           ("D", "A", 0, 0), ("D", "A", 1, 1)])
            sage: fsm.equivalence_classes()
            [['A', 'C'], ['B', 'D']]
        """

        # Two states a and b are said to be 0-equivalent, if their output
        # labels agree and if they are both final or non-final.
        #
        # For some j >= 1, two states a and b are said to be j-equivalent, if
        # they are j-1 equivalent and if for each element letter letter_in of
        # the input alphabet and transitions t_a from a with input label
        # letter_in, output label word_out_a to a' and t_b from b with input
        # label letter_in, output label word_out_b to b', we have
        # word_out_a=word_out_b and a' and b' are j-1 equivalent.

        # If for some j the relations j-1 equivalent and j-equivalent
        # coincide, then they are equal to the equivalence relation described
        # in the docstring.

        # classes_current holds the equivalence classes of j-equivalence,
        # classes_previous holds the equivalence classes of j-1 equivalence.

        if not self.is_deterministic():
            raise NotImplementedError, "Minimization via Moore's Algorithm is only implemented for deterministic finite state machines"

        # initialize with 0-equivalence
        classes_previous = []
        key_0 = lambda state: (state.is_final, state.word_out)
        states_grouped = full_group_by(self.states(), key=key_0)
        classes_current = [equivalence_class for
                           (key,equivalence_class) in states_grouped]

        while len(classes_current) != len(classes_previous):
            class_of = {}
            classes_previous = classes_current
            classes_current = []

            for k in range(len(classes_previous)):
                for state in classes_previous[k]:
                    class_of[state] = k

            key_current = lambda state: sorted(
                [(transition.word_in,
                  transition.word_out,
                  class_of[transition.to_state])
                 for transition in state.transitions])

            for class_previous in classes_previous:
                states_grouped = full_group_by(class_previous, key=key_current)
                classes_current.extend([equivalence_class for
                                       (key,equivalence_class) in states_grouped])

        return classes_current


    def quotient(self, classes):
        """
        Constructs the quotient with respect to the equivalence
        classes.

        INPUT:

        - ``classes`` is a list of equivalence classes of states.

        OUTPUT:

        A finite state machine.

        Assume that `c` is a class and `s`, `s'` are states in `c`. If
        there is a transition from `s` to some `t` with input label
        ``word_in`` and output label ``word_out``, then there has to
        be a transition from `s'` to some `t'` with input label
        ``word_in`` and output label ``word_out`` such that `s'` and
        `t'` are states of the same class `c'`. Then there is a
        transition from `c` to `c'` in the quotient with input label
        ``word_in`` and output label ``word_out``.

        Non-initial states may be merged with initial states, the
        resulting state is an initial state.

        All states in a class must have the same ``is_final`` and
        ``word_out`` values.

        EXAMPLES::

            sage: fsm = FiniteStateMachine([("A", "B", 0, 1), ("A", "B", 1, 0),
            ....:                           ("B", "C", 0, 0), ("B", "C", 1, 1),
            ....:                           ("C", "D", 0, 1), ("C", "D", 1, 0),
            ....:                           ("D", "A", 0, 0), ("D", "A", 1, 1)])
            sage: fsmq = fsm.quotient([[fsm.state("A"), fsm.state("C")],
            ....:                      [fsm.state("B"), fsm.state("D")]])
            sage: fsmq.transitions()
            [Transition from ('A', 'C')
                          to ('B', 'D'): 0|1,
             Transition from ('A', 'C')
                          to ('B', 'D'): 1|0,
             Transition from ('B', 'D')
                          to ('A', 'C'): 0|0,
             Transition from ('B', 'D')
                          to ('A', 'C'): 1|1]
            sage: fsmq.relabeled().transitions()
            [Transition from 0 to 1: 0|1,
             Transition from 0 to 1: 1|0,
             Transition from 1 to 0: 0|0,
             Transition from 1 to 0: 1|1]
            sage: fsmq1 = fsm.quotient(fsm.equivalence_classes())
            sage: fsmq1 == fsmq
            True
            sage: fsm.quotient([[fsm.state("A"), fsm.state("B"), fsm.state("C"), fsm.state("D")]])
            Traceback (most recent call last):
                ...
            ValueError: There is a transition Transition from 'B' to 'C': 0|0 in the original transducer, but no corresponding transition in the new transducer.
        """
        new = self.empty_copy()
        state_mapping = {}

        # Create new states and build state_mapping
        for c in classes:
            new_state = new.add_state(tuple(c))
            for state in c:
                state_mapping[state] = new_state

        # Copy data from old transducer
        for c in classes:
            new_state = state_mapping[c[0]]
            # copy all data from first class member
            new_state.is_initial = c[0].is_initial
            new_state.is_final = c[0].is_final
            new_state.word_out = c[0].word_out
            for transition in self.iter_transitions(c[0]):
                new.add_transition(
                    from_state=new_state,
                    to_state = state_mapping[transition.to_state],
                    word_in  = transition.word_in,
                    word_out = transition.word_out)

            # check that all class members have the same information (modulo classes)
            for state in c:
                new_state.is_initial = new_state.is_initial or state.is_initial
                assert new_state.is_final == state.is_final, "Class %s mixes final and non-final states" % (c,)
                assert new_state.word_out == state.word_out, "Class %s mixes different word_out" % (c,)
                assert len(self.transitions(state)) == len(new.transitions(new_state)), \
                    "Class %s has %d outgoing transitions, but class member %s has %d outgoing transitions" %  \
                    (c, len(new.transitions(new_state)), state, len(self.transitions(state)))
                for transition in self.transitions(state):
                    try:
                        new.transition((new_state, state_mapping[transition.to_state], transition.word_in, transition.word_out))
                    except LookupError:
                        raise ValueError, "There is a transition %s in the original transducer, but no corresponding transition in the new transducer." % transition
        return new


    # *************************************************************************
    # other
    # *************************************************************************


    def graph(self, edge_labels='words_in_out'):
        """
        Returns the graph of the finite state machine with labeled
        vertices and labeled edges.

        INPUT:

        - ``edge_label``: (default: ``'words_in_out'``) can be
             - ``'words_in_out'`` (labels will be strings ``'i|o'``)
             - a function with which takes as input a transition
               and outputs (returns) the label

        OUTPUT:

        A graph.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A')
            sage: T = Transducer()
            sage: T.graph()
            Digraph on 0 vertices
            sage: T.add_state(A)
            'A'
            sage: T.graph()
            Digraph on 1 vertex
            sage: T.add_transition(('A', 'A', 0, 1))
            Transition from 'A' to 'A': 0|1
            sage: T.graph()
            Looped digraph on 1 vertex
        """
        if edge_labels == 'words_in_out':
            label_fct = lambda t:t._in_out_label_()
        elif hasattr(edge_labels, '__call__'):
            label_fct = edge_labels
        else:
            raise TypeError, 'Wrong argument for edge_labels.'

        graph_data = []
        isolated_vertices = []
        for state in self.iter_states():
            transitions = state.transitions
            if len(transitions) == 0:
                isolated_vertices.append(state.label())
            for t in transitions:
                graph_data.append((t.from_state.label(), t.to_state.label(),
                                   label_fct(t)))

        G = DiGraph(graph_data)
        G.add_vertices(isolated_vertices)
        return G


    digraph = graph


    def plot(self):
        """
        Plots a graph of the finite state machine with labeled
        vertices and labeled edges.

        INPUT:

        Nothing.

        OUTPUT:

        A plot of the graph of the finite state machine.

        TESTS::

            sage: FiniteStateMachine([('A', 'A', 0)]).plot()
        """
        return self.graph(edge_labels='words_in_out').plot()


    def predecessors(self, state, valid_input=None):
        """
        Lists all predecessors of a state.

        INPUT:

        - ``state`` -- the state from which the predecessors should be
          listed.

        - ``valid_input`` -- If ``valid_input`` is a list, then we
          only consider transitions whose input labels are contained
          in ``valid_input``. ``state`` has to be a :class:`FSMState`
          (not a label of a state). If input labels of length larger
          than `1` are used, then ``valid_input`` has to be a list of
          lists.

        OUTPUT:

        A list of states.

        EXAMPLES::

            sage: A = Transducer([('I', 'A', 'a', 'b'), ('I', 'B', 'b', 'c'),
            ....:                 ('I', 'C', 'c', 'a'), ('A', 'F', 'b', 'a'),
            ....:                 ('B', 'F', ['c', 'b'], 'b'), ('C', 'F', 'a', 'c')],
            ....:                initial_states=['I'], final_states=['F'])
            sage: A.predecessors(A.state('A'))
            ['A', 'I']
            sage: A.predecessors(A.state('F'), valid_input=['b', 'a'])
            ['F', 'C', 'A', 'I']
            sage: A.predecessors(A.state('F'), valid_input=[['c', 'b'], 'a'])
            ['F', 'C', 'B']
        """
        if valid_input != None:
            valid_list = list()
            for input in valid_input:
                input_list = input
                if not isinstance(input_list, list):
                    input_list = [input]
                valid_list.append(input_list)
            valid_input = valid_list

        unhandeled_direct_predecessors = {s:[] for s in self.states() }
        for t in self.transitions():
            if valid_input is None or t.word_in in valid_input:
                unhandeled_direct_predecessors[t.to_state].append(t.from_state)
        done = []
        open = [state]
        while len(open) > 0:
            s = open.pop()
            candidates = unhandeled_direct_predecessors[s]
            if candidates is not None:
                open.extend(candidates)
                unhandeled_direct_predecessors[s] = None
                done.append(s)
        return(done)


#*****************************************************************************


def is_Automaton(FSM):
    """
    Tests whether or not ``FSM`` inherits from :class:`Automaton`.

    TESTS::

        sage: from sage.combinat.finite_state_machine import is_FiniteStateMachine, is_Automaton
        sage: is_Automaton(FiniteStateMachine())
        False
        sage: is_Automaton(Automaton())
        True
        sage: is_FiniteStateMachine(Automaton())
        True
    """
    return isinstance(FSM, Automaton)


class Automaton(FiniteStateMachine):
    """
    This creates an automaton, which is a finite state machine, whose
    transitions have input labels.

    An automaton has additional features like creating a deterministic
    and a minimized automaton.

    See class :class:`FiniteStateMachine` for more information.

    EXAMPLES:

    We can create an automaton recognizing even numbers (given in
    binary and read from left to right) in the following way::

        sage: A = Automaton([('P', 'Q', 0), ('P', 'P', 1),
        ....:                ('Q', 'P', 1), ('Q', 'Q', 0)],
        ....:               initial_states=['P'], final_states=['Q'])
        sage: A
        Automaton with 2 states
        sage: A([0])[0]
        True
        sage: A([1,1,0])[0]
        True
        sage: A([1,0,1])[0]
        False

    Note that the full output of the commands above gives more
    information and looks like this::

        sage: A([1,0,1])
        (False, 'P', [])

    TESTS::

        sage: Automaton()
        Automaton with 0 states
    """

    def _repr_(self):
        """
        Represents the finite state machine as "Automaton with n
        states" where n is the number of states.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: Automaton()._repr_()
            'Automaton with 0 states'
        """
        return "Automaton with %s states" % len(self._states_)

    def _latex_transition_label_(self, transition, format_function=latex):
        r"""
        Returns the proper transition label.

        INPUT:

        - ``transition`` - a transition

        - ``format_function'' - a function formatting the labels

        OUTPUT:

        A string.

        EXAMPLES::

            sage: F = Automaton([('A', 'B', 1)])
            sage: F._latex_()
            '\\begin{tikzpicture}[auto]\n\\node[state] (v0) at (3.000000,0.000000) {\\text{\\texttt{A}}}\n;\\node[state] (v1) at (-3.000000,0.000000) {\\text{\\texttt{B}}}\n;\\path[->] (v0) edge node {$\\left[1\\right]$} (v1);\n\\end{tikzpicture}'

        TESTS::

            sage: F = Automaton([('A', 'B', 0, 1)])
            sage: t = F.transitions()[0]
            sage: F._latex_transition_label_(t)
            \left[0\right]
        """
        return format_function(transition.word_in)

    def determinisation(self):
        """
        Returns a deterministic automaton which accepts the same input
        words as the original one.

        INPUT:

        Nothing.

        OUTPUT:

        A new automaton, which is deterministic.

        The labels of the states of the new automaton are frozensets of
        states of ``self``.

        The input alphabet must be specified. It is restricted to nice
        cases: input words have to have length at most `1`.

        EXAMPLES::

            sage: aut = Automaton([('A', 'A', 0), ('A', 'B', 1), ('B', 'B', 1)],
            ....:                 initial_states=['A'], final_states=['B'])
            sage: aut.determinisation().transitions()
            [Transition from frozenset(['A'])
                          to frozenset(['A']): 0|-,
             Transition from frozenset(['A'])
                          to frozenset(['B']): 1|-,
             Transition from frozenset(['B'])
                          to frozenset([]): 0|-,
             Transition from frozenset(['B'])
                          to frozenset(['B']): 1|-,
             Transition from frozenset([])
                          to frozenset([]): 0|-,
             Transition from frozenset([])
                          to frozenset([]): 1|-]

        ::

            sage: A = Automaton([('A', 'A', 1), ('A', 'A', 0), ('A', 'B', 1),
            ....:                ('B', 'C', 0), ('C', 'C', 1), ('C', 'C', 0)],
            ....:               initial_states=['A'], final_states=['C'])
            sage: A.determinisation().states()
            [frozenset(['A']), frozenset(['A', 'B']),
            frozenset(['A', 'C']), frozenset(['A', 'C', 'B'])]

        TESTS:

        This is from #15078, comment 13.

        ::

            sage: D = {'A': [('A', 'a'), ('B', 'a'), ('A', 'b')],
            ....:      'C': [], 'B': [('C', 'b')]}
            sage: auto = Automaton(D, initial_states=['A'], final_states=['C'])
            sage: auto.is_deterministic()
            False
            sage: auto.process(list('aaab'))
            (False, 'A', [])
            sage: auto.states()
            ['A', 'C', 'B']
            sage: auto.determinisation()
            Automaton with 3 states
        """
        for transition in self.transitions():
            assert len(transition.word_in) <= 1, "%s has input label of length > 1, which we cannot handle" % (transition,)

        epsilon_successors = {}
        direct_epsilon_successors = {}
        for state in self.states():
            direct_epsilon_successors[state] = set(map(lambda t:t.to_state,
                                                       filter(lambda transition: len(transition.word_in) == 0,
                                                              self.transitions(state)
                                                              )
                                                       )
                                                   )
            epsilon_successors[state] = set([state])

        old_count_epsilon_successors = 0
        count_epsilon_successors = len(epsilon_successors)

        while old_count_epsilon_successors < count_epsilon_successors:
            old_count_epsilon_successors = count_epsilon_successors
            count_epsilon_successors = 0
            for state in self.states():
                for direct_successor in direct_epsilon_successors[state]:
                    epsilon_successors[state] = epsilon_successors[state].union(epsilon_successors[direct_successor])
                count_epsilon_successors += len(epsilon_successors[state])


        def set_transition(states, letter):
            result = set()
            for state in states:
                for transition in self.transitions(state):
                    if transition.word_in == [letter]:
                        result.add(transition.to_state)
            result = result.union(*map(lambda s:epsilon_successors[s], result))
            return (frozenset(result), [])

        result = self.empty_copy()
        new_initial_states = [frozenset([state for state in self.initial_states()])]
        result.add_from_transition_function(set_transition,
                                            initial_states=new_initial_states)

        for state in result.states():
            if any(map(lambda s: s.is_final, state.label())):
                state.is_final = True


        return result


    def minimization(self, algorithm=None):
        """
        Returns the minimization of the input automaton as a new automaton.

        INPUT:

        - ``algorithm`` -- Either Moore's algorithm is used (default
          or ``algorithm='Moore'``), or Brzozowski's algorithm when
          ``algorithm='Brzozowski'``.

        OUTPUT:

        A new automaton.

        The resulting automaton is deterministic and has a minimal
        number of states.

        EXAMPLES::

            sage: A = Automaton([('A', 'A', 1), ('A', 'A', 0), ('A', 'B', 1),
            ....:                ('B', 'C', 0), ('C', 'C', 1), ('C', 'C', 0)],
            ....:               initial_states=['A'], final_states=['C'])
            sage: B = A.minimization(algorithm='Brzozowski')
            sage: B.transitions(B.states()[1])
            [Transition from frozenset([frozenset(['A', 'C', 'B']),
            frozenset(['C', 'B']), frozenset(['A', 'C'])]) to
            frozenset([frozenset(['A', 'C', 'B']), frozenset(['C', 'B']),
            frozenset(['A', 'C']), frozenset(['C'])]): 0|-,
            Transition from frozenset([frozenset(['A', 'C', 'B']),
            frozenset(['C', 'B']), frozenset(['A', 'C'])]) to
            frozenset([frozenset(['A', 'C', 'B']), frozenset(['C', 'B']),
            frozenset(['A', 'C'])]): 1|-]
            sage: len(B.states())
            3
            sage: C = A.minimization(algorithm='Brzozowski')
            sage: C.transitions(C.states()[1])
            [Transition from frozenset([frozenset(['A', 'C', 'B']),
            frozenset(['C', 'B']), frozenset(['A', 'C'])]) to
            frozenset([frozenset(['A', 'C', 'B']), frozenset(['C', 'B']),
            frozenset(['A', 'C']), frozenset(['C'])]): 0|-,
            Transition from frozenset([frozenset(['A', 'C', 'B']),
            frozenset(['C', 'B']), frozenset(['A', 'C'])]) to
            frozenset([frozenset(['A', 'C', 'B']), frozenset(['C', 'B']),
            frozenset(['A', 'C'])]): 1|-]
            sage: len(C.states())
            3

        ::

            sage: aut = Automaton([('1', '2', 'a'), ('2', '3', 'b'),
            ....:                  ('3', '2', 'a'), ('2', '1', 'b'),
            ....:                  ('3', '4', 'a'), ('4', '3', 'b')],
            ....:                  initial_states=['1'], final_states=['1'])
            sage: min = aut.minimization(algorithm='Brzozowski')
            sage: [len(min.states()), len(aut.states())]
            [3, 4]
            sage: min = aut.minimization(algorithm='Moore')
            Traceback (most recent call last):
            ...
            NotImplementedError: Minimization via Moore's Algorithm is only
            implemented for deterministic finite state machines
        """
        if algorithm is None or algorithm == "Moore":
            return self._minimization_Moore_()
        elif algorithm == "Brzozowski":
            return self._minimization_Brzozowski_()
        else:
            raise NotImplementedError, "Algorithm '%s' is not implemented. Choose 'Moore' or 'Brzozowski'" % algorithm


    def _minimization_Brzozowski_(self):
        """
        Returns a minimized automaton by using Brzozowski's algorithm.

        See also :meth:`.minimization`.

        TESTS::

            sage: A = Automaton([('A', 'A', 1), ('A', 'A', 0), ('A', 'B', 1),
            ....:                ('B', 'C', 0), ('C', 'C', 1), ('C', 'C', 0)],
            ....:               initial_states=['A'], final_states=['C'])
            sage: B = A._minimization_Brzozowski_()
            sage: len(B.states())
            3
        """
        return self.transposition().determinisation().transposition().determinisation()


    def _minimization_Moore_(self):
        """
        Returns a minimized automaton by using Brzozowski's algorithm.

        See also :meth:`.minimization`.

        TESTS::

            sage: aut = Automaton([('1', '2', 'a'), ('2', '3', 'b'),
            ....:                  ('3', '2', 'a'), ('2', '1', 'b'),
            ....:                  ('3', '4', 'a'), ('4', '3', 'b')],
            ....:                  initial_states=['1'], final_states=['1'])
            sage: min = aut._minimization_Moore_()
            Traceback (most recent call last):
            ...
            NotImplementedError: Minimization via Moore's Algorithm is only
            implemented for deterministic finite state machines
        """
        return self.quotient(self.equivalence_classes())


#*****************************************************************************


def is_Transducer(FSM):
    """
    Tests whether or not ``FSM`` inherits from :class:`Transducer`.

    TESTS::

        sage: from sage.combinat.finite_state_machine import is_FiniteStateMachine, is_Transducer
        sage: is_Transducer(FiniteStateMachine())
        False
        sage: is_Transducer(Transducer())
        True
        sage: is_FiniteStateMachine(Transducer())
        True
    """
    return isinstance(FSM, Transducer)


class Transducer(FiniteStateMachine):
    """
    This creates a transducer, which is a finite state machine, whose
    transitions have input and output labels.

    An transducer has additional features like creating a simplified
    transducer.

    See class :class:`FiniteStateMachine` for more information.

    EXAMPLES:

    We can create a transducer performing the addition of 1 (for
    numbers given in binary and read from right to left) in the
    following way::

        sage: T = Transducer([('C', 'C', 1, 0), ('C', 'N', 0, 1),
        ....:                 ('N', 'N', 0, 0), ('N', 'N', 1, 1)],
        ....:                initial_states=['C'], final_states=['N'])
        sage: T
        Transducer with 2 states
        sage: T([0])
        (True, 'N', [1])
        sage: T([1,1,0])
        (True, 'N', [0, 0, 1])
        sage: ZZ(T(15.digits(base=2)+[0])[2], base=2)
        16

    Note that we have padded the binary input sequence by a `0` so
    that the transducer can reach its final state.

    TESTS::

        sage: Transducer()
        Transducer with 0 states
    """

    def _repr_(self):
        """
        Represents the transducer as "Transducer with n states" where
        n is the number of states.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: Transducer()._repr_()
            'Transducer with 0 states'
        """
        return "Transducer with %s states" % len(self._states_)

    def _latex_transition_label_(self, transition, format_function=latex):
        r"""
        Returns the proper transition label.

        INPUT:

        - ``transition`` - a transition

        - ``format_function'' - a function formatting the labels

        OUTPUT:

        A string.

            sage: F = Transducer([('A', 'B', 1, 2)])
            sage: F._latex_()
            '\\begin{tikzpicture}[auto]\n\\node[state] (v0) at (3.000000,0.000000) {\\text{\\texttt{A}}}\n;\\node[state] (v1) at (-3.000000,0.000000) {\\text{\\texttt{B}}}\n;\\path[->] (v0) edge node {$\\left[1\\right] \\mid \\left[2\\right]$} (v1);\n\\end{tikzpicture}'

        TESTS::

            sage: F = Transducer([('A', 'B', 0, 1)])
            sage: t = F.transitions()[0]
            sage: F._latex_transition_label_(t)
            \left[0\right] \mid \left[1\right]
        """
        return (format_function(transition.word_in) + "\\mid"
                + format_function(transition.word_out))


    def simplification(self):
        """
        Returns a simplified transducer.

        INPUT:

        Nothing.

        OUTPUT:

        A new transducer.

        This function simplifies a transducer by Moore's algorithm,
        first moving common output labels of transitions leaving a
        state to output labels of transitions entering the state
        (cf. :meth:`.prepone_output`).

        The resulting transducer implements the same function as the
        original transducer.

        EXAMPLES::

            sage: fsm = Transducer([("A", "B", 0, 1), ("A", "B", 1, 0),
            ....:                           ("B", "C", 0, 0), ("B", "C", 1, 1),
            ....:                           ("C", "D", 0, 1), ("C", "D", 1, 0),
            ....:                           ("D", "A", 0, 0), ("D", "A", 1, 1)])
            sage: fsms = fsm.simplification()
            sage: fsms
            Transducer with 2 states
            sage: fsms.transitions()
            [Transition from ('A', 'C')
                          to ('B', 'D'): 0|1,
             Transition from ('A', 'C')
                          to ('B', 'D'): 1|0,
             Transition from ('B', 'D')
                          to ('A', 'C'): 0|0,
             Transition from ('B', 'D')
                          to ('A', 'C'): 1|1]
            sage: fsms.relabeled().transitions()
            [Transition from 0 to 1: 0|1,
             Transition from 0 to 1: 1|0,
             Transition from 1 to 0: 0|0,
             Transition from 1 to 0: 1|1]
        """
        fsm = deepcopy(self)
        fsm.prepone_output()
        return fsm.quotient(fsm.equivalence_classes())


#*****************************************************************************


def is_FSMProcessIterator(PI):
    """
    Tests whether or not ``PI`` inherits from :class:`FSMProcessIterator`.

    TESTS::

        sage: from sage.combinat.finite_state_machine import is_FSMProcessIterator, FSMProcessIterator
        sage: is_FSMProcessIterator(FSMProcessIterator(FiniteStateMachine()))
        Traceback (most recent call last):
        ...
        ValueError: No state is initial.
    """
    return isinstance(PI, FSMProcessIterator)


class FSMProcessIterator:
    """
    This class is for processing an input string on a finite state
    machine.

    An instance of this class is generated when
    :meth:`FiniteStateMachine.process` or
    :meth:`FiniteStateMachine.iter_process` of the finite state
    machine is invoked. It behaves like an iterator which, in each
    step, takes one letter of the input and runs (one step on) the
    finite state machine with this input. More precisely, in each
    step, the process iterator takes an outgoing transition of the
    current state, whose input label equals the input letter of the
    tape. The output label of the transition, if present, is written
    on the output tape.

    INPUT:

    - ``fsm`` -- The finite state machine on which the input should be
      processed.

    - ``input_tape`` -- The input tape. It can be anything that is
      iterable.

    - ``initial_state`` -- The initial state in which the machine
      starts. If this is ``None``, the unique inital state of the
      finite state machine is takes. If there are several, an error is
      reported.

    The process (iteration) stops if there are no more input letters
    on the tape. In this case a StopIteration exception is thrown. As
    result the following attributes are available:

    - ``accept_input`` -- Is True if the reached state is a final state.

    - ``current_state`` -- The current/reached state in the process.

    - ``output_tape`` -- The written output.

    Current values of those attribtes (except ``accept_input``) are
    (also) available during the iteration.

    OUTPUT:

    An iterator.

    EXAMPLES:

    The following transducer reads binary words and outputs a word,
    where blocks of ones are replaced by just a single one. Further
    only words that end with a zero are accepted.

    ::

        sage: T = Transducer({'A': [('A', 0, 0), ('B', 1, None)],
        ....:                 'B': [('B', 1, None), ('A', 0, [1, 0])]},
        ....:     initial_states=['A'], final_states=['A'])
        sage: input = [1, 1, 0, 0, 1, 0, 1, 1, 1, 0]
        sage: T.process(input)
        (True, 'A', [1, 0, 0, 1, 0, 1, 0])

    The function :meth:`FiniteStateMachine.process` created a new
    ``FSMProcessIterator``. We can do that manually, too, and get full
    access to the iteration process::

        sage: from sage.combinat.finite_state_machine import FSMProcessIterator
        sage: it = FSMProcessIterator(T, input_tape=input)
        sage: for _ in it:
        ....:     print (it.current_state, it.output_tape)
        ('B', [])
        ('B', [])
        ('A', [1, 0])
        ('A', [1, 0, 0])
        ('B', [1, 0, 0])
        ('A', [1, 0, 0, 1, 0])
        ('B', [1, 0, 0, 1, 0])
        ('B', [1, 0, 0, 1, 0])
        ('B', [1, 0, 0, 1, 0])
        ('A', [1, 0, 0, 1, 0, 1, 0])
        sage: it.accept_input
        True
    """
    def __init__(self, fsm, input_tape=None, initial_state=None):
        """
        See :class:`FSMProcessIterator` for more information.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMProcessIterator
            sage: inverter = Transducer({'A': [('A', 0, 1), ('A', 1, 0)]},
            ....:     initial_states=['A'], final_states=['A'])
            sage: it = FSMProcessIterator(inverter, input_tape=[0, 1])
            sage: for _ in it:
            ....:     pass
            sage: it.output_tape
            [1, 0]
        """
        self.fsm = fsm
        if initial_state is None:
            fsm_initial_states = self.fsm.initial_states()
            try:
                self.current_state = fsm_initial_states[0]
            except IndexError:
                raise ValueError, "No state is initial."
            if len(fsm_initial_states) > 1:
                raise ValueError, "Several initial states."
        else:
            self.current_state = initial_state
        self.output_tape = []
        if input_tape is None:
            self._input_tape_iter_ = iter([])
        else:
            if hasattr(input_tape, '__iter__'):
                self._input_tape_iter_ = iter(input_tape)
            else:
                raise ValueError, "Given input tape is not iterable."

    def __iter__(self):
        """
        Returns ``self``.

        TESTS::

            sage: from sage.combinat.finite_state_machine import FSMProcessIterator
            sage: inverter = Transducer({'A': [('A', 0, 1), ('A', 1, 0)]},
            ....:     initial_states=['A'], final_states=['A'])
            sage: it = FSMProcessIterator(inverter, input_tape=[0, 1])
            sage: id(it) == id(iter(it))
            True
        """
        return self

    def next(self):
        """
        Makes one step in processing the input tape.

        INPUT:

        Nothing.

        OUTPUT:

        It returns the taken transition. A ``StopIteration`` exception is
        thrown when there is nothing more to read.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMProcessIterator
            sage: inverter = Transducer({'A': [('A', 0, 1), ('A', 1, 0)]},
            ....:     initial_states=['A'], final_states=['A'])
            sage: it = FSMProcessIterator(inverter, input_tape=[0, 1])
            sage: it.next()
            Transition from 'A' to 'A': 0|1
            sage: it.next()
            Transition from 'A' to 'A': 1|0
            sage: it.next()
            Traceback (most recent call last):
            ...
            StopIteration
        """
        if hasattr(self, 'accept_input'):
            raise StopIteration
        try:
            # process current state
            transition = None
            try:
                transition = self.current_state.hook(
                    self.current_state, self)
            except AttributeError:
                pass
            self.write_word(self.current_state.word_out)

            # get next
            if not isinstance(transition, FSMTransition):
                next_word = []
                found = False

                try:
                    while not found:
                        next_word.append(self.read_letter())
                        try:
                            transition = self.get_next_transition(
                                next_word)
                            found = True
                        except ValueError:
                            pass
                except StopIteration:
                    # this means input tape is finished
                    if len(next_word) > 0:
                        self.accept_input = False
                    raise StopIteration

            # process transition
            try:
                transition.hook(transition, self)
            except AttributeError:
                pass
            self.write_word(transition.word_out)

            # go to next state
            self.current_state = transition.to_state

        except StopIteration:
            # this means, either input tape is finished or
            # someone has thrown StopIteration manually (in one
            # of the hooks)
            if not self.current_state.is_final:
                self.accept_input = False
            if not hasattr(self, 'accept_input'):
                self.accept_input = True
            raise StopIteration

        # return
        return transition

    def read_letter(self):
        """
        Reads a letter from the input tape.

        INPUT:

        Nothing.

        OUTPUT:

        A letter.

        Exception ``StopIteration`` is thrown if tape has reached
        the end.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMProcessIterator
            sage: inverter = Transducer({'A': [('A', 0, 1), ('A', 1, 0)]},
            ....:     initial_states=['A'], final_states=['A'])
            sage: it = FSMProcessIterator(inverter, input_tape=[0, 1])
            sage: it.read_letter()
            0
        """
        return self._input_tape_iter_.next()

    def write_letter(self, letter):
        """
        Writes a letter on the output tape.

        INPUT:

        - ``letter`` -- the letter to be written.

        OUTPUT:

        Nothing.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMProcessIterator
            sage: inverter = Transducer({'A': [('A', 0, 1), ('A', 1, 0)]},
            ....:     initial_states=['A'], final_states=['A'])
            sage: it = FSMProcessIterator(inverter, input_tape=[0, 1])
            sage: it.write_letter(42)
            sage: it.output_tape
            [42]
        """
        self.output_tape.append(letter)

    def write_word(self, word):
        """
        Writes a word on the output tape.

        INPUT:

        - ``word`` -- the word to be written.

        OUTPUT:

        Nothing.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMProcessIterator
            sage: inverter = Transducer({'A': [('A', 0, 1), ('A', 1, 0)]},
            ....:     initial_states=['A'], final_states=['A'])
            sage: it = FSMProcessIterator(inverter, input_tape=[0, 1])
            sage: it.write_word([4, 2])
            sage: it.output_tape
            [4, 2]
        """
        for letter in word:
            self.write_letter(letter)

    def get_next_transition(self, word_in):
        """
        Returns the next transition according to ``word_in``. It is
        assumed that we are in state ``self.current_state``.

        INPUT:

        - ``word_in`` -- the input word.

        OUTPUT:

        The next transition according to ``word_in``. It is assumed
        that we are in state ``self.current_state``.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMProcessIterator
            sage: inverter = Transducer({'A': [('A', 0, 1), ('A', 1, 0)]},
            ....:     initial_states=['A'], final_states=['A'])
            sage: it = FSMProcessIterator(inverter, input_tape=[0, 1])
            sage: it.get_next_transition([0])
            Transition from 'A' to 'A': 0|1
        """
        for transition in self.current_state.transitions:
            if transition.word_in == word_in:
                return transition
        raise ValueError


#*****************************************************************************


def setup_latex_preamble():
    """
    This function adds the package ``tikz`` with support for automata
    to the preamble of Latex so that the finite state machines can be
    drawn nicely.

    INPUT:

    Nothing.

    OUTPUT:

    Nothing.

    TESTS::

        sage: from sage.combinat.finite_state_machine import setup_latex_preamble
        sage: setup_latex_preamble()
    """
    latex.add_package_to_preamble_if_available('tikz')
    latex.add_to_preamble('\\usetikzlibrary{automata}')


#*****************************************************************************
