# sage_setup: distribution = sagemath-graphs
# sage.doctest: needs sage.graphs sage.modules
r"""
Finite state machines, automata, transducers

This module adds support for finite state machines, automata and
transducers.

For creating automata and transducers you can use classes

- :class:`Automaton` and :class:`Transducer`
  (or the more general class :class:`FiniteStateMachine`)

or the generators

- :class:`automata <sage.combinat.finite_state_machine_generators.AutomatonGenerators>` and
  :class:`transducers <sage.combinat.finite_state_machine_generators.TransducerGenerators>`

which contain :doc:`preconstructed and commonly used automata and transducers
<finite_state_machine_generators>`. See also the
:ref:`examples <finite_state_machine_examples>` below.


Contents
========

:class:`FiniteStateMachine` and derived classes :class:`Transducer` and :class:`Automaton`
------------------------------------------------------------------------------------------


Accessing parts of a finite state machine
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |


    :meth:`~FiniteStateMachine.state` | Get a state by its label
    :meth:`~FiniteStateMachine.states` | List of states
    :meth:`~FiniteStateMachine.iter_states` | Iterator over the states
    :meth:`~FiniteStateMachine.initial_states` | List of initial states
    :meth:`~FiniteStateMachine.iter_initial_states` | Iterator over initial states
    :meth:`~FiniteStateMachine.final_states` | List of final states
    :meth:`~FiniteStateMachine.iter_final_states` | Iterator over final states
    :meth:`~FiniteStateMachine.transition` | Get a transition by its states and labels
    :meth:`~FiniteStateMachine.transitions` | List of transitions
    :meth:`~FiniteStateMachine.iter_transitions` | Iterator over the transitions
    :meth:`~FiniteStateMachine.predecessors` | List of predecessors of a state
    :meth:`~FiniteStateMachine.induced_sub_finite_state_machine` | Induced sub-machine
    :meth:`~FiniteStateMachine.accessible_components` | Accessible components
    :meth:`~FiniteStateMachine.coaccessible_components` | Coaccessible components
    :meth:`~FiniteStateMachine.final_components` | Final components (connected components which cannot be left again)


(Modified) Copies
^^^^^^^^^^^^^^^^^

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~FiniteStateMachine.empty_copy` | Return an empty deep copy
    :meth:`~FiniteStateMachine.deepcopy` | Return a deep copy
    :meth:`~FiniteStateMachine.relabeled` | Return a relabeled deep copy
    :meth:`Automaton.with_output` | Extend an automaton to a transducer


Manipulation
^^^^^^^^^^^^

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~FiniteStateMachine.add_state` | Add a state
    :meth:`~FiniteStateMachine.add_states` | Add states
    :meth:`~FiniteStateMachine.delete_state` | Delete a state
    :meth:`~FiniteStateMachine.add_transition` | Add a transition
    :meth:`~FiniteStateMachine.add_transitions_from_function` | Add transitions
    :attr:`~FiniteStateMachine.input_alphabet` | Input alphabet
    :attr:`~FiniteStateMachine.output_alphabet` | Output alphabet
    :attr:`~FiniteStateMachine.on_duplicate_transition` | Hook for handling duplicate transitions
    :meth:`~FiniteStateMachine.add_from_transition_function` | Add transitions by a transition function
    :meth:`~FiniteStateMachine.delete_transition` | Delete a transition
    :meth:`~FiniteStateMachine.remove_epsilon_transitions` | Remove epsilon transitions (not implemented)
    :meth:`~FiniteStateMachine.split_transitions` | Split transitions with input words of length ``> 1``
    :meth:`~FiniteStateMachine.determine_alphabets` | Determine input and output alphabets
    :meth:`~FiniteStateMachine.determine_input_alphabet` | Determine input alphabet
    :meth:`~FiniteStateMachine.determine_output_alphabet` | Determine output alphabet
    :meth:`~FiniteStateMachine.construct_final_word_out` | Construct final output by implicitly reading trailing letters; cf. :meth:`~FiniteStateMachine.with_final_word_out`


Properties
^^^^^^^^^^

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~FiniteStateMachine.has_state` | Check for a state
    :meth:`~FiniteStateMachine.has_initial_state` | Check for an initial state
    :meth:`~FiniteStateMachine.has_initial_states` | Check for initial states
    :meth:`~FiniteStateMachine.has_final_state` | Check for a final state
    :meth:`~FiniteStateMachine.has_final_states` | Check for final states
    :meth:`~FiniteStateMachine.has_transition` | Check for a transition
    :meth:`~FiniteStateMachine.is_deterministic` | Check for a deterministic machine
    :meth:`~FiniteStateMachine.is_complete` | Check for a complete machine
    :meth:`~FiniteStateMachine.is_connected` | Check for a connected machine
    :meth:`Automaton.is_equivalent` | Check for equivalent automata
    :meth:`~FiniteStateMachine.is_Markov_chain` | Check for a Markov chain
    :meth:`~FiniteStateMachine.is_monochromatic` | Check whether the colors of all states are equal
    :meth:`~FiniteStateMachine.number_of_words` | Determine the number of successful paths
    :meth:`~FiniteStateMachine.asymptotic_moments` | Main terms of expectation and variance of sums of labels
    :meth:`~FiniteStateMachine.moments_waiting_time` | Moments of the waiting time for first true output
    :meth:`~FiniteStateMachine.epsilon_successors` | Epsilon successors of a state
    :meth:`Automaton.shannon_parry_markov_chain` | Compute Markov chain with Parry measure

Operations
^^^^^^^^^^

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~FiniteStateMachine.disjoint_union` | Disjoint union
    :meth:`~FiniteStateMachine.concatenation` | Concatenation
    :meth:`~FiniteStateMachine.kleene_star` | Kleene star
    :meth:`Automaton.complement` | Complement of an automaton
    :meth:`Automaton.intersection` | Intersection of automata
    :meth:`Transducer.intersection` | Intersection of transducers
    :meth:`Transducer.cartesian_product` | Cartesian product of a transducer with another finite state machine
    :meth:`~FiniteStateMachine.product_FiniteStateMachine` | Product of finite state machines
    :meth:`~FiniteStateMachine.composition` | Composition (output of other is input of self)
    :meth:`~FiniteStateMachine.input_projection` | Input projection (output is deleted)
    :meth:`~FiniteStateMachine.output_projection` | Output projection (old output is new input)
    :meth:`~FiniteStateMachine.projection` | Input or output projection
    :meth:`~FiniteStateMachine.transposition` | Transposition (all transitions are reversed)
    :meth:`~FiniteStateMachine.with_final_word_out` | Machine with final output constructed by implicitly reading trailing letters, cf. :meth:`~FiniteStateMachine.construct_final_word_out` for inplace version
    :meth:`Automaton.determinisation` | Determinisation of an automaton
    :meth:`~FiniteStateMachine.completion` | Completion of a finite state machine
    :meth:`~FiniteStateMachine.process` | Process input
    :meth:`~FiniteStateMachine.__call__` | Process input with shortened output
    :meth:`Automaton.process` | Process input of an automaton (output differs from general case)
    :meth:`Transducer.process` | Process input of a transducer (output differs from general case)
    :meth:`~FiniteStateMachine.iter_process` | Return process iterator
    :meth:`~FiniteStateMachine.language` | Return all possible output words
    :meth:`Automaton.language` | Return all possible accepted words


Simplification
^^^^^^^^^^^^^^

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |


    :meth:`~FiniteStateMachine.prepone_output` | Prepone output where possible
    :meth:`~FiniteStateMachine.equivalence_classes` | List of equivalent states
    :meth:`~FiniteStateMachine.quotient` | Quotient with respect to equivalence classes
    :meth:`~FiniteStateMachine.merged_transitions` | Merge transitions while adding input
    :meth:`~FiniteStateMachine.markov_chain_simplification` | Simplification of a Markov chain
    :meth:`Automaton.minimization` | Minimization of an automaton
    :meth:`Transducer.simplification` | Simplification of a transducer


Conversion
^^^^^^^^^^

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~FiniteStateMachine.adjacency_matrix` | (Weighted) adjacency :func:`matrix`
    :meth:`~FiniteStateMachine.graph` | Underlying :class:`DiGraph`
    :meth:`~FiniteStateMachine.plot` | Plot


LaTeX output
++++++++++++

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~FiniteStateMachine.latex_options` | Set options
    :meth:`~FiniteStateMachine.set_coordinates` | Set coordinates of the states
    :meth:`~FiniteStateMachine.default_format_transition_label` | Default formatting of words in transition labels
    :meth:`~FiniteStateMachine.format_letter_negative` | Format negative numbers as overlined number
    :meth:`~FiniteStateMachine.format_transition_label_reversed` | Format words in transition labels in reversed order

.. SEEALSO::

    :ref:`finite_state_machine_LaTeX_output`


:class:`FSMState`
-----------------

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :attr:`~FSMState.final_word_out` | Final output of a state
    :attr:`~FSMState.is_final` | Describe whether a state is final or not
    :attr:`~FSMState.is_initial` | Describe whether a state is initial or not
    :attr:`~FSMState.initial_probability` | Probability of starting in this state as part of a Markov chain
    :meth:`~FSMState.label` | Label of a state
    :meth:`~FSMState.relabeled` | Return a relabeled deep copy of a state
    :meth:`~FSMState.fully_equal` | Check whether two states are fully equal (including all attributes)


:class:`FSMTransition`
----------------------

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :attr:`~FSMTransition.from_state` | State in which transition starts
    :attr:`~FSMTransition.to_state` | State in which transition ends
    :attr:`~FSMTransition.word_in` | Input word of the transition
    :attr:`~FSMTransition.word_out` | Output word of the transition
    :meth:`~FSMTransition.deepcopy` | Return a deep copy of the transition


:class:`FSMProcessIterator`
---------------------------

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~FSMProcessIterator.next` | Make one step in processing the input tape
    :meth:`~FSMProcessIterator.preview_word` | Read a word from the input tape
    :meth:`~FSMProcessIterator.result` | Return the finished branches during process


Helper Functions
----------------

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :func:`equal` | Check whether all elements of ``iterator`` are equal
    :func:`full_group_by` | Group iterable by values of some key
    :func:`startswith` | Determine whether list starts with the given prefix
    :func:`FSMLetterSymbol` | Return a string associated to the input letter
    :func:`FSMWordSymbol` | Return a string associated to a word
    :func:`is_FSMState` | Test whether an object inherits from :class:`FSMState`
    :func:`is_FSMTransition` | Test whether an object inherits from :class:`FSMTransition`
    :func:`is_FiniteStateMachine` | Test whether an object inherits from :class:`FiniteStateMachine`
    :func:`duplicate_transition_ignore` |  Default function for handling duplicate transitions
    :func:`duplicate_transition_raise_error` | Raise error when inserting a duplicate transition
    :func:`duplicate_transition_add_input` | Add input when inserting a duplicate transition


.. _finite_state_machine_examples:

Examples
========

We start with a general :class:`FiniteStateMachine`. Later there will
be also an :class:`Automaton` and a :class:`Transducer`.

A simple finite state machine
-----------------------------

We can easily create a finite state machine by

::

    sage: fsm = FiniteStateMachine()
    sage: fsm
    Empty finite state machine

By default this is the empty finite state machine, so not very
interesting. Let's create and add some states and transitions::

    sage: day = fsm.add_state('day')
    sage: night = fsm.add_state('night')
    sage: sunrise = fsm.add_transition(night, day)
    sage: sunset = fsm.add_transition(day, night)

Let us look at ``sunset`` more closely::

    sage: sunset
    Transition from 'day' to 'night': -|-

Note that could also have created and added the transitions directly
by::

    sage: fsm.add_transition('day', 'night')
    Transition from 'day' to 'night': -|-

This would have had added the states automatically, since they are
present in the transitions.

Anyhow, we got the following finite state machine::

    sage: fsm
    Finite state machine with 2 states

We can also obtain the underlying :class:`directed graph <DiGraph>` by

::

    sage: fsm.graph()
    Looped multi-digraph on 2 vertices

To visualize a finite state machine, we can use
:func:`~sage.misc.latex.latex` and run the result through LaTeX,
see the section on :ref:`finite_state_machine_LaTeX_output`
below.

Alternatively, we could have created the finite state machine above
simply by

::

    sage: FiniteStateMachine([('night', 'day'), ('day', 'night')])
    Finite state machine with 2 states

See :class:`FiniteStateMachine` for a lot of possibilities to create
finite state machines.

.. _finite_state_machine_recognizing_NAFs_example:

A simple Automaton (recognizing NAFs)
---------------------------------------

We want to build an automaton which recognizes non-adjacent forms
(NAFs), i.e., sequences which have no adjacent nonzeros.
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

    sage: NAF([0])
    True
    sage: NAF([0, 1])
    True
    sage: NAF([1, -1])
    False
    sage: NAF([0, -1, 0, 1])
    True
    sage: NAF([0, -1, -1, -1, 0])
    False
    sage: NAF([-1, 0, 0, 1, 1])
    False

Alternatively, we could call that by

::

    sage: NAF.process([0, -1, 0, 1])
    (True, 'B')

which gives additionally the state in which we arrived.

We can also let an automaton act on a :doc:`word <words/words>`::

    sage: # needs sage.combinat
    sage: W = Words([-1, 0, 1]); W
    Finite and infinite words over {-1, 0, 1}
    sage: w = W([1, 0, 1, 0, -1]); w
    word: 1,0,1,0,-1
    sage: NAF(w)
    True

Recognizing NAFs via Automata Operations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Alternatively, we can use automata operations to recognize NAFs; for
simplicity, we only use the input alphabet ``[0, 1]``. On the one
hand, we can construct such an automaton by forbidding the word
``11``::

    sage: forbidden = automata.ContainsWord([1, 1], input_alphabet=[0, 1])
    sage: NAF_negative = forbidden.complement()
    sage: NAF_negative([1, 1, 0, 1])
    False
    sage: NAF_negative([1, 0, 1, 0, 1])
    True

On the other hand, we can write this as a regular expression and
translate that into automata operations::

    sage: zero = automata.Word([0])
    sage: one = automata.Word([1])
    sage: epsilon = automata.EmptyWord(input_alphabet=[0, 1])
    sage: NAF_positive = (zero + one*zero).kleene_star() * (epsilon + one)

We check that the two approaches are equivalent::

    sage: NAF_negative.is_equivalent(NAF_positive)
    True

.. SEEALSO::

    :meth:`~sage.combinat.finite_state_machine_generators.AutomatonGenerators.ContainsWord`,
    :meth:`~sage.combinat.finite_state_machine_generators.AutomatonGenerators.Word`,
    :meth:`~Automaton.complement`,
    :meth:`~FiniteStateMachine.kleene_star`,
    :meth:`~sage.combinat.finite_state_machine_generators.AutomatonGenerators.EmptyWord`,
    :meth:`~Automaton.is_equivalent`.

.. _finite_state_machine_LaTeX_output:

LaTeX output
------------

We can visualize a finite state machine by converting it to LaTeX by
using the usual function :func:`~sage.misc.latex.latex`. Within LaTeX,
TikZ is used for typesetting the graphics, see the
:wikipedia:`PGF/TikZ`.

::

    sage: print(latex(NAF))  # abs tol 1e-3
    \begin{tikzpicture}[auto, initial text=, >=latex]
    \node[state, accepting, initial] (v0) at (3.000000, 0.000000) {$\text{\texttt{A}}$};
    \node[state, accepting] (v1) at (-3.000000, 0.000000) {$\text{\texttt{B}}$};
    \path[->] (v0) edge[loop above] node {$0$} ();
    \path[->] (v0.185.00) edge node[rotate=360.00, anchor=north] {$1, -1$} (v1.355.00);
    \path[->] (v1.5.00) edge node[rotate=0.00, anchor=south] {$0$} (v0.175.00);
    \end{tikzpicture}

We can turn this into a graphical representation.

::

    sage: view(NAF) # not tested

To actually see this, use the live documentation in the Sage notebook
and execute the cells in this and the previous section.

Several options can be set to customize the output, see
:meth:`~FiniteStateMachine.latex_options` for details. In particular,
we use :meth:`~FiniteStateMachine.format_letter_negative` to format
`-1` as `\overline{1}`.

::

    sage: NAF.latex_options(
    ....:     coordinates={'A': (0, 0),
    ....:                  'B': (6, 0)},
    ....:     initial_where={'A': 'below'},
    ....:     format_letter=NAF.format_letter_negative,
    ....:     format_state_label=lambda x:
    ....:         r'\mathcal{%s}' % x.label()
    ....: )
    sage: print(latex(NAF))
    \begin{tikzpicture}[auto, initial text=, >=latex]
    \node[state, accepting, initial, initial where=below] (v0) at (0.000000, 0.000000) {$\mathcal{A}$};
    \node[state, accepting] (v1) at (6.000000, 0.000000) {$\mathcal{B}$};
    \path[->] (v0) edge[loop above] node {$0$} ();
    \path[->] (v0.5.00) edge node[rotate=0.00, anchor=south] {$1, \overline{1}$} (v1.175.00);
    \path[->] (v1.185.00) edge node[rotate=360.00, anchor=north] {$0$} (v0.355.00);
    \end{tikzpicture}
    sage: view(NAF) # not tested

To use the output of :func:`~sage.misc.latex.latex` in your own
`\LaTeX` file, you have to include

.. code-block:: latex

    \usepackage{tikz}
    \usetikzlibrary{automata}

into the preamble of your file.

A simple transducer (binary inverter)
-------------------------------------

Let's build a simple transducer, which rewrites a binary word by
inverting each bit::

    sage: inverter = Transducer({'A': [('A', 0, 1), ('A', 1, 0)]},
    ....:     initial_states=['A'], final_states=['A'])

We can look at the states and transitions::

    sage: inverter.states()
    ['A']
    sage: for t in inverter.transitions():
    ....:     print(t)
    Transition from 'A' to 'A': 0|1
    Transition from 'A' to 'A': 1|0

Now we apply a word to it and see what the transducer does::

    sage: inverter([0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1])
    [1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0]

``True`` means, that we landed in a final state, that state is labeled
``'A'``, and we also got an output.

.. _finite_state_machine_division_by_3_example:


Transducers and (in)finite Words
--------------------------------

A transducer can also act on everything iterable, in particular, on
Sage's :doc:`words <words/words>`.

::

    sage: W = Words([0, 1]); W
    Finite and infinite words over {0, 1}

Let us take the inverter from the previous section and feed some
finite word into it::

    sage: w = W([1, 1, 0, 1]); w
    word: 1101
    sage: inverter(w)
    word: 0010

We see that the output is again a word (this is a consequence of
calling :meth:`~Transducer.process` with ``automatic_output_type``).

We can even input something infinite like an infinite word::

    sage: tm = words.ThueMorseWord(); tm
    word: 0110100110010110100101100110100110010110...
    sage: inverter(tm)
    word: 1001011001101001011010011001011001101001...


A transducer which performs division by `3` in binary
-----------------------------------------------------

Now we build a transducer, which divides a binary number by `3`.
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

which assumes reading a binary number from left to right.
We get the transducer with

::

    sage: D = Transducer(f, initial_states=[0], final_states=[0],
    ....:                input_alphabet=[0, 1])

Let us try to divide `12` by `3`::

    sage: D([1, 1, 0, 0])
    [0, 1, 0, 0]

Now we want to divide `13` by `3`::

    sage: D([1, 1, 0, 1])
    Traceback (most recent call last):
    ...
    ValueError: Invalid input sequence.

The raised :exc:`ValueError` means `13` is not divisible by `3`.

.. _finite_state_machine_gray_code_example:

Gray Code
---------

The Gray code is a binary :wikipedia:`numeral system <Numeral_system>`
where two successive values differ in only one bit, cf. the
:wikipedia:`Gray_code`. The Gray code of an integer `n` is obtained by
a bitwise xor between the binary expansion of `n` and the binary
expansion of `\lfloor n/2\rfloor`; the latter corresponds to a
shift by one position in binary.

The purpose of this example is to construct a transducer converting the
standard binary expansion to the Gray code by translating this
construction into operations with transducers.

For this construction, the least significant digit is at
the left-most position.
Note that it is easier to shift everything to
the right first, i.e., multiply by `2` instead of building
`\lfloor n/2\rfloor`. Then, we take the input xor with the right
shift of the input and forget the first letter.

We first construct a transducer shifting the binary expansion to the
right. This requires storing the previously read digit in a state.

::

    sage: def shift_right_transition(state, digit):
    ....:     if state == 'I':
    ....:         return (digit, None)
    ....:     else:
    ....:         return (digit, state)
    sage: shift_right_transducer = Transducer(
    ....:     shift_right_transition,
    ....:     initial_states=['I'],
    ....:     input_alphabet=[0, 1],
    ....:     final_states=[0])
    sage: shift_right_transducer.transitions()
    [Transition from 'I' to 0: 0|-,
     Transition from 'I' to 1: 1|-,
     Transition from 0 to 0: 0|0,
     Transition from 0 to 1: 1|0,
     Transition from 1 to 0: 0|1,
     Transition from 1 to 1: 1|1]
    sage: shift_right_transducer([0, 1, 1, 0])
    [0, 1, 1]
    sage: shift_right_transducer([1, 0, 0])
    [1, 0]

The output of the shifts above look a bit weird (from a right-shift
transducer, we would expect, for example, that ``[1, 0, 0]`` was
mapped to ``[0, 1, 0]``), since we write ``None`` instead of the zero
at the left. Further, note that only `0` is listed as a final state
as we have to enforce that a most significant zero is read as the last
input letter in order to flush the last digit::

    sage: shift_right_transducer([0, 1, 0, 1])
    Traceback (most recent call last):
    ...
    ValueError: Invalid input sequence.

Next, we construct the transducer performing the xor operation. We also
have to take ``None`` into account as our ``shift_right_transducer``
waits one iteration until it starts writing output. This corresponds
with our intention to forget the first letter.

::

    sage: def xor_transition(state, digits):
    ....:    if digits[0] is None or digits[1] is None:
    ....:        return (0, None)
    ....:    else:
    ....:        return (0, digits[0].__xor__(digits[1]))
    sage: from itertools import product
    sage: xor_transducer = Transducer(
    ....:    xor_transition,
    ....:    initial_states=[0],
    ....:    final_states=[0],
    ....:    input_alphabet=list(product([None, 0, 1], [0, 1])))
    sage: xor_transducer.transitions()
    [Transition from 0 to 0: (None, 0)|-,
     Transition from 0 to 0: (None, 1)|-,
     Transition from 0 to 0: (0, 0)|0,
     Transition from 0 to 0: (0, 1)|1,
     Transition from 0 to 0: (1, 0)|1,
     Transition from 0 to 0: (1, 1)|0]
    sage: xor_transducer([(None, 0), (None, 1), (0, 0), (0, 1), (1, 0), (1, 1)])
    [0, 1, 1, 0]
    sage: xor_transducer([(0, None)])
    Traceback (most recent call last):
    ...
    ValueError: Invalid input sequence.

The transducer computing the Gray code is then constructed as a
:meth:`Cartesian product <Transducer.cartesian_product>` between the
shifted version and the original input (represented here by the
``shift_right_transducer`` and the :meth:`identity transducer
<sage.combinat.finite_state_machine_generators.TransducerGenerators.Identity>`,
respectively). This Cartesian product is then fed into the
``xor_transducer`` as a :meth:`composition
<FiniteStateMachine.composition>` of transducers.

::

    sage: product_transducer = shift_right_transducer.cartesian_product(transducers.Identity([0, 1]))
    sage: Gray_transducer = xor_transducer(product_transducer)

We use :meth:`~FiniteStateMachine.construct_final_word_out` to make sure that all output
is written; otherwise, we would have to make sure that a sufficient number of trailing
zeros is read.

::

    sage: Gray_transducer.construct_final_word_out([0])
    sage: Gray_transducer.transitions()
    [Transition from (('I', 0), 0) to ((0, 0), 0): 0|-,
     Transition from (('I', 0), 0) to ((1, 0), 0): 1|-,
     Transition from ((0, 0), 0) to ((0, 0), 0): 0|0,
     Transition from ((0, 0), 0) to ((1, 0), 0): 1|1,
     Transition from ((1, 0), 0) to ((0, 0), 0): 0|1,
     Transition from ((1, 0), 0) to ((1, 0), 0): 1|0]

There is a :meth:`prepackaged transducer
<sage.combinat.finite_state_machine_generators.TransducerGenerators.GrayCode>`
for Gray code, let's see whether they agree. We have to use
:meth:`~FiniteStateMachine.relabeled` to relabel our states with
integers.

::

    sage: constructed = Gray_transducer.relabeled()
    sage: packaged = transducers.GrayCode()
    sage: constructed == packaged
    True

Finally, we check that this indeed computes the Gray code of the first
10 nonnegative integers.

::

    sage: for n in srange(10):
    ....:     Gray_transducer(n.bits())
    []
    [1]
    [1, 1]
    [0, 1]
    [0, 1, 1]
    [1, 1, 1]
    [1, 0, 1]
    [0, 0, 1]
    [0, 0, 1, 1]
    [1, 0, 1, 1]


Using the hook-functions
------------------------

Let's use the :ref:`previous example "division by
3" <finite_state_machine_division_by_3_example>` to demonstrate the optional
state and transition parameters ``hook``.

First, we define what those functions should do. In our case, this is
just saying in which state we are and which transition we take

::

    sage: def state_hook(process, state, output):
    ....:     print("We are now in State %s." % (state.label(),))
    sage: from sage.combinat.finite_state_machine import FSMWordSymbol
    sage: def transition_hook(transition, process):
    ....:     print("Currently we go from %s to %s, "
    ....:           "reading %s and writing %s." % (
    ....:                transition.from_state, transition.to_state,
    ....:                FSMWordSymbol(transition.word_in),
    ....:                FSMWordSymbol(transition.word_out)))

Now, let's add these hook-functions to the existing transducer::

    sage: for s in D.iter_states():
    ....:     s.hook = state_hook
    sage: for t in D.iter_transitions():
    ....:     t.hook = transition_hook

Rerunning the process again now gives the following output::

    sage: D.process([1, 1, 0, 1], check_epsilon_transitions=False)
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
hook-functions. In the following, we will use those hooks more
seriously.

.. WARNING::

   The hooks of the states are also called while exploring the epsilon
   successors of a state (during processing). In the example above, we
   used ``check_epsilon_transitions=False`` to avoid this (and also
   therefore got a cleaner output).

.. WARNING::

   The arguments used when calling a hook have changed in
   :issue:`16538` from ``hook(state, process)`` to
   ``hook(process, state, output)``.

Detecting sequences with same number of `0` and `1`
---------------------------------------------------

Suppose we have a binary input and want to accept all sequences with
the same number of `0` and `1`. This cannot be done with a finite
automaton. Anyhow, we can make usage of the hook functions to extend
our finite automaton by a counter::

    sage: from sage.combinat.finite_state_machine import FSMState, FSMTransition
    sage: C = FiniteStateMachine()
    sage: def update_counter(process, state, output):
    ....:     try:
    ....:         l = process.preview_word()
    ....:     except RuntimeError:
    ....:         raise StopIteration
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

See also methods :meth:`Automaton.process` and
:meth:`Transducer.process` (or even
:meth:`FiniteStateMachine.process`), the explanation of the parameter
``hook`` and the examples in :class:`FSMState` and
:class:`FSMTransition`, and the description and examples in
:class:`FSMProcessIterator` for more information on processing and
hooks.

REFERENCES:

.. [HKP2015] Clemens Heuberger, Sara Kropf, and Helmut Prodinger,
   *Output sum of transducers: Limiting distribution and periodic
   fluctuation*,
   `Electron. J. Combin. 22 (2015), #P2.19 <http://www.combinatorics.org/ojs/index.php/eljc/article/view/v22i2p19>`_.

.. [HKW2015] Clemens Heuberger, Sara Kropf and Stephan Wagner,
   *Variances and Covariances in the Central Limit Theorem for the Output
   of a Transducer*, European J. Combin. 49 (2015), 167-187,
   :doi:`10.1016/j.ejc.2015.03.004`.

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
- Daniel Krenn (2013-11-11): comments from :issue:`15078` included:
    docstring of FiniteStateMachine rewritten, Automaton and Transducer
    inherited from FiniteStateMachine
- Daniel Krenn (2013-11-25): documentation improved according to
    comments from :issue:`15078`
- Clemens Heuberger, Daniel Krenn, Sara Kropf (2014-02-21--2014-07-18):
  A huge bunch of improvements. Details see
  :issue:`15841`, :issue:`15847`, :issue:`15848`, :issue:`15849`, :issue:`15850`, :issue:`15922`, :issue:`15923`, :issue:`15924`,
  :issue:`15925`, :issue:`15928`, :issue:`15960`, :issue:`15961`, :issue:`15962`, :issue:`15963`, :issue:`15975`, :issue:`16016`,
  :issue:`16024`, :issue:`16061`, :issue:`16128`, :issue:`16132`, :issue:`16138`, :issue:`16139`, :issue:`16140`, :issue:`16143`,
  :issue:`16144`, :issue:`16145`, :issue:`16146`, :issue:`16191`, :issue:`16200`, :issue:`16205`, :issue:`16206`, :issue:`16207`,
  :issue:`16229`, :issue:`16253`, :issue:`16254`, :issue:`16255`, :issue:`16266`, :issue:`16355`, :issue:`16357`, :issue:`16387`,
  :issue:`16425`, :issue:`16539`, :issue:`16555`, :issue:`16557`, :issue:`16588`, :issue:`16589`, :issue:`16666`, :issue:`16668`,
  :issue:`16674`, :issue:`16675`, :issue:`16677`.
- Daniel Krenn (2015-09-14): cleanup :issue:`18227`

ACKNOWLEDGEMENT:

- Clemens Heuberger, Daniel Krenn and Sara Kropf are supported by the
  Austrian Science Fund (FWF): P 24644-N26.

Methods
=======
"""
# ****************************************************************************
# Copyright (C) 2012--2015 Clemens Heuberger <clemens.heuberger@aau.at>
#               2012--2015 Daniel Krenn <dev@danielkrenn.at>
#               2012--2015 Sara Kropf <sara.kropf@aau.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from IPython.lib.pretty import pretty
import itertools
from collections import defaultdict, deque, namedtuple, OrderedDict
# Use isinstance(x, Iterable) to test whether x is iterable, and
# use isinstance(x, Mapping) to test whether x is a dict.
from collections.abc import Iterator, Iterable, Mapping
from copy import copy, deepcopy

from sage.graphs.digraph import DiGraph
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_function
from sage.misc.latex import latex
from sage.misc.verbose import verbose
from sage.misc.sageinspect import sage_getargspec
from sage.rings.integer_ring import ZZ
from sage.rings.real_mpfr import RR
from sage.structure.sage_object import SageObject


def full_group_by(l, key=None):
    """
    Group iterable ``l`` by values of ``key``.

    INPUT:

    - iterable ``l``
    - key function ``key``

    OUTPUT:

    A list of pairs ``(k, elements)`` such that ``key(e)=k`` for all
    ``e`` in ``elements``.

    This is similar to :func:`itertools.groupby` except that lists are
    returned instead of iterables and no prior sorting is required.

    We do not require

    - that the keys are sortable (in contrast to the
      approach via :func:`sorted` and :func:`itertools.groupby`) and
    - that the keys are hashable (in contrast to the
      implementation proposed in `<https://stackoverflow.com/a/15250161>`_).

    However, it is required

    - that distinct keys have distinct ``str``-representations.

    The implementation is inspired by
    `<https://stackoverflow.com/a/15250161>`_, but non-hashable keys are
    allowed.

    EXAMPLES::

        sage: # needs sage.symbolic
        sage: from sage.combinat.finite_state_machine import full_group_by
        sage: t = [2/x, 1/x, 2/x]
        sage: r = full_group_by([0, 1, 2], key=lambda i: t[i])
        sage: sorted(r, key=lambda p: p[1])
        [(2/x, [0, 2]), (1/x, [1])]
        sage: from itertools import groupby
        sage: for k, elements in groupby(sorted([0, 1, 2],
        ....:                            key=lambda i:t[i]),
        ....:                            key=lambda i:t[i]):
        ....:     print("{} {}".format(k, list(elements)))
        2/x [0]
        1/x [1]
        2/x [2]

    Note that the behavior is different from :func:`itertools.groupby`
    because neither `1/x<2/x` nor `2/x<1/x` does hold.

    Here, the result ``r`` has been sorted in order to guarantee a
    consistent order for the doctest suite.
    """
    if key is None:
        key = lambda x: x
    elements = defaultdict(list)
    original_keys = {}
    for item in l:
        k = key(item)
        s = str(k)
        if s in original_keys:
            if original_keys[s] != k:
                raise ValueError(
                    'two distinct elements with representation {}'.format(s))
        else:
            original_keys[s] = k
        elements[s].append(item)
    return [(original_keys[s], elements[s]) for s in sorted(elements)]


def equal(iterator):
    """
    Check whether all elements of ``iterator`` are equal.

    INPUT:

    - ``iterator`` -- an iterator of the elements to check

    OUTPUT: boolean

    This implements `<https://stackoverflow.com/a/3844832/1052778>`_.

    EXAMPLES::

        sage: from sage.combinat.finite_state_machine import equal
        sage: equal([0, 0, 0])
        True
        sage: equal([0, 1, 0])
        False
        sage: equal([])
        True
        sage: equal(iter([None, None]))
        True

    We can test other properties of the elements than the elements
    themselves. In the following example, we check whether all tuples
    have the same lengths::

        sage: equal(len(x) for x in [(1, 2), (2, 3), (3, 1)])
        True
        sage: equal(len(x) for x in [(1, 2), (1, 2, 3), (3, 1)])
        False
    """
    try:
        iterator = iter(iterator)
        first = next(iterator)
        return all(first == rest for rest in iterator)
    except StopIteration:
        return True


def startswith(list_, prefix):
    """
    Determine whether list starts with the given prefix.

    INPUT:

    - ``list_`` -- list
    - ``prefix`` -- list representing the prefix

    OUTPUT: boolean

    Similar to :meth:`str.startswith`.

    EXAMPLES::

        sage: from sage.combinat.finite_state_machine import startswith
        sage: startswith([1, 2, 3], [1, 2])
        True
        sage: startswith([1], [1, 2])
        False
        sage: startswith([1, 3, 2], [1, 2])
        False
    """
    if len(prefix) > len(list_):
        return False
    return list_[:len(prefix)] == prefix

# ****************************************************************************


FSMEmptyWordSymbol = '-'
EmptyWordLaTeX = r'\varepsilon'
EndOfWordLaTeX = r'\$'
tikz_automata_where = {"right": 0,
                       "above": 90,
                       "left": 180,
                       "below": 270}


def FSMLetterSymbol(letter):
    """
    Return a string associated to the input letter.

    INPUT:

    - ``letter`` -- the input letter or ``None`` (representing the
      empty word)

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
    Return a string of ``word``. It may returns the symbol of the
    empty word ``FSMEmptyWordSymbol``.

    INPUT:

    - ``word`` -- the input word

    OUTPUT: string of ``word``

    EXAMPLES::

        sage: from sage.combinat.finite_state_machine import FSMWordSymbol
        sage: FSMWordSymbol([0, 1, 1])
        '0,1,1'
    """
    if not isinstance(word, list):
        return FSMLetterSymbol(word)
    if not word:
        return FSMEmptyWordSymbol
    return ','.join(FSMLetterSymbol(letter) for letter in word)


# ****************************************************************************


def is_FSMState(S):
    """
    Test whether or not ``S`` inherits from :class:`FSMState`.

    TESTS::

        sage: from sage.combinat.finite_state_machine import is_FSMState, FSMState
        sage: is_FSMState(FSMState('A'))
        doctest:warning...
        DeprecationWarning: The function is_FSMState is deprecated; use 'isinstance(..., FSMState)' instead.
        See https://github.com/sagemath/sage/issues/38032 for details.
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(38032, "The function is_FSMState is deprecated; use 'isinstance(..., FSMState)' instead.")
    return isinstance(S, FSMState)


class FSMState(SageObject):
    """
    Class for a state of a finite state machine.

    INPUT:

    - ``label`` -- the label of the state

    - ``word_out`` -- (default: ``None``) a word that is written when
      the state is reached

    - ``is_initial`` -- (default: ``False``)

    - ``is_final`` -- (default: ``False``)

    - ``final_word_out`` -- (default: ``None``) a word that is written when
      the state is reached as the last state of some input; only for final
      states.

    - ``initial_probability`` -- (default: ``None``) the probability of
      starting in this state if it is a state of a Markov chain

    - ``hook`` -- (default: ``None``) a function which is called when
      the state is reached during processing input. It takes two input
      parameters: the first is the current state (to allow using the same
      hook for several states), the second is the current process
      iterator object (to have full access to everything; e.g. the
      next letter from the input tape can be read in). It can output
      the next transition, i.e. the transition to take next. If it
      returns ``None`` the process iterator chooses. Moreover, this
      function can raise a ``StopIteration`` exception to stop
      processing of a finite state machine the input immediately. See
      also the example below.

    - ``color`` -- (default: ``None``) in order to distinguish states,
      they can be given an arbitrary "color" (an arbitrary object).
      This is used in :meth:`FiniteStateMachine.equivalence_classes`:
      states of different colors are never considered to be
      equivalent. Note that :meth:`Automaton.determinisation` requires
      that ``color`` is hashable.

    - ``allow_label_None`` -- boolean (default: ``False``); if ``True`` allows
      also ``None`` as label. Note that a state with label ``None`` is used in
      :class:`FSMProcessIterator`.

    OUTPUT: a state of a finite state machine

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

    We can also define a final output word of a final state which is
    used if the input of a transducer leads to this state. Such final
    output words are used in subsequential transducers. ::

        sage: C = FSMState('state 3', is_final=True, final_word_out='end')
        sage: C.final_word_out
        ['end']

    The final output word can be a single letter, ``None`` or a list of
    letters::

        sage: A = FSMState('A')
        sage: A.is_final = True
        sage: A.final_word_out = 2
        sage: A.final_word_out
        [2]
        sage: A.final_word_out = [2, 3]
        sage: A.final_word_out
        [2, 3]

    Only final states can have a final output word which is not
    ``None``::

        sage: B = FSMState('B')
        sage: B.final_word_out is None
        True
        sage: B.final_word_out = 2
        Traceback (most recent call last):
        ...
        ValueError: Only final states can have a final output word,
        but state B is not final.

    Setting the ``final_word_out`` of a final state to ``None`` is the
    same as setting it to ``[]`` and is also the default for a final
    state::

        sage: C = FSMState('C', is_final=True)
        sage: C.final_word_out
        []
        sage: C.final_word_out = None
        sage: C.final_word_out
        []
        sage: C.final_word_out = []
        sage: C.final_word_out
        []

    It is not allowed to use ``None`` as a label::

        sage: from sage.combinat.finite_state_machine import FSMState
        sage: FSMState(None)
        Traceback (most recent call last):
        ...
        ValueError: Label None reserved for a special state,
        choose another label.

    This can be overridden by::

        sage: FSMState(None, allow_label_None=True)
        None

    Note that :meth:`Automaton.determinisation` requires that ``color``
    is hashable::

        sage: A = Automaton([[0, 0, 0]], initial_states=[0])
        sage: A.state(0).color = []
        sage: A.determinisation()
        Traceback (most recent call last):
        ...
        TypeError: unhashable type: 'list'
        sage: A.state(0).color = ()
        sage: A.determinisation()
        Automaton with 1 state

    We can use a hook function of a state to stop processing. This is
    done by raising a ``StopIteration`` exception. The following code
    demonstrates this::

        sage: T = Transducer([(0, 1, 9, 'a'), (1, 2, 9, 'b'),
        ....:                 (2, 3, 9, 'c'), (3, 4, 9, 'd')],
        ....:                initial_states=[0],
        ....:                final_states=[4],
        ....:                input_alphabet=[9])
        sage: def stop(process, state, output):
        ....:     raise StopIteration()
        sage: T.state(3).hook = stop
        sage: T.process([9, 9, 9, 9])
        (False, 3, ['a', 'b', 'c'])

    TESTS:

    Test for ``is_initial``::

        sage: T = Automaton([(0,0,0)])
        sage: T.initial_states()
        []
        sage: T.state(0).is_initial = True
        sage: T.initial_states()
        [0]

    Test for ``initial_probability``::

        sage: from sage.combinat.finite_state_machine import FSMState
        sage: S = FSMState('state', initial_probability=1/3)
        sage: S.initial_probability
        1/3
    """
    is_initial = False  # Describes whether the state is initial
    initial_probability = None
    # The probability of starting in this state if it is part of a Markov chain

    def __init__(self, label, word_out=None,
                 is_initial=False, is_final=False, final_word_out=None,
                 initial_probability=None,
                 hook=None, color=None, allow_label_None=False):
        """
        See :class:`FSMState` for more information.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: FSMState('final', is_final=True)
            'final'

        TESTS::

            sage: A = FSMState('A', is_final=True)
            sage: A.final_word_out
            []
            sage: A.is_final = True
            sage: A = FSMState('A', is_final=True, final_word_out='end')
            sage: A.final_word_out
            ['end']
            sage: A = FSMState('A', is_final=True,
            ....:              final_word_out=['e', 'n', 'd'])
            sage: A.final_word_out
            ['e', 'n', 'd']
            sage: A = FSMState('A', is_final=True, final_word_out=[])
            sage: A.final_word_out
            []
            sage: A = FSMState('A', is_final=True, final_word_out=None)
            sage: A.final_word_out
            []
            sage: A = FSMState('A', is_final=False)
            sage: A.final_word_out is None
            True
            sage: A.is_final = False
            sage: A = FSMState('A', is_final=False, final_word_out='end')
            Traceback (most recent call last):
            ...
            ValueError: Only final states can have a final output word,
            but state A is not final.
            sage: A = FSMState('A', is_final=False,
            ....:              final_word_out=['e', 'n', 'd'])
            Traceback (most recent call last):
            ...
            ValueError: Only final states can have a final output word,
            but state A is not final.
            sage: A = FSMState('A', is_final=False, final_word_out=None)
            sage: A.final_word_out is None
            True
            sage: A = FSMState('A', is_final=False, final_word_out=[])
            Traceback (most recent call last):
            ...
            ValueError: Only final states can have a final output word,
            but state A is not final.
        """
        if not allow_label_None and label is None:
            raise ValueError("Label None reserved for a special state, "
                             "choose another label.")
        self._label_ = label

        if isinstance(word_out, list):
            self.word_out = word_out
        elif word_out is not None:
            self.word_out = [word_out]
        else:
            self.word_out = []

        self.is_initial = is_initial
        self._final_word_out_ = None
        self.is_final = is_final
        self.final_word_out = final_word_out
        self.initial_probability = initial_probability

        if hook is not None:
            if callable(hook):
                self.hook = hook
            else:
                raise TypeError('Wrong argument for hook.')

        self.color = color

    def __lt__(self, other):
        """
        Return ``True`` if label of ``self`` is less than label of ``other``.

        INPUT:

        - ``other`` -- a state

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: FSMState(0) < FSMState(1)
            True
        """
        return self.label() < other.label()

    @property
    def final_word_out(self):
        """
        The final output word of a final state which is written if the
        state is reached as the last state of the input of the finite
        state machine. For a non-final state, the value is ``None``.

        ``final_word_out`` can be a single letter, a list or ``None``,
        but for a final-state, it is always saved as a list.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A', is_final=True, final_word_out=2)
            sage: A.final_word_out
            [2]
            sage: A.final_word_out = 3
            sage: A.final_word_out
            [3]
            sage: A.final_word_out = [3, 4]
            sage: A.final_word_out
            [3, 4]
            sage: A.final_word_out = None
            sage: A.final_word_out
            []
            sage: B = FSMState('B')
            sage: B.final_word_out is None
            True

        A non-final state cannot have a final output word::

            sage: B.final_word_out = [3, 4]
            Traceback (most recent call last):
            ...
            ValueError: Only final states can have a final
            output word, but state B is not final.
        """
        return self._final_word_out_

    @final_word_out.setter
    def final_word_out(self, final_word_out):
        """
        Set the value of the final output word of a final state.

        INPUT:

        - ``final_word_out`` -- list; any element or ``None``

        OUTPUT: nothing

        TESTS::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: B = FSMState('B')
            sage: B.final_word_out = []
            Traceback (most recent call last):
            ...
            ValueError: Only final states can have a final
            output word, but state B is not final.
            sage: B.final_word_out = None
            sage: B.final_word_out is None
            True

        The exception is raised also when the initial state is a tuple
        (see :issue:`18990`)::

            sage: A = Transducer(initial_states=[(0, 0)])
            sage: A.state((0, 0)).final_word_out = []
            Traceback (most recent call last):
            ...
            ValueError: Only final states can have a final output word,
            but state (0, 0) is not final.

        No exception is raised if we set the state to be a final one::

            sage: A.state((0, 0)).is_final=True
            sage: A.state((0, 0)).final_word_out = []
            sage: A.state((0, 0)).final_word_out == []
            True
        """
        if not self.is_final:
            if final_word_out is not None:
                raise ValueError("Only final states can have a "
                                 "final output word, but state %s is not final."
                                 % (self.label(),))
            else:
                self._final_word_out_ = None
        elif isinstance(final_word_out, list):
            self._final_word_out_ = final_word_out
        elif final_word_out is not None:
            self._final_word_out_ = [final_word_out]
        else:
            self._final_word_out_ = []

    @property
    def is_final(self):
        """
        Describe whether the state is final or not.

        ``True`` if the state is final and ``False`` otherwise.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A', is_final=True, final_word_out=3)
            sage: A.is_final
            True
            sage: A.is_final = False
            Traceback (most recent call last):
            ...
            ValueError: State A cannot be non-final, because it has a
            final output word. Only final states can have a final output
            word.
            sage: A.final_word_out = None
            sage: A.is_final = False
            sage: A.is_final
            False
        """
        return self.final_word_out is not None

    @is_final.setter
    def is_final(self, is_final):
        """
        Define the state as a final state or a non-final state.

        INPUT:

        - ``is_final`` -- ``True`` if the state should be final and
          ``False`` otherwise

        OUTPUT: nothing

        TESTS::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A', is_final=True)
            sage: A.final_word_out
            []
            sage: A.is_final = False
            sage: A.final_word_out is None
            True
            sage: A = FSMState('A', is_final=True, final_word_out='a')
            sage: A.is_final = False
            Traceback (most recent call last):
            ...
            ValueError: State A cannot be non-final, because it has a
            final output word. Only final states can have a final output
            word.

        The exception is raised also when the final state is a tuple
        (see :issue:`18990`)::

            sage: A = Transducer(final_states=[(0, 0)])
            sage: A.state((0, 0)).final_word_out = [1]
            sage: A.state((0, 0)).is_final = False
            Traceback (most recent call last):
            ...
            ValueError: State (0, 0) cannot be non-final, because it has
            a final output word. Only final states can have a final
            output word.

        No exception is raised if we empty the final_word_out of the
        state::

            sage: A.state((0, 0)).final_word_out = []
            sage: A.state((0, 0)).is_final = False
            sage: A.state((0, 0)).is_final
            False

            sage: A = FSMState('A', is_final=True, final_word_out=[])
            sage: A.is_final = False
            sage: A.final_word_out is None
            True
        """
        if is_final and self.final_word_out is None:
            self._final_word_out_ = []
        elif not is_final:
            if not self.final_word_out:
                self._final_word_out_ = None
            else:
                raise ValueError("State %s cannot be non-final, because it "
                                 "has a final output word. Only final states "
                                 "can have a final output word. "
                                 % (self.label(),))

    def label(self):
        """
        Return the label of the state.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('state')
            sage: A.label()
            'state'
        """
        return self._label_

    def __copy__(self):
        """
        Return a (shallow) copy of the state.

        OUTPUT: a new state

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A')
            sage: A.is_initial = True
            sage: A.is_final = True
            sage: A.final_word_out = [1]
            sage: A.color = 'green'
            sage: A.initial_probability = 1/2
            sage: B = copy(A)
            sage: B.fully_equal(A)
            True
            sage: A.label() is B.label()
            True
            sage: A.is_initial is B.is_initial
            True
            sage: A.is_final is B.is_final
            True
            sage: A.final_word_out is B.final_word_out
            True
            sage: A.color is B.color
            True
            sage: A.initial_probability is B.initial_probability
            True
        """
        new = FSMState(self.label(), self.word_out,
                       self.is_initial, self.is_final,
                       color=self.color,
                       final_word_out=self.final_word_out,
                       initial_probability=self.initial_probability)
        if hasattr(self, 'hook'):
            new.hook = self.hook
        return new

    copy = __copy__

    def __deepcopy__(self, memo):
        """
        Return a deep copy of the state.

        INPUT:

        - ``memo`` -- dictionary storing already processed elements

        OUTPUT: a new state

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
        new.color = deepcopy(self.color, memo)
        new.final_word_out = deepcopy(self.final_word_out, memo)
        new.initial_probability = deepcopy(self.initial_probability, memo)
        return new

    def deepcopy(self, memo=None):
        """
        Return a deep copy of the state.

        INPUT:

        - ``memo`` -- (default: ``None``) a dictionary storing already
          processed elements

        OUTPUT: a new state

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState((1, 3), color=[1, 2],
            ....:              is_final=True, final_word_out=3,
            ....:              initial_probability=1/3)
            sage: B = deepcopy(A)
            sage: B
            (1, 3)
            sage: B.label() == A.label()
            True
            sage: B.label is A.label
            False
            sage: B.color == A.color
            True
            sage: B.color is A.color
            False
            sage: B.is_final == A.is_final
            True
            sage: B.is_final is A.is_final
            True
            sage: B.final_word_out == A.final_word_out
            True
            sage: B.final_word_out is A.final_word_out
            False
            sage: B.initial_probability == A.initial_probability
            True
        """
        return deepcopy(self, memo)

    def relabeled(self, label, memo=None):
        """
        Return a deep copy of the state with a new label.

        INPUT:

        - ``label`` -- the label of new state

        - ``memo`` -- (default: ``None``) a dictionary storing already
          processed elements

        OUTPUT: a new state

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

    def __getstate__(self):
        """
        Return state for pickling excluding outgoing transitions.

        OUTPUT: a dictionary

        Outgoing transitions are in fact stored in states,
        but must be pickled by the finite state machine
        in order to avoid deep recursion.

        TESTS::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: S = FSMState('S')
            sage: S == loads(dumps(S))
            True
            sage: A = Automaton([[0, 1, 0]])
            sage: T = A.state(0)
            sage: T1 = loads(dumps(T))
            sage: T == T1
            True
            sage: T.transitions
            [Transition from 0 to 1: 0|-]
            sage: T1.transitions
            Traceback (most recent call last):
            ...
            AttributeError: 'FSMState' object has no attribute 'transitions'...
            sage: A1 = loads(dumps(A))
            sage: all(A.state(j) == A1.state(j) for j in [0, 1])
            True
            sage: all(A.state(j).transitions == A1.state(j).transitions
            ....:     for j in [0, 1])
            True
        """
        odict = self.__dict__.copy()  # copy the dict since we change it
        try:
            del odict['transitions']  # remove transitions entry
        except KeyError:
            pass                      # Standalone FSMState has no transitions
        return odict

    def __hash__(self):
        """
        Return a hash value for the object.

        OUTPUT: the hash of this state

        TESTS::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A')
            sage: hash(A) #random
            -269909568
        """
        return hash(self.label())

    def _repr_(self):
        """
        Return the string "label".

        OUTPUT: string

        TESTS::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: FSMState('A')
            'A'
        """
        return pretty(self.label())

    def __eq__(self, other):
        """
        Return ``True`` if two states are the same, i.e., if they have
        the same labels.

        INPUT:

        - ``self`` -- a state

        - ``other`` -- a state

        OUTPUT: boolean

        Note that the hooks and whether the states are initial or
        final are not checked. To fully compare two states (including
        these attributes), use :meth:`.fully_equal`.

        As only the labels are used when hashing a state, only the
        labels can actually be compared by the equality relation.
        Note that the labels are unique within one finite state machine,
        so this may only lead to ambiguities when comparing states
        belonging to different finite state machines.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A')
            sage: B = FSMState('A', is_initial=True)
            sage: A == B
            True
        """
        if not isinstance(other, FSMState):
            return False
        return self.label() == other.label()

    def __ne__(self, other):
        """
        Test for inequality, complement of __eq__.

        INPUT:

        - ``self`` -- a state

        - ``other`` -- a state

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A', is_initial=True)
            sage: B = FSMState('A', is_final=True)
            sage: A != B
            False
        """
        return not (self == other)

    def fully_equal(self, other, compare_color=True):
        """
        Check whether two states are fully equal, i.e., including all
        attributes except ``hook``.

        INPUT:

        - ``self`` -- a state

        - ``other`` -- a state

        - ``compare_color`` -- boolean (default: ``True``); if ``True`` colors
          are compared as well, otherwise not

        OUTPUT: boolean

        Note that usual comparison by ``==`` does only compare the labels.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A')
            sage: B = FSMState('A', is_initial=True)
            sage: A.fully_equal(B)
            False
            sage: A == B
            True
            sage: A.is_initial = True; A.color = 'green'
            sage: A.fully_equal(B)
            False
            sage: A.fully_equal(B, compare_color=False)
            True
        """
        color = not compare_color or self.color == other.color
        return (self == other and
                self.is_initial == other.is_initial and
                self.is_final == other.is_final and
                self.final_word_out == other.final_word_out and
                self.word_out == other.word_out and
                color and
                self.initial_probability == other.initial_probability)

    def __bool__(self):
        """
        Return ``True``.

        TESTS::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: bool(FSMState('A'))
            True
        """
        return True  # A state cannot be zero (see __init__)

    def _epsilon_successors_(self, fsm=None):
        """
        Return the dictionary with states reachable from ``self``
        without reading anything from an input tape as keys. The
        values are lists of outputs.

        INPUT:

        - ``fsm`` -- the finite state machine to which ``self``
          belongs

        OUTPUT: a dictionary mapping states to a list of output words

        The states in the output are the epsilon successors of
        ``self``. Each word of the list of words is an output word
        written when taking a path from ``self`` to the corresponding
        state.

        TESTS::

            sage: T = Transducer([(0, 1, None, 'a'), (1, 2, None, 'b')])
            sage: T.state(0)._epsilon_successors_(T)
            {1: [['a']], 2: [['a', 'b']]}
            sage: T.state(1)._epsilon_successors_(T)
            {2: [['b']]}
            sage: T.state(2)._epsilon_successors_(T)
            {}

        ::

            sage: T.state(0)._epsilon_successors_()
            {1: [['a']], 2: [['a', 'b']]}

        ::

            sage: T.add_transition(2, 0, None, 'c')
            Transition from 2 to 0: -|'c'
            sage: T.state(0)._epsilon_successors_()
            {0: [['a', 'b', 'c']], 1: [['a']], 2: [['a', 'b']]}

        ::

            sage: T.add_transition(0, 2, None, ['a', 'b'])
            Transition from 0 to 2: -|'a','b'
            sage: T.state(0)._epsilon_successors_()
            {0: [['a', 'b', 'c']], 1: [['a']], 2: [['a', 'b']]}
        """
        if not hasattr(self, 'transitions'):
            raise ValueError('State %s does not belong to a '
                             'finite state machine.' % (self,))

        it = _FSMProcessIteratorEpsilon_(fsm, input_tape=[],
                                         initial_state=self)
        # TODO: optimize the following lines (use already calculated
        # epsilon successors)
        for _ in it:
            pass
        _epsilon_successors_dict_ = it.visited_states
        _epsilon_successors_dict_[self].remove([])  # delete starting state
        if not _epsilon_successors_dict_[self]:
            del _epsilon_successors_dict_[self]
        for s, outputs in _epsilon_successors_dict_.items():
            _epsilon_successors_dict_[s] = [t for t, _ in
                                            itertools.groupby(sorted(outputs))]
        return _epsilon_successors_dict_

    def _in_epsilon_cycle_(self, fsm=None):
        """
        Return whether ``self`` is in an epsilon-cycle or not.

        INPUT:

        - ``fsm`` -- the finite state machine to which ``self``
          belongs

        OUTPUT: boolean

        TESTS::

            sage: A = Automaton([(0, 1, None, 'a'), (1, 2, None, 'b'),
            ....:                (2, 0, None, 'c'), (4, 1, None, 'd')])
            sage: A.state(0)._epsilon_successors_(A)
            {0: [['a', 'b', 'c']], 1: [['a']], 2: [['a', 'b']]}
            sage: A.state(0)._in_epsilon_cycle_(A)
            True
            sage: A.state(4)._epsilon_successors_(A)
            {0: [['d', 'b', 'c']], 1: [['d'], ['d', 'b', 'c', 'a']],
             2: [['d', 'b']]}
            sage: A.state(4)._in_epsilon_cycle_(A)
            False
        """
        return self in self._epsilon_successors_(fsm)

    def _epsilon_cycle_output_empty_(self, fsm=None):
        """
        Return whether all epsilon-cycles in which ``self`` is
        contained have an empty output (i.e., do not write any output
        word).

        INPUT:

        - ``fsm`` -- the finite state machine to which ``self``
          belongs

        OUTPUT: boolean

        A :exc:`ValueError` is raised when ``self`` is not in an epsilon
        cycle.

        TESTS::

            sage: A = Automaton([(0, 1, None, 'a'), (1, 2, None, None),
            ....:                (2, 0, None, None), (4, 1, None, None)])
            sage: A.state(0)._epsilon_successors_(A)
            {0: [['a']], 1: [['a']], 2: [['a']]}
            sage: A.state(0)._epsilon_cycle_output_empty_(A)
            False
            sage: A.state(4)._epsilon_cycle_output_empty_(A)
            Traceback (most recent call last):
            ...
            ValueError: State 4 is not in an epsilon cycle.
            sage: A = Automaton([(0, 1, None, None), (1, 2, None, None),
            ....:                (2, 0, None, None), (4, 1, None, None)])
            sage: A.state(0)._epsilon_successors_(A)
            {0: [[]], 1: [[]], 2: [[]]}
            sage: A.state(0)._epsilon_cycle_output_empty_(A)
            True
            sage: A.process([], initial_state=A.state(0))
            [(False, 0), (False, 1), (False, 2)]
            sage: A.add_transition(0, 0, None, 'x')
            Transition from 0 to 0: -|'x'
            sage: A.state(0)._epsilon_successors_(A)
            {0: [[], ['x']], 1: [[]], 2: [[]]}
            sage: A.state(0)._epsilon_cycle_output_empty_(A)
            False
            sage: A.process([], initial_state=A.state(0))
            Traceback (most recent call last):
            ...
            RuntimeError: State 0 is in an epsilon cycle (no input),
            but output is written.
            sage: T = Transducer([(0, 1, None, None), (1, 2, None, None),
            ....:                 (2, 0, None, None), (0, 0, None, None)])
            sage: T.state(0)._epsilon_successors_(T)
            {0: [[]], 1: [[]], 2: [[]]}
            sage: T.state(0)._epsilon_cycle_output_empty_(T)
            True
        """
        try:
            return not any(self._epsilon_successors_(fsm)[self])
        except KeyError:
            raise ValueError("State %s is not in an epsilon cycle." % (self,))


# ****************************************************************************


def is_FSMTransition(T):
    """
    Test whether or not ``T`` inherits from :class:`FSMTransition`.

    TESTS::

        sage: from sage.combinat.finite_state_machine import is_FSMTransition, FSMTransition
        sage: is_FSMTransition(FSMTransition('A', 'B'))
        doctest:warning...
        DeprecationWarning: The function is_FSMTransition is deprecated; use 'isinstance(..., FSMTransition)' instead.
        See https://github.com/sagemath/sage/issues/38032 for details.
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(38032, "The function is_FSMTransition is deprecated; use 'isinstance(..., FSMTransition)' instead.")
    return isinstance(T, FSMTransition)


class FSMTransition(SageObject):
    """
    Class for a transition of a finite state machine.

    INPUT:

    - ``from_state`` -- state from which transition starts

    - ``to_state`` -- state in which transition ends

    - ``word_in`` -- the input word of the transitions (when the
      finite state machine is used as automaton)

    - ``word_out`` -- the output word of the transitions (when the
      finite state machine is used as transducer)

    OUTPUT: a transition of a finite state machine

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

    from_state = None
    """State from which the transition starts. Read-only."""

    to_state = None
    """State in which the transition ends. Read-only."""

    word_in = None
    """Input word of the transition. Read-only."""

    word_out = None
    """Output word of the transition. Read-only."""

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
        if isinstance(from_state, FSMState):
            self.from_state = from_state
        else:
            self.from_state = FSMState(from_state)
        if isinstance(to_state, FSMState):
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
            if callable(hook):
                self.hook = hook
            else:
                raise TypeError('Wrong argument for hook.')

    def __lt__(self, other):
        """
        Return ``True`` if ``self`` is less than ``other`` with respect to the
        key ``(self.from_state, self.word_in, self.to_state, self.word_out)``.

        INPUT:

        - ``other`` -- a transition

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: FSMTransition(0,1,0,0) < FSMTransition(1,0,0,0)
            True
        """
        return (self.from_state, self.word_in, self.to_state, self.word_out) < \
            (other.from_state, other.word_in, other.to_state, other.word_out)

    def __copy__(self):
        """
        Return a (shallow) copy of the transition.

        OUTPUT: a new transition

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
        Return a deep copy of the transition.

        INPUT:

        - ``memo`` -- dictionary storing already processed elements

        OUTPUT: a new transition

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
        Return a deep copy of the transition.

        INPUT:

        - ``memo`` -- (default: ``None``) a dictionary storing already
          processed elements

        OUTPUT: a new transition

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: t = FSMTransition('A', 'B', 0)
            sage: deepcopy(t)
            Transition from 'A' to 'B': 0|-
        """
        return deepcopy(self, memo)

    __hash__ = None
    # Since transitions are mutable, they should not be hashable

    def _repr_(self):
        """
        Represent a transitions as from state to state and input, output.

        OUTPUT: string

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
        Return the input and output of a transition as "word_in|word_out".

        OUTPUT: string of the input and output labels

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: FSMTransition('A', 'B', 0, 1)._in_out_label_()
            '0|1'
        """
        return "%s|%s" % (FSMWordSymbol(self.word_in),
                          FSMWordSymbol(self.word_out))

    def __eq__(self, other):
        """
        Return ``True`` if the two transitions are the same, i.e., if the
        both go from the same states to the same states and read and
        write the same words.

        Note that the hooks are not checked.

        INPUT:

        - ``self`` -- a transition

        - ``other`` -- a transition

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState, FSMTransition
            sage: A = FSMState('A', is_initial=True)
            sage: t1 = FSMTransition('A', 'B', 0, 1)
            sage: t2 = FSMTransition(A, 'B', 0, 1)
            sage: t1 == t2
            True
        """
        if not isinstance(other, FSMTransition):
            return False
        return self.from_state == other.from_state \
            and self.to_state == other.to_state \
            and self.word_in == other.word_in \
            and self.word_out == other.word_out

    def __ne__(self, other):
        """
        Test for inequality, complement of __eq__.

        INPUT:

        - ``self`` -- a transition

        - ``other`` -- a transition

        OUTPUT: boolean

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState, FSMTransition
            sage: A = FSMState('A', is_initial=True)
            sage: t1 = FSMTransition('A', 'B', 0, 1)
            sage: t2 = FSMTransition(A, 'B', 0, 1)
            sage: t1 != t2
            False
        """
        return not (self == other)

    def __bool__(self):
        """
        Return ``True``.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: bool(FSMTransition('A', 'B', 0))
            True
        """
        return True  # A transition cannot be zero (see __init__)


# ****************************************************************************


def is_FiniteStateMachine(FSM):
    """
    Test whether or not ``FSM`` inherits from :class:`FiniteStateMachine`.

    TESTS::

        sage: from sage.combinat.finite_state_machine import is_FiniteStateMachine
        sage: is_FiniteStateMachine(FiniteStateMachine())
        doctest:warning...
        DeprecationWarning: The function is_FiniteStateMachine is deprecated; use 'isinstance(..., FiniteStateMachine)' instead.
        See https://github.com/sagemath/sage/issues/38032 for details.
        True
        sage: is_FiniteStateMachine(Automaton())
        True
        sage: is_FiniteStateMachine(Transducer())
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(38032, "The function is_FiniteStateMachine is deprecated; use 'isinstance(..., FiniteStateMachine)' instead.")
    return isinstance(FSM, FiniteStateMachine)


def duplicate_transition_ignore(old_transition, new_transition):
    """
    Default function for handling duplicate transitions in finite
    state machines. This implementation ignores the occurrence.

    See the documentation of the ``on_duplicate_transition`` parameter
    of :class:`FiniteStateMachine`.

    INPUT:

    - ``old_transition`` -- a transition in a finite state machine

    - ``new_transition`` -- a transition, identical to ``old_transition``,
      which is to be inserted into the finite state machine

    OUTPUT: the same transition, unchanged

    EXAMPLES::

        sage: from sage.combinat.finite_state_machine import duplicate_transition_ignore
        sage: from sage.combinat.finite_state_machine import FSMTransition
        sage: duplicate_transition_ignore(FSMTransition(0, 0, 1),
        ....:                             FSMTransition(0, 0, 1))
        Transition from 0 to 0: 1|-
    """
    return old_transition


def duplicate_transition_raise_error(old_transition, new_transition):
    """
    Alternative function for handling duplicate transitions in finite
    state machines.

    This implementation raises a :exc:`ValueError`.

    See the documentation of the ``on_duplicate_transition`` parameter
    of :class:`FiniteStateMachine`.

    INPUT:

    - ``old_transition`` -- a transition in a finite state machine

    - ``new_transition`` -- a transition, identical to ``old_transition``,
      which is to be inserted into the finite state machine

    OUTPUT: nothing. A :exc:`ValueError` is raised

    EXAMPLES::

        sage: from sage.combinat.finite_state_machine import duplicate_transition_raise_error
        sage: from sage.combinat.finite_state_machine import FSMTransition
        sage: duplicate_transition_raise_error(FSMTransition(0, 0, 1),
        ....:                                  FSMTransition(0, 0, 1))
        Traceback (most recent call last):
        ...
        ValueError: Attempting to re-insert transition Transition from 0 to 0: 1|-
    """
    raise ValueError("Attempting to re-insert transition %s" % old_transition)


def duplicate_transition_add_input(old_transition, new_transition):
    """
    Alternative function for handling duplicate transitions in finite
    state machines. This implementation adds the input label of the
    new transition to the input label of the old transition.  This is
    intended for the case where a Markov chain is modelled by a finite
    state machine using the input labels as transition probabilities.

    See the documentation of the ``on_duplicate_transition`` parameter
    of :class:`FiniteStateMachine`.

    INPUT:

    - ``old_transition`` -- a transition in a finite state machine

    - ``new_transition`` -- a transition, identical to ``old_transition``,
      which is to be inserted into the finite state machine

    OUTPUT:

    A transition whose input weight is the sum of the input
    weights of ``old_transition`` and ``new_transition``.

    EXAMPLES::

        sage: from sage.combinat.finite_state_machine import duplicate_transition_add_input
        sage: from sage.combinat.finite_state_machine import FSMTransition
        sage: duplicate_transition_add_input(FSMTransition('a', 'a', 1/2),
        ....:                                FSMTransition('a', 'a', 1/2))
        Transition from 'a' to 'a': 1|-

    Input labels must be lists of length 1::

        sage: duplicate_transition_add_input(FSMTransition('a', 'a', [1, 1]),
        ....:                                FSMTransition('a', 'a', [1, 1]))
        Traceback (most recent call last):
        ...
        TypeError: Trying to use duplicate_transition_add_input on
        "Transition from 'a' to 'a': 1,1|-" and
        "Transition from 'a' to 'a': 1,1|-",
        but input words are assumed to be lists of length 1
    """
    if (isinstance(old_transition.word_in, Iterable)
        and len(old_transition.word_in) == 1
        and isinstance(new_transition.word_in, Iterable)
        and len(new_transition.word_in) == 1):
        old_transition.word_in = [old_transition.word_in[0]
                                  + new_transition.word_in[0]]
    else:
        raise TypeError('Trying to use duplicate_transition_add_input on ' +
                        '"%s" and "%s", ' % (old_transition, new_transition) +
                        'but input words are assumed to be lists of length 1')
    return old_transition


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

    - ``determine_alphabets`` -- if ``True``, then the function
      :meth:`.determine_alphabets` is called after ``data`` was read and
      processed, if ``False``, then not. If it is ``None``, then it is
      decided during the construction of the finite state machine
      whether :meth:`.determine_alphabets` should be called.

    - ``with_final_word_out`` -- if given (not ``None``), then the
      function :meth:`.with_final_word_out` (more precisely, its inplace
      pendant :meth:`.construct_final_word_out`) is called with input
      ``letters=with_final_word_out`` at the end of the creation
      process.

    - ``store_states_dict`` -- if ``True``, then additionally the states
      are stored in an internal dictionary for speed up

    - ``on_duplicate_transition`` -- a function which is called when a
      transition is inserted into ``self`` which already existed (same
      ``from_state``, same ``to_state``, same ``word_in``, same ``word_out``).

      This function is assumed to take two arguments, the first being
      the already existing transition, the second being the new
      transition (as an :class:`FSMTransition`). The function must
      return the (possibly modified) original transition.

      By default, we have ``on_duplicate_transition=None``, which is
      interpreted as
      ``on_duplicate_transition=duplicate_transition_ignore``, where
      ``duplicate_transition_ignore`` is a predefined function
      ignoring the occurrence. Other such predefined functions are
      ``duplicate_transition_raise_error`` and
      ``duplicate_transition_add_input``.

    OUTPUT: a finite state machine

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

        If you want to use the arguments of :class:`FSMTransition`
        directly, you can use a dictionary::

            sage: FiniteStateMachine({'a':{'b':{'word_in':0, 'word_out':1},
            ....:                          'c':{'word_in':1, 'word_out':1}}})
            Finite state machine with 3 states

        In the case you already have instances of
        :class:`FSMTransition`, it is possible to use them directly::

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
        :exc:`LookupError` or return an empty list.

        When constructing a finite state machine in this way, some
        initial states and an input alphabet have to be specified.

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

            sage: F = FiniteStateMachine()
            sage: G = Transducer(F)
            sage: G == F
            True

        The other parameters cannot be specified in that case. If you
        want to change these, use the attributes
        :attr:`FSMState.is_initial`, :attr:`FSMState.is_final`,
        :attr:`input_alphabet`, :attr:`output_alphabet`,
        :attr:`on_duplicate_transition` and methods
        :meth:`.determine_alphabets`,
        :meth:`.construct_final_word_out` on the new machine,
        respectively.

    The following examples demonstrate the use of ``on_duplicate_transition``::

        sage: F = FiniteStateMachine([['a', 'a', 1/2], ['a', 'a', 1/2]])
        sage: F.transitions()
        [Transition from 'a' to 'a': 1/2|-]

    ::

        sage: from sage.combinat.finite_state_machine import duplicate_transition_raise_error
        sage: F1 = FiniteStateMachine([['a', 'a', 1/2], ['a', 'a', 1/2]],
        ....:                         on_duplicate_transition=duplicate_transition_raise_error)
        Traceback (most recent call last):
        ...
        ValueError: Attempting to re-insert transition Transition from 'a' to 'a': 1/2|-

    Use ``duplicate_transition_add_input`` to emulate a Markov chain,
    the input labels are considered as transition probabilities::

        sage: from sage.combinat.finite_state_machine import duplicate_transition_add_input
        sage: F = FiniteStateMachine([['a', 'a', 1/2], ['a', 'a', 1/2]],
        ....:                        on_duplicate_transition=duplicate_transition_add_input)
        sage: F.transitions()
        [Transition from 'a' to 'a': 1|-]

    Use ``with_final_word_out`` to construct final output::

        sage: T = Transducer([(0, 1, 0, 0), (1, 0, 0, 0)],
        ....:                initial_states=[0],
        ....:                final_states=[0],
        ....:                with_final_word_out=0)
        sage: for s in T.iter_final_states():
        ....:     print("{} {}".format(s, s.final_word_out))
        0 []
        1 [0]

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

    We test that no input parameter is allowed when creating a finite
    state machine from an existing instance::

        sage: F = FiniteStateMachine()
        sage: FiniteStateMachine(F, initial_states=[1])
        Traceback (most recent call last):
        ...
        ValueError: initial_states cannot be specified when
        copying another finite state machine.
        sage: FiniteStateMachine(F, final_states=[1])
        Traceback (most recent call last):
        ...
        ValueError: final_states cannot be specified when
        copying another finite state machine.
        sage: FiniteStateMachine(F, input_alphabet=[1])
        Traceback (most recent call last):
        ...
        ValueError: input_alphabet cannot be specified when
        copying another finite state machine.
        sage: FiniteStateMachine(F, output_alphabet=[1])
        Traceback (most recent call last):
        ...
        ValueError: output_alphabet cannot be specified when
        copying another finite state machine.
        sage: from sage.combinat.finite_state_machine import (
        ....:     duplicate_transition_add_input)
        sage: FiniteStateMachine(F,
        ....:     on_duplicate_transition=duplicate_transition_add_input)
        Traceback (most recent call last):
        ...
        ValueError: on_duplicate_transition cannot be specified when
        copying another finite state machine.
        sage: FiniteStateMachine(F, determine_alphabets=False)
        Traceback (most recent call last):
        ...
        ValueError: determine_alphabets cannot be specified when
        copying another finite state machine.
        sage: FiniteStateMachine(F, with_final_word_out=[1])
        Traceback (most recent call last):
        ...
        ValueError: with_final_word_out cannot be specified when
        copying another finite state machine.

    :issue:`19454` rewrote automatic detection of the alphabets::

        sage: def transition_function(state, letter):
        ....:     return (0, 3 + letter)
        sage: T1 = Transducer(transition_function,
        ....:     input_alphabet=[0, 1],
        ....:     initial_states=[0],
        ....:     final_states=[0])
        sage: T1.output_alphabet
        [3, 4]
        sage: T2 = Transducer([(0, 0, 0, 3), (0, 0, 0, 4)],
        ....:     initial_states=[0],
        ....:     final_states=[0])
        sage: T2.output_alphabet
        [3, 4]
        sage: T = Transducer([(0, 0, 1, 2)])
        sage: (T.input_alphabet, T.output_alphabet)
        ([1], [2])
        sage: T = Transducer([(0, 0, 1, 2)], determine_alphabets=False)
        sage: (T.input_alphabet, T.output_alphabet)
        (None, None)
        sage: T = Transducer([(0, 0, 1, 2)], input_alphabet=[0, 1])
        sage: (T.input_alphabet, T.output_alphabet)
        ([0, 1], [2])
        sage: T = Transducer([(0, 0, 1, 2)], output_alphabet=[2, 3])
        sage: (T.input_alphabet, T.output_alphabet)
        ([1], [2, 3])
        sage: T = Transducer([(0, 0, 1, 2)], input_alphabet=[0, 1],
        ....:     output_alphabet=[2, 3])
        sage: (T.input_alphabet, T.output_alphabet)
        ([0, 1], [2, 3])

    .. automethod:: __call__
    """

    on_duplicate_transition = duplicate_transition_ignore
    """
    Which function to call when a duplicate transition is inserted.

    It can be set by the parameter ``on_duplicate_transition`` when
    initializing a finite state machine, see
    :class:`FiniteStateMachine`.

    .. SEEALSO::

        :class:`FiniteStateMachine`, :meth:`is_Markov_chain`,
        :meth:`markov_chain_simplification`.
    """

    input_alphabet = None
    """
    A list of letters representing the input alphabet of the finite
    state machine.

    It can be set by the parameter ``input_alphabet`` when initializing
    a finite state machine, see :class:`FiniteStateMachine`.

    It can also be set by the method :meth:`determine_alphabets`.

    .. SEEALSO::

        :class:`FiniteStateMachine`, :meth:`determine_alphabets`,
        :attr:`output_alphabet`.
    """

    output_alphabet = None
    """
    A list of letters representing the output alphabet of the finite
    state machine.

    It can be set by the parameter ``output_alphabet`` when initializing
    a finite state machine, see :class:`FiniteStateMachine`.

    It can also be set by the method :meth:`determine_alphabets`.

    .. SEEALSO::

        :class:`FiniteStateMachine`,
        :meth:`determine_alphabets`,
        :attr:`input_alphabet`.
    """

    # ************************************************************************
    # init
    # ************************************************************************

    def __init__(self,
                 data=None,
                 initial_states=None, final_states=None,
                 input_alphabet=None, output_alphabet=None,
                 determine_alphabets=None,
                 with_final_word_out=None,
                 store_states_dict=True,
                 on_duplicate_transition=None):
        """
        See :class:`FiniteStateMachine` for more information.

        TESTS::

            sage: FiniteStateMachine()
            Empty finite state machine
        """
        self._states_ = []
        # List of states in the finite state
        # machine.  Each state stores a list of
        # outgoing transitions.
        if store_states_dict:
            self._states_dict_ = {}

        self._allow_composition_ = True

        if isinstance(data, FiniteStateMachine):
            if initial_states is not None:
                raise ValueError(
                    "initial_states cannot be specified when copying "
                    "another finite state machine.")
            if final_states is not None:
                raise ValueError(
                    "final_states cannot be specified when copying "
                    "another finite state machine.")
            if input_alphabet is not None:
                raise ValueError(
                    "input_alphabet cannot be specified when copying "
                    "another finite state machine.")
            if output_alphabet is not None:
                raise ValueError(
                    "output_alphabet cannot be specified when copying "
                    "another finite state machine.")
            if on_duplicate_transition is not None:
                raise ValueError(
                    "on_duplicate_transition cannot be specified when "
                    "copying another finite state machine.")
            if determine_alphabets is not None:
                raise ValueError(
                    "determine_alphabets cannot be specified when "
                    "copying another finite state machine.")
            if with_final_word_out is not None:
                raise ValueError(
                    "with_final_word_out cannot be specified when "
                    "copying another finite state machine.")

            self._copy_from_other_(data)
            return

        if initial_states is not None:
            if not isinstance(initial_states, Iterable):
                raise TypeError('Initial states must be iterable '
                                '(e.g. a list of states).')
            for s in initial_states:
                state = self.add_state(s)
                state.is_initial = True

        if final_states is not None:
            if not isinstance(final_states, Iterable):
                raise TypeError('Final states must be iterable '
                                '(e.g. a list of states).')
            for s in final_states:
                state = self.add_state(s)
                state.is_final = True

        self.input_alphabet = input_alphabet
        self.output_alphabet = output_alphabet

        if on_duplicate_transition is None:
            on_duplicate_transition = duplicate_transition_ignore
        if callable(on_duplicate_transition):
            self.on_duplicate_transition = on_duplicate_transition
        else:
            raise TypeError('on_duplicate_transition must be callable')

        if data is None:
            pass
        elif isinstance(data, Mapping):
            # data is a dict (or something similar),
            # format: key = from_state, value = iterator of transitions
            for (sf, iter_transitions) in data.items():
                self.add_state(sf)
                if isinstance(iter_transitions, Mapping):
                    for (st, transition) in iter_transitions.items():
                        self.add_state(st)
                        if isinstance(transition, FSMTransition):
                            self.add_transition(transition)
                        elif isinstance(transition, Mapping):
                            self.add_transition(sf, st, **transition)
                        elif isinstance(transition, Iterable):
                            self.add_transition(sf, st, *transition)
                        else:
                            self.add_transition(sf, st, transition)
                elif isinstance(iter_transitions, Iterable):
                    for transition in iter_transitions:
                        if isinstance(transition, Iterable):
                            L = [sf]
                            L.extend(transition)
                        elif isinstance(transition, FSMTransition):
                            L = transition
                        else:
                            L = [sf, transition]
                        self.add_transition(L)
                else:
                    raise TypeError('Wrong input data for transition.')
        elif isinstance(data, Iterable):
            # data is a something that is iterable,
            # items are transitions
            for transition in data:
                if isinstance(transition, FSMTransition):
                    self.add_transition(transition)
                elif isinstance(transition, Mapping):
                    self.add_transition(transition)
                elif isinstance(transition, Iterable):
                    self.add_transition(transition)
                else:
                    raise TypeError('Wrong input data for transition.')
        elif callable(data):
            self.add_from_transition_function(data)
        else:
            raise TypeError('Cannot decide what to do with data.')

        if determine_alphabets is None and data is not None:
            if input_alphabet is None:
                self.determine_input_alphabet()
            if output_alphabet is None:
                self.determine_output_alphabet()
        elif determine_alphabets:
            self.determine_alphabets()

        if with_final_word_out is not None:
            self.construct_final_word_out(with_final_word_out)

    # ************************************************************************
    # copy and hash
    # ************************************************************************

    def __copy__(self):
        """
        Return a (shallow) copy of the finite state machine.

        OUTPUT: a new finite state machine

        TESTS::

            sage: copy(FiniteStateMachine())
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    copy = __copy__

    def empty_copy(self, memo=None, new_class=None):
        """
        Return an empty deep copy of the finite state machine, i.e.,
        ``input_alphabet``, ``output_alphabet``, ``on_duplicate_transition``
        are preserved, but states and transitions are not.

        INPUT:

        - ``memo`` -- dictionary storing already processed elements

        - ``new_class`` -- a class for the copy; by default
          (``None``), the class of ``self`` is used

        OUTPUT: a new finite state machine

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import duplicate_transition_raise_error
            sage: F = FiniteStateMachine([('A', 'A', 0, 2), ('A', 'A', 1, 3)],
            ....:                        input_alphabet=[0, 1],
            ....:                        output_alphabet=[2, 3],
            ....:                        on_duplicate_transition=duplicate_transition_raise_error)
            sage: FE = F.empty_copy(); FE
            Empty finite state machine
            sage: FE.input_alphabet
            [0, 1]
            sage: FE.output_alphabet
            [2, 3]
            sage: FE.on_duplicate_transition == duplicate_transition_raise_error
            True

        TESTS::

            sage: T = Transducer()
            sage: type(T.empty_copy())
            <class 'sage.combinat.finite_state_machine.Transducer'>
            sage: type(T.empty_copy(new_class=Automaton))
            <class 'sage.combinat.finite_state_machine.Automaton'>
        """
        if new_class is None:
            new = self.__class__()
        else:
            new = new_class()
        new._copy_from_other_(self, memo=memo, empty=True)
        return new

    def __deepcopy__(self, memo):
        """
        Return a deep copy of the finite state machine.

        INPUT:

        - ``memo`` -- dictionary storing already processed elements

        OUTPUT: a new finite state machine

        EXAMPLES::

            sage: F = FiniteStateMachine([('A', 'A', 0, 1), ('A', 'A', 1, 0)])
            sage: deepcopy(F)
            Finite state machine with 1 state
        """
        new = self.__class__()
        new._copy_from_other_(self)
        return new

    def deepcopy(self, memo=None):
        """
        Return a deep copy of the finite state machine.

        INPUT:

        - ``memo`` -- (default: ``None``) a dictionary storing already
          processed elements

        OUTPUT: a new finite state machine

        EXAMPLES::

            sage: F = FiniteStateMachine([('A', 'A', 0, 1), ('A', 'A', 1, 0)])
            sage: deepcopy(F)
            Finite state machine with 1 state

        TESTS:

        Make sure that the links between transitions and states
        are still intact::

            sage: C = deepcopy(F)
            sage: C.transitions()[0].from_state is C.state('A')
            True
            sage: C.transitions()[0].to_state is C.state('A')
            True
        """
        return deepcopy(self, memo)

    def _copy_from_other_(self, other, memo=None, empty=False):
        """
        Copy all data from other to self, to be used in the constructor.

        INPUT:

        - ``other`` -- a :class:`FiniteStateMachine`

        OUTPUT: nothing

        EXAMPLES::

            sage: A = Automaton([(0, 0, 0)],
            ....:               initial_states=[0],
            ....:               final_states=[0])
            sage: B = Automaton()
            sage: B._copy_from_other_(A)
            sage: A == B
            True
        """
        if memo is None:
            memo = {}
        self.input_alphabet = deepcopy(other.input_alphabet, memo)
        self.output_alphabet = deepcopy(other.output_alphabet, memo)
        self.on_duplicate_transition = other.on_duplicate_transition

        if not empty:
            relabel = hasattr(other, '_deepcopy_relabel_')
            relabel_iter = itertools.count(0)
            for state in other.iter_states():
                if relabel:
                    if other._deepcopy_labels_ is None:
                        state._deepcopy_relabel_ = next(relabel_iter)
                    elif callable(other._deepcopy_labels_):
                        state._deepcopy_relabel_ = \
                            other._deepcopy_labels_(state.label())
                    elif hasattr(other._deepcopy_labels_, '__getitem__'):
                        state._deepcopy_relabel_ = \
                            other._deepcopy_labels_[state.label()]
                    else:
                        raise TypeError("labels must be None, a callable "
                                        "or a dictionary.")
                s = deepcopy(state, memo)
                if relabel:
                    del state._deepcopy_relabel_
                self.add_state(s)
            for transition in other.iter_transitions():
                self.add_transition(deepcopy(transition, memo))

    def __getstate__(self):
        """
        Return state for pickling excluding outgoing transitions.

        OUTPUT: a dictionary

        Outgoing transitions are in fact stored in states,
        but must be pickled by the finite state machine
        in order to avoid deep recursion.

        TESTS::

            sage: A = Automaton([(0, 1, 0)])
            sage: loads(dumps(A)) == A
            True
        """
        odict = self.__dict__.copy()  # copy the dict since we change it
        odict.update({'transitions': self.transitions()})
        return odict

    def __setstate__(self, d):
        """
        Set state from pickling.

        INPUT:

        - ``d`` -- dictionary

        OUTPUT: none

        As transitions are in fact stored in states but not saved
        by states in order to avoid deep recursion, transitions
        have to be manually restored here.

        TESTS::

            sage: A = Automaton([(0, 1, 0)])
            sage: loads(dumps(A)) == A
            True
        """
        transitions = d.pop('transitions')
        self.__dict__.update(d)
        for state in self.iter_states():
            state.transitions = []  # clean outgoing transitions
        for transition in transitions:
            self.add_transition(transition)

    def relabeled(self, memo=None, labels=None):
        """
        Return a deep copy of the finite state machine, but the
        states are relabeled.

        INPUT:

        - ``memo`` -- (default: ``None``) a dictionary storing already
          processed elements

        - ``labels`` -- (default: ``None``) a dictionary or callable
          mapping old labels to new labels. If ``None``, then the new
          labels are integers starting with 0.

        OUTPUT: a new finite state machine

        EXAMPLES::

            sage: FSM1 = FiniteStateMachine([('A', 'B'), ('B', 'C'), ('C', 'A')])
            sage: FSM1.states()
            ['A', 'B', 'C']
            sage: FSM2 = FSM1.relabeled()
            sage: FSM2.states()
            [0, 1, 2]
            sage: FSM3 = FSM1.relabeled(labels={'A': 'a', 'B': 'b', 'C': 'c'})
            sage: FSM3.states()
            ['a', 'b', 'c']
            sage: FSM4 = FSM2.relabeled(labels=lambda x: 2*x)
            sage: FSM4.states()
            [0, 2, 4]

        TESTS::

            sage: FSM2.relabeled(labels=1)
            Traceback (most recent call last):
            ...
            TypeError: labels must be None, a callable or a dictionary.
        """
        self._deepcopy_relabel_ = True
        self._deepcopy_labels_ = labels
        new = deepcopy(self, memo)
        del self._deepcopy_relabel_
        del self._deepcopy_labels_
        return new

    def induced_sub_finite_state_machine(self, states):
        """
        Return a sub-finite-state-machine of the finite state machine
        induced by the given states.

        INPUT:

        - ``states`` -- list (or an iterator) of states (either labels or
          instances of :class:`FSMState`) of the sub-finite-state-machine

        OUTPUT:

        A new finite state machine. It consists (of deep copies) of
        the given states and (deep copies) of all transitions of ``self``
        between these states.

        EXAMPLES::

            sage: FSM = FiniteStateMachine([(0, 1, 0), (0, 2, 0),
            ....:                           (1, 2, 0), (2, 0, 0)])
            sage: sub_FSM = FSM.induced_sub_finite_state_machine([0, 1])
            sage: sub_FSM.states()
            [0, 1]
            sage: sub_FSM.transitions()
            [Transition from 0 to 1: 0|-]
            sage: FSM.induced_sub_finite_state_machine([3])
            Traceback (most recent call last):
            ...
            ValueError: 3 is not a state of this finite state machine.

        TESTS:

        Make sure that the links between transitions and states
        are still intact::

            sage: sub_FSM.transitions()[0].from_state is sub_FSM.state(0)
            True
        """
        good_states = set()
        for state in states:
            if not self.has_state(state):
                raise ValueError("%s is not a state of this finite state machine." % state)
            good_states.add(self.state(state))

        memo = {}
        new = self.empty_copy(memo=memo)
        for state in good_states:
            s = deepcopy(state, memo)
            new.add_state(s)

        for state in good_states:
            for transition in self.iter_transitions(state):
                if transition.to_state in good_states:
                    new.add_transition(deepcopy(transition, memo))

        return new

    def __hash__(self):
        """
        Since finite state machines are mutable, they should not be
        hashable, so we return a type error.

        OUTPUT: the hash of this finite state machine

        EXAMPLES::

            sage: hash(FiniteStateMachine())
            Traceback (most recent call last):
            ...
            TypeError: Finite state machines are mutable, and thus not hashable.
        """
        if getattr(self, "_immutable", False):
            return hash((tuple(self.states()), tuple(self.transitions())))
        raise TypeError("Finite state machines are mutable, "
                        "and thus not hashable.")

    # ************************************************************************
    # operators
    # ************************************************************************

    def __or__(self, other):
        """
        Return the disjoint union of this and another finite state machine.

        INPUT:

        - ``other`` -- a finite state machine

        OUTPUT: a new finite state machine

        .. SEEALSO::

            :meth:`.disjoint_union`, :meth:`.__and__`,
            :meth:`Automaton.intersection`,
            :meth:`Transducer.intersection`.

        TESTS::

            sage: FiniteStateMachine() | FiniteStateMachine([('A', 'B')])
            Finite state machine with 2 states
            sage: FiniteStateMachine() | 42
            Traceback (most recent call last):
            ...
            TypeError: Can only add finite state machine
        """
        if isinstance(other, FiniteStateMachine):
            return self.disjoint_union(other)
        else:
            raise TypeError("Can only add finite state machine")

    __add__ = __or__

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

    def __and__(self, other):
        """
        Return the intersection of ``self`` with ``other``.

        TESTS::

            sage: FiniteStateMachine() & FiniteStateMachine([('A', 'B')])
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if isinstance(other, FiniteStateMachine):
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
        Call either method :meth:`.composition` or :meth:`.process`
        (with ``full_output=False``). If the input is not finite
        (``is_finite`` of input is ``False``), then
        :meth:`.iter_process` (with ``iterator_type='simple'``) is
        called. Moreover, the flag ``automatic_output_type`` is set
        (unless ``format_output`` is specified).
        See the documentation of these functions for possible
        parameters.

        EXAMPLES:

        The following code performs a :meth:`composition`::

            sage: F = Transducer([('A', 'B', 1, 0), ('B', 'B', 1, 1),
            ....:                 ('B', 'B', 0, 0)],
            ....:                initial_states=['A'], final_states=['B'])
            sage: G = Transducer([(1, 1, 0, 0), (1, 2, 1, 0),
            ....:                 (2, 2, 0, 1), (2, 1, 1, 1)],
            ....:                initial_states=[1], final_states=[1])
            sage: H = G(F)
            sage: H.states()
            [('A', 1), ('B', 1), ('B', 2)]

        An automaton or transducer can also act on an input (a list
        or other iterable of letters)::

            sage: binary_inverter = Transducer({'A': [('A', 0, 1), ('A', 1, 0)]},
            ....:                              initial_states=['A'], final_states=['A'])
            sage: binary_inverter([0, 1, 0, 0, 1, 1])
            [1, 0, 1, 1, 0, 0]

        We can also let them act on :doc:`words <words/words>`::

            sage: # needs sage.combinat
            sage: W = Words([0, 1]); W
            Finite and infinite words over {0, 1}
            sage: binary_inverter(W([0, 1, 1, 0, 1, 1]))
            word: 100100

        Infinite words work as well::

            sage: # needs sage.combinat
            sage: words.FibonacciWord()
            word: 0100101001001010010100100101001001010010...
            sage: binary_inverter(words.FibonacciWord())
            word: 1011010110110101101011011010110110101101...

        When only one successful path is found in a non-deterministic
        transducer, the result of that path is returned.

        ::

            sage: T = Transducer([(0, 1, 0, 1), (0, 2, 0, 2)],
            ....:                initial_states=[0], final_states=[1])
            sage: T.process([0])
            [(False, 2, [2]), (True, 1, [1])]
            sage: T([0])
            [1]

        .. SEEALSO::

            :meth:`.composition`,
            :meth:`~FiniteStateMachine.process`,
            :meth:`~FiniteStateMachine.iter_process`,
            :meth:`Automaton.process`,
            :meth:`Transducer.process`.

        TESTS::

            sage: F = FiniteStateMachine([(0, 1, 1, 'a'), (0, 2, 2, 'b')],
            ....:                        initial_states=[0],
            ....:                        final_states=[1])
            sage: A = Automaton([(0, 1, 1), (0, 2, 2)],
            ....:               initial_states=[0],
            ....:               final_states=[1])
            sage: T = Transducer([(0, 1, 1, 'a'), (0, 2, 2, 'b')],
            ....:                initial_states=[0],
            ....:                final_states=[1])
            sage: F([1])
            (True, 1, ['a'])
            sage: A([1])
            True
            sage: T([1])
            ['a']
            sage: F([2])
            (False, 2, ['b'])
            sage: A([2])
            False
            sage: T([2])
            Traceback (most recent call last):
            ...
            ValueError: Invalid input sequence.
            sage: F([3])
            (False, None, None)
            sage: A([3])
            False
            sage: T([3])
            Traceback (most recent call last):
            ...
            ValueError: Invalid input sequence.

        ::

            sage: F = FiniteStateMachine([(11, 11, 1, 'a'), (11, 12, 2, 'b'),
            ....:                         (11, 13, 3, 'c'), (11, 14, 4, 'd'),
            ....:                         (12, 13, 3, 'e'), (12, 13, 3, 'f'),
            ....:                         (12, 14, 4, 'g'), (12, 14, 4, 'h'),
            ....:                         (12, 13, 2, 'i'), (12, 14, 2, 'j')],
            ....:                        initial_states=[11],
            ....:                        final_states=[13])
            sage: def f(o):
            ....:     return ''.join(o)
            sage: F([0], format_output=f)
            (False, None, None)
            sage: F([3], format_output=f)
            (True, 13, 'c')
            sage: F([4], format_output=f)
            (False, 14, 'd')
            sage: F([2, 2], format_output=f)
            Traceback (most recent call last):
            ...
            ValueError: Got more than one output, but only allowed to show
            one. Change list_of_outputs option.
            sage: F([2, 2], format_output=f, list_of_outputs=True)
            [(False, 14, 'bj'), (True, 13, 'bi')]
            sage: F([2, 3], format_output=f)
            Traceback (most recent call last):
            ...
            ValueError: Got more than one output, but only allowed to show
            one. Change list_of_outputs option.
            sage: F([2, 3], format_output=f, list_of_outputs=True)
            [(True, 13, 'be'), (True, 13, 'bf')]
            sage: F([2, 4], format_output=f)
            Traceback (most recent call last):
            ...
            ValueError: Got more than one output, but only allowed to show
            one. Change list_of_outputs option.
            sage: F([2, 4], format_output=f, list_of_outputs=True)
            [(False, 14, 'bg'), (False, 14, 'bh')]

        ::

            sage: A = Automaton([(11, 11, 1), (11, 12, 2),
            ....:                (11, 13, 3), (11, 14, 4),
            ....:                (12, 13, 3), (12, 14, 4),
            ....:                (12, 32, 3), (12, 42, 4),
            ....:                (12, 13, 2), (12, 14, 2)],
            ....:               initial_states=[11],
            ....:               final_states=[13, 32])
            sage: def f(o):
            ....:     return ''.join(o)
            sage: A([0], format_output=f)
            False
            sage: A([3], format_output=f)
            True
            sage: A([4], format_output=f)
            False
            sage: A([2, 2], format_output=f)
            True
            sage: A([2, 2], format_output=f, list_of_outputs=True)
            [False, True]
            sage: A([2, 3], format_output=f)
            True
            sage: A([2, 3], format_output=f, list_of_outputs=True)
            [True, True]
            sage: A([2, 4], format_output=f)
            False
            sage: A([2, 4], format_output=f, list_of_outputs=True)
            [False, False]

        ::

            sage: T = Transducer([(11, 11, 1, 'a'), (11, 12, 2, 'b'),
            ....:                 (11, 13, 3, 'c'), (11, 14, 4, 'd'),
            ....:                 (12, 13, 3, 'e'), (12, 13, 3, 'f'),
            ....:                 (12, 14, 4, 'g'), (12, 14, 4, 'h'),
            ....:                 (12, 13, 2, 'i'), (12, 14, 2, 'j')],
            ....:                initial_states=[11],
            ....:                final_states=[13])
            sage: def f(o):
            ....:     return ''.join(o)
            sage: T([0], format_output=f)
            Traceback (most recent call last):
            ...
            ValueError: Invalid input sequence.
            sage: T([3], format_output=f)
            'c'
            sage: T([4], format_output=f)
            Traceback (most recent call last):
            ...
            ValueError: Invalid input sequence.
            sage: T([2, 2], format_output=f)
            'bi'
            sage: T([2, 2], format_output=f, list_of_outputs=True)
            [None, 'bi']
            sage: T([2, 2], format_output=f,
            ....:   list_of_outputs=True, only_accepted=True)
            ['bi']
            sage: T.process([2, 2], format_output=f, list_of_outputs=True)
            [(False, 14, 'bj'), (True, 13, 'bi')]
            sage: T([2, 3], format_output=f)
            Traceback (most recent call last):
            ...
            ValueError: Found more than one accepting path.
            sage: T([2, 3], format_output=f, list_of_outputs=True)
            ['be', 'bf']
            sage: T([2, 4], format_output=f)
            Traceback (most recent call last):
            ...
            ValueError: Invalid input sequence.
            sage: T([2, 4], format_output=f, list_of_outputs=True)
            [None, None]

        ::

            sage: from itertools import islice
            sage: inverter = Transducer({'A': [('A', 0, 1), ('A', 1, 0)]},
            ....:     initial_states=['A'], final_states=['A'])
            sage: inverter(words.FibonacciWord())                                       # needs sage.combinat
            word: 1011010110110101101011011010110110101101...
            sage: inverter(words.FibonacciWord(), automatic_output_type=True)           # needs sage.combinat
            word: 1011010110110101101011011010110110101101...
            sage: tuple(islice(inverter(words.FibonacciWord(),                          # needs sage.combinat
            ....:                       automatic_output_type=False), 10r))
            (1, 0, 1, 1, 0, 1, 0, 1, 1, 0)
            sage: type(inverter((1, 0, 1, 1, 0, 1, 0, 1, 1, 0),
            ....:               automatic_output_type=False))
            <class 'list'>
            sage: type(inverter((1, 0, 1, 1, 0, 1, 0, 1, 1, 0),
            ....:               automatic_output_type=True))
            <class 'tuple'>
        """
        if not args:
            raise TypeError("Called with too few arguments.")
        if isinstance(args[0], FiniteStateMachine):
            return self.composition(*args, **kwargs)
        if isinstance(args[0], Iterable):
            if 'full_output' not in kwargs:
                kwargs['full_output'] = False
            if 'list_of_outputs' not in kwargs:
                kwargs['list_of_outputs'] = False
            if 'automatic_output_type' not in kwargs:
                kwargs['automatic_output_type'] = 'format_output' not in kwargs
            input_tape = args[0]
            if hasattr(input_tape, 'is_finite') and not input_tape.is_finite():
                if 'iterator_type' not in kwargs:
                    kwargs['iterator_type'] = 'simple'
                return self.iter_process(*args, **kwargs)
            return self.process(*args, **kwargs)
        raise TypeError("Do not know what to do with that arguments.")

    # ************************************************************************
    # tests
    # ************************************************************************

    def __bool__(self):
        """
        Return ``True`` if the finite state machine consists of at least
        one state.

        TESTS::

            sage: bool(FiniteStateMachine())
            False
        """
        return bool(self._states_)

    def __eq__(self, other):
        """
        Return ``True`` if the two finite state machines are equal,
        i.e., if they have the same states and the same transitions.

        INPUT:

        - ``self`` -- a finite state machine

        - ``other`` -- a finite state machine

        OUTPUT: boolean

        Note that this function compares all attributes of a state (by
        using :meth:`FSMState.fully_equal`) except for colors. Colors
        are handled as follows: If the colors coincide, then the
        finite state machines are also considered equal. If not, then
        they are considered as equal if both finite state machines are
        monochromatic.

        EXAMPLES::

            sage: F = FiniteStateMachine([('A', 'B', 1)])
            sage: F == FiniteStateMachine()
            False
            sage: G = FiniteStateMachine([('A', 'B', 1)],
            ....:                        initial_states=['A'])
            sage: F == G
            False
            sage: F.state('A').is_initial = True
            sage: F == G
            True

        This shows the behavior when the states have colors::

            sage: F.state('A').color = 'red'
            sage: G.state('A').color = 'red'
            sage: F == G
            True
            sage: G.state('A').color = 'blue'
            sage: F == G
            False
            sage: F.state('B').color = 'red'
            sage: F.is_monochromatic()
            True
            sage: G.state('B').color = 'blue'
            sage: G.is_monochromatic()
            True
            sage: F == G
            True
        """
        if not isinstance(other, FiniteStateMachine):
            return False
        if len(self._states_) != len(other._states_):
            return False
        colors_equal = True
        for state in self.iter_states():
            try:
                other_state = other.state(state.label())
            except LookupError:
                return False

            # we handle colors separately
            if not state.fully_equal(other_state, compare_color=False):
                return False
            if state.color != other_state.color:
                colors_equal = False

            self_transitions = state.transitions
            other_transitions = other.state(state).transitions
            if len(self_transitions) != len(other_transitions):
                return False
            for t in self_transitions:
                if t not in other_transitions:
                    return False

        # handle colors
        if colors_equal:
            return True
        return self.is_monochromatic() and other.is_monochromatic()

    def __ne__(self, other):
        """
        Test for inequality, complement of :meth:`.__eq__`.

        INPUT:

        - ``self`` -- a finite state machine

        - ``other`` -- a finite state machine

        OUTPUT: boolean

        EXAMPLES::

            sage: E = FiniteStateMachine([('A', 'B', 0)])
            sage: F = Automaton([('A', 'B', 0)])
            sage: G = Transducer([('A', 'B', 0, 1)])
            sage: E == F
            True
            sage: E == G
            False
        """
        return not (self == other)

    def __contains__(self, item):
        """
        Return ``True``, if the finite state machine contains the
        state or transition item.

        Note that only the labels of the
        states and the input and output words are tested.

        INPUT:

        - ``item`` -- a state or a transition

        OUTPUT: boolean

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
        if isinstance(item, FSMState):
            return self.has_state(item)
        if isinstance(item, FSMTransition):
            return self.has_transition(item)
        return False

    def is_Markov_chain(self, is_zero=None):
        """
        Check whether ``self`` is a Markov chain where the transition
        probabilities are modeled as input labels.

        INPUT:

        - ``is_zero`` -- by default (``is_zero=None``), checking for
          zero is simply done by
          :meth:`~sage.structure.element.Element.is_zero`.  This
          parameter can be used to provide a more sophisticated check
          for zero, e.g. in the case of symbolic probabilities, see
          the examples below.

        OUTPUT: boolean

        :attr:`on_duplicate_transition` must be
        :func:`duplicate_transition_add_input`, the sum of the input weights
        of the transitions leaving a state must add up to 1 and the sum of
        initial probabilities must add up to 1 (or all be ``None``).

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import duplicate_transition_add_input
            sage: F = Transducer([[0, 0, 1/4, 0], [0, 1, 3/4, 1],
            ....:                 [1, 0, 1/2, 0], [1, 1, 1/2, 1]],
            ....:                on_duplicate_transition=duplicate_transition_add_input)
            sage: F.is_Markov_chain()
            True

        :attr:`on_duplicate_transition` must be
        :func:`duplicate_transition_add_input`::

            sage: F = Transducer([[0, 0, 1/4, 0], [0, 1, 3/4, 1],
            ....:                 [1, 0, 1/2, 0], [1, 1, 1/2, 1]])
            sage: F.is_Markov_chain()
            False

        Sum of input labels of the transitions leaving states must be 1::

            sage: F = Transducer([[0, 0, 1/4, 0], [0, 1, 3/4, 1],
            ....:                 [1, 0, 1/2, 0]],
            ....:                on_duplicate_transition=duplicate_transition_add_input)
            sage: F.is_Markov_chain()
            False

        The initial probabilities of all states must be ``None`` or they must
        sum up to 1. The initial probabilities of all states have to be set in the latter case::

            sage: F = Transducer([[0, 0, 1/4, 0], [0, 1, 3/4, 1],
            ....:                 [1, 0, 1, 0]],
            ....:                on_duplicate_transition=duplicate_transition_add_input)
            sage: F.is_Markov_chain()
            True
            sage: F.state(0).initial_probability = 1/4
            sage: F.is_Markov_chain()
            False
            sage: F.state(1).initial_probability = 7
            sage: F.is_Markov_chain()
            False
            sage: F.state(1).initial_probability = 3/4
            sage: F.is_Markov_chain()
            True

        If the probabilities are variables in the symbolic ring,
        :func:`~sage.symbolic.assumptions.assume` will do the trick::

            sage: # needs sage.symbolic
            sage: var('p q')
            (p, q)
            sage: F = Transducer([(0, 0, p, 1), (0, 0, q, 0)],
            ....:                on_duplicate_transition=duplicate_transition_add_input)
            sage: assume(p + q == 1)
            sage: (p + q - 1).is_zero()
            True
            sage: F.is_Markov_chain()
            True
            sage: forget()
            sage: del(p, q)

        If the probabilities are variables in some polynomial ring,
        the parameter ``is_zero`` can be used::

            sage: R.<p, q> = PolynomialRing(QQ)
            sage: def is_zero_polynomial(polynomial):
            ....:     return polynomial in (p + q - 1)*R
            sage: F = Transducer([(0, 0, p, 1), (0, 0, q, 0)],
            ....:                on_duplicate_transition=duplicate_transition_add_input)
            sage: F.state(0).initial_probability = p + q
            sage: F.is_Markov_chain()
            False
            sage: F.is_Markov_chain(is_zero_polynomial)                                 # needs sage.libs.singular
            True
        """
        def default_is_zero(expression):
            return expression.is_zero()

        is_zero_function = default_is_zero
        if is_zero is not None:
            is_zero_function = is_zero

        if self.on_duplicate_transition != duplicate_transition_add_input:
            return False

        if any(s.initial_probability is not None for s in self.iter_states()) and \
               any(s.initial_probability is None for s in self.iter_states()):
            return False

        if any(s.initial_probability is not None for s in self.iter_states()) and \
               not is_zero_function(sum(s.initial_probability for s
                                        in self.iter_states()) - 1):
            return False

        return all(is_zero_function(sum(t.word_in[0] for t in state.transitions) - 1)
                   for state in self.iter_states())

    # ************************************************************************
    # representations / LaTeX
    # ************************************************************************

    def _repr_(self):
        """
        Represent the finite state machine as "Finite state machine
        with n states" where n is the number of states.

        OUTPUT: string

        EXAMPLES::

            sage: FiniteStateMachine()._repr_()
            'Empty finite state machine'

        TESTS::

            sage: F = FiniteStateMachine()
            sage: F
            Empty finite state machine
            sage: F.add_state(42)
            42
            sage: F
            Finite state machine with 1 state
            sage: F.add_state(43)
            43
            sage: F
            Finite state machine with 2 states
        """
        if not self._states_:
            return "Empty finite state machine"
        if len(self._states_) == 1:
            return "Finite state machine with 1 state"
        else:
            return "Finite state machine with %s states" % len(self._states_)

    default_format_letter = latex
    format_letter = default_format_letter

    def format_letter_negative(self, letter):
        r"""
        Format negative numbers as overlined numbers, everything
        else by standard LaTeX formatting.

        INPUT:

        - ``letter`` -- anything

        OUTPUT:

        Overlined absolute value if letter is a negative integer,
        :func:`latex(letter) <sage.misc.latex.latex>` otherwise.

        EXAMPLES::

            sage: A = Automaton([(0, 0, -1)])
            sage: list(map(A.format_letter_negative, [-1, 0, 1, 'a', None]))
            ['\\overline{1}', 0, 1, \text{\texttt{a}}, \mathrm{None}]
            sage: A.latex_options(format_letter=A.format_letter_negative)
            sage: print(latex(A))
            \begin{tikzpicture}[auto, initial text=, >=latex]
            \node[state] (v0) at (3.000000, 0.000000) {$0$};
            \path[->] (v0) edge[loop above] node {$\overline{1}$} ();
            \end{tikzpicture}
        """
        if letter in ZZ and letter < 0:
            return r'\overline{%d}' % -letter
        else:
            return latex(letter)

    def format_transition_label_reversed(self, word):
        r"""
        Format words in transition labels in reversed order.

        INPUT:

        - ``word`` -- list of letters

        OUTPUT:

        String representation of ``word`` suitable to be typeset in
        mathematical mode, letters are written in reversed order.

        This is the reversed version of
        :meth:`.default_format_transition_label`.

        In digit expansions, digits are frequently processed from the
        least significant to the most significant position, but it is
        customary to write the least significant digit at the
        right-most position. Therefore, the labels have to be
        reversed.

        EXAMPLES::

            sage: T = Transducer([(0, 0, 0, [1, 2, 3])])
            sage: T.format_transition_label_reversed([1, 2, 3])
            '3 2 1'
            sage: T.latex_options(format_transition_label=T.format_transition_label_reversed)
            sage: print(latex(T))
            \begin{tikzpicture}[auto, initial text=, >=latex]
            \node[state] (v0) at (3.000000, 0.000000) {$0$};
            \path[->] (v0) edge[loop above] node {$0\mid 3 2 1$} ();
            \end{tikzpicture}

        TESTS:

        Check that :issue:`16357` is fixed::

            sage: T = Transducer()
            sage: T.format_transition_label_reversed([])
            '\\varepsilon'
        """
        return self.default_format_transition_label(reversed(word))

    def default_format_transition_label(self, word):
        r"""
        Default formatting of words in transition labels for LaTeX output.

        INPUT:

        - ``word`` -- list of letters

        OUTPUT:

        String representation of ``word`` suitable to be typeset in
        mathematical mode.

        -   For a non-empty word: Concatenation of the letters, piped through
            ``self.format_letter`` and separated by blanks.
        -   For an empty word:
            ``sage.combinat.finite_state_machine.EmptyWordLaTeX``.

        There is also a variant :meth:`.format_transition_label_reversed`
        writing the words in reversed order.

        EXAMPLES:

        #.  Example of a non-empty word::

                sage: T = Transducer()
                sage: print(T.default_format_transition_label(
                ....:    ['a', 'alpha', 'a_1', '0', 0, (0, 1)]))
                \text{\texttt{a}} \text{\texttt{alpha}}
                \text{\texttt{a{\char`\_}1}} 0 0 \left(0, 1\right)

        #.  In the example above, ``'a'`` and ``'alpha'`` should perhaps
            be symbols::

                sage: var('a alpha a_1')                                                # needs sage.symbolic
                (a, alpha, a_1)
                sage: print(T.default_format_transition_label([a, alpha, a_1]))         # needs sage.symbolic
                a \alpha a_{1}

        #.  Example of an empty word::

                sage: print(T.default_format_transition_label([]))
                \varepsilon

            We can change this by setting
            ``sage.combinat.finite_state_machine.EmptyWordLaTeX``::

                sage: sage.combinat.finite_state_machine.EmptyWordLaTeX = ''
                sage: T.default_format_transition_label([])
                ''

            Finally, we restore the default value::

                sage: sage.combinat.finite_state_machine.EmptyWordLaTeX = r'\varepsilon'

        #.  This method is the default value for
            ``FiniteStateMachine.format_transition_label``. That can be changed to be
            any other function::

                sage: A = Automaton([(0, 1, 0)])
                sage: def custom_format_transition_label(word):
                ....:     return "t"
                sage: A.latex_options(format_transition_label=custom_format_transition_label)
                sage: print(latex(A))
                \begin{tikzpicture}[auto, initial text=, >=latex]
                \node[state] (v0) at (3.000000, 0.000000) {$0$};
                \node[state] (v1) at (-3.000000, 0.000000) {$1$};
                \path[->] (v0) edge node[rotate=360.00, anchor=south] {$t$} (v1);
                \end{tikzpicture}

        TESTS:

        Check that :issue:`16357` is fixed::

            sage: T = Transducer()
            sage: T.default_format_transition_label([])
            '\\varepsilon'
            sage: T.default_format_transition_label(iter([]))
            '\\varepsilon'
        """
        result = " ".join(self.format_letter(u) for u in word)
        if result:
            return result
        else:
            return EmptyWordLaTeX

    format_transition_label = default_format_transition_label

    def latex_options(self,
                      coordinates=None,
                      format_state_label=None,
                      format_letter=None,
                      format_transition_label=None,
                      loop_where=None,
                      initial_where=None,
                      accepting_style=None,
                      accepting_distance=None,
                      accepting_where=None,
                      accepting_show_empty=None):
        r"""
        Set options for LaTeX output via
        :func:`~sage.misc.latex.latex` and therefore
        :func:`~sage.misc.latex.view`.

        INPUT:

        - ``coordinates`` -- dictionary or a function mapping labels
          of states to pairs interpreted as coordinates. If no
          coordinates are given, states a placed equidistantly on a
          circle of radius `3`. See also :meth:`.set_coordinates`.

        - ``format_state_label`` -- a function mapping labels of
          states to a string suitable for typesetting in LaTeX's
          mathematics mode. If not given, :func:`~sage.misc.latex.latex`
          is used.

        - ``format_letter`` -- a function mapping letters of the input
          and output alphabets to a string suitable for typesetting in
          LaTeX's mathematics mode. If not given,
          :meth:`.default_format_transition_label` uses
          :func:`~sage.misc.latex.latex`.

        - ``format_transition_label`` -- a function mapping words over
          the input and output alphabets to a string suitable for
          typesetting in LaTeX's mathematics mode. If not given,
          :meth:`.default_format_transition_label` is used.

        - ``loop_where`` -- dictionary or a function mapping labels of
          initial states to one of ``'above'``, ``'left'``, ``'below'``,
          ``'right'``. If not given, ``'above'`` is used.

        - ``initial_where`` -- dictionary or a function mapping
          labels of initial states to one of ``'above'``, ``'left'``,
          ``'below'``, ``'right'``. If not given, TikZ' default
          (currently ``'left'``) is used.

        - ``accepting_style`` -- one of ``'accepting by double'`` and
          ``'accepting by arrow'``. If not given, ``'accepting by
          double'`` is used unless there are non-empty final output
          words.

        - ``accepting_distance`` -- string giving a LaTeX length
          used for the length of the arrow leading from a final state.
          If not given, TikZ' default (currently ``'3ex'``) is used
          unless there are non-empty final output words, in which case
          ``'7ex'`` is used.

        - ``accepting_where`` -- dictionary or a function mapping
          labels of final states to one of ``'above'``, ``'left'``,
          ``'below'``, ``'right'``. If not given, TikZ' default
          (currently ``'right'``) is used. If the final state has a
          final output word, it is also possible to give an angle
          in degrees.

        - ``accepting_show_empty`` -- if ``True`` the arrow of an
          empty final output word is labeled as well. Note that this
          implicitly implies ``accepting_style='accepting by
          arrow'``. If not given, the default ``False`` is used.

        OUTPUT: nothing

        As TikZ (cf. the :wikipedia:`PGF/TikZ`) is used to typeset
        the graphics, the syntax is oriented on TikZ' syntax.

        This is a convenience function collecting all options for
        LaTeX output. All of its functionality can also be achieved by
        directly setting the attributes

        - ``coordinates``, ``format_label``, ``loop_where``,
          ``initial_where``, and ``accepting_where`` of
          :class:`FSMState` (here, ``format_label`` is a callable
          without arguments, everything else is a specific value);

        - ``format_label`` of :class:`FSMTransition` (``format_label``
          is a callable without arguments);

        - ``format_state_label``, ``format_letter``,
          ``format_transition_label``, ``accepting_style``,
          ``accepting_distance``, and ``accepting_show_empty``
          of :class:`FiniteStateMachine`.

        This function, however, also (somewhat) checks its input and
        serves to collect documentation on all these options.

        The function can be called several times, only those arguments
        which are not ``None`` are taken into account. By the same
        means, it can be combined with directly setting some
        attributes as outlined above.

        EXAMPLES:

        See also the section on :ref:`finite_state_machine_LaTeX_output`
        in the introductory examples of this module.

        ::

            sage: T = Transducer(initial_states=[4],
            ....:     final_states=[0, 3])
            sage: for j in srange(4):
            ....:     T.add_transition(4, j, 0, [0, j])
            ....:     T.add_transition(j, 4, 0, [0, -j])
            ....:     T.add_transition(j, j, 0, 0)
            Transition from 4 to 0: 0|0,0
            Transition from 0 to 4: 0|0,0
            Transition from 0 to 0: 0|0
            Transition from 4 to 1: 0|0,1
            Transition from 1 to 4: 0|0,-1
            Transition from 1 to 1: 0|0
            Transition from 4 to 2: 0|0,2
            Transition from 2 to 4: 0|0,-2
            Transition from 2 to 2: 0|0
            Transition from 4 to 3: 0|0,3
            Transition from 3 to 4: 0|0,-3
            Transition from 3 to 3: 0|0
            sage: T.add_transition(4, 4, 0, 0)
            Transition from 4 to 4: 0|0
            sage: T.state(3).final_word_out = [0, 0]
            sage: T.latex_options(
            ....:     coordinates={4: (0, 0),
            ....:                  0: (-6, 3),
            ....:                  1: (-2, 3),
            ....:                  2: (2, 3),
            ....:                  3: (6, 3)},
            ....:     format_state_label=lambda x: r'\mathbf{%s}' % x,
            ....:     format_letter=lambda x: r'w_{%s}' % x,
            ....:     format_transition_label=lambda x:
            ....:         r"{\scriptstyle %s}" % T.default_format_transition_label(x),
            ....:     loop_where={4: 'below', 0: 'left', 1: 'above',
            ....:                 2: 'right', 3:'below'},
            ....:     initial_where=lambda x: 'above',
            ....:     accepting_style='accepting by double',
            ....:     accepting_distance='10ex',
            ....:     accepting_where={0: 'left', 3: 45}
            ....:     )
            sage: T.state(4).format_label=lambda: r'\mathcal{I}'
            sage: latex(T)
            \begin{tikzpicture}[auto, initial text=, >=latex]
            \node[state, initial, initial where=above] (v0) at (0.000000, 0.000000) {$\mathcal{I}$};
            \node[state, accepting, accepting where=left] (v1) at (-6.000000, 3.000000) {$\mathbf{0}$};
            \node[state, accepting, accepting where=45] (v2) at (6.000000, 3.000000) {$\mathbf{3}$};
            \path[->] (v2.45.00) edge node[rotate=45.00, anchor=south] {$\$ \mid {\scriptstyle w_{0} w_{0}}$} ++(45.00:10ex);
            \node[state] (v3) at (-2.000000, 3.000000) {$\mathbf{1}$};
            \node[state] (v4) at (2.000000, 3.000000) {$\mathbf{2}$};
            \path[->] (v1) edge[loop left] node[rotate=90, anchor=south] {${\scriptstyle w_{0}}\mid {\scriptstyle w_{0}}$} ();
            \path[->] (v1.-21.57) edge node[rotate=-26.57, anchor=south] {${\scriptstyle w_{0}}\mid {\scriptstyle w_{0} w_{0}}$} (v0.148.43);
            \path[->] (v3) edge[loop above] node {${\scriptstyle w_{0}}\mid {\scriptstyle w_{0}}$} ();
            \path[->] (v3.-51.31) edge node[rotate=-56.31, anchor=south] {${\scriptstyle w_{0}}\mid {\scriptstyle w_{0} w_{-1}}$} (v0.118.69);
            \path[->] (v4) edge[loop right] node[rotate=90, anchor=north] {${\scriptstyle w_{0}}\mid {\scriptstyle w_{0}}$} ();
            \path[->] (v4.-118.69) edge node[rotate=56.31, anchor=north] {${\scriptstyle w_{0}}\mid {\scriptstyle w_{0} w_{-2}}$} (v0.51.31);
            \path[->] (v2) edge[loop below] node {${\scriptstyle w_{0}}\mid {\scriptstyle w_{0}}$} ();
            \path[->] (v2.-148.43) edge node[rotate=26.57, anchor=north] {${\scriptstyle w_{0}}\mid {\scriptstyle w_{0} w_{-3}}$} (v0.21.57);
            \path[->] (v0.158.43) edge node[rotate=333.43, anchor=north] {${\scriptstyle w_{0}}\mid {\scriptstyle w_{0} w_{0}}$} (v1.328.43);
            \path[->] (v0.128.69) edge node[rotate=303.69, anchor=north] {${\scriptstyle w_{0}}\mid {\scriptstyle w_{0} w_{1}}$} (v3.298.69);
            \path[->] (v0.61.31) edge node[rotate=56.31, anchor=south] {${\scriptstyle w_{0}}\mid {\scriptstyle w_{0} w_{2}}$} (v4.231.31);
            \path[->] (v0.31.57) edge node[rotate=26.57, anchor=south] {${\scriptstyle w_{0}}\mid {\scriptstyle w_{0} w_{3}}$} (v2.201.57);
            \path[->] (v0) edge[loop below] node {${\scriptstyle w_{0}}\mid {\scriptstyle w_{0}}$} ();
            \end{tikzpicture}
            sage: view(T) # not tested

        To actually see this, use the live documentation in the Sage notebook
        and execute the cells.

        By changing some of the options, we get the following output::

            sage: T.latex_options(
            ....:     format_transition_label=T.default_format_transition_label,
            ....:     accepting_style='accepting by arrow',
            ....:     accepting_show_empty=True
            ....:     )
            sage: latex(T)
            \begin{tikzpicture}[auto, initial text=, >=latex, accepting text=, accepting/.style=accepting by arrow, accepting distance=10ex]
            \node[state, initial, initial where=above] (v0) at (0.000000, 0.000000) {$\mathcal{I}$};
            \node[state] (v1) at (-6.000000, 3.000000) {$\mathbf{0}$};
            \path[->] (v1.180.00) edge node[rotate=360.00, anchor=south] {$\$ \mid \varepsilon$} ++(180.00:10ex);
            \node[state] (v2) at (6.000000, 3.000000) {$\mathbf{3}$};
            \path[->] (v2.45.00) edge node[rotate=45.00, anchor=south] {$\$ \mid w_{0} w_{0}$} ++(45.00:10ex);
            \node[state] (v3) at (-2.000000, 3.000000) {$\mathbf{1}$};
            \node[state] (v4) at (2.000000, 3.000000) {$\mathbf{2}$};
            \path[->] (v1) edge[loop left] node[rotate=90, anchor=south] {$w_{0}\mid w_{0}$} ();
            \path[->] (v1.-21.57) edge node[rotate=-26.57, anchor=south] {$w_{0}\mid w_{0} w_{0}$} (v0.148.43);
            \path[->] (v3) edge[loop above] node {$w_{0}\mid w_{0}$} ();
            \path[->] (v3.-51.31) edge node[rotate=-56.31, anchor=south] {$w_{0}\mid w_{0} w_{-1}$} (v0.118.69);
            \path[->] (v4) edge[loop right] node[rotate=90, anchor=north] {$w_{0}\mid w_{0}$} ();
            \path[->] (v4.-118.69) edge node[rotate=56.31, anchor=north] {$w_{0}\mid w_{0} w_{-2}$} (v0.51.31);
            \path[->] (v2) edge[loop below] node {$w_{0}\mid w_{0}$} ();
            \path[->] (v2.-148.43) edge node[rotate=26.57, anchor=north] {$w_{0}\mid w_{0} w_{-3}$} (v0.21.57);
            \path[->] (v0.158.43) edge node[rotate=333.43, anchor=north] {$w_{0}\mid w_{0} w_{0}$} (v1.328.43);
            \path[->] (v0.128.69) edge node[rotate=303.69, anchor=north] {$w_{0}\mid w_{0} w_{1}$} (v3.298.69);
            \path[->] (v0.61.31) edge node[rotate=56.31, anchor=south] {$w_{0}\mid w_{0} w_{2}$} (v4.231.31);
            \path[->] (v0.31.57) edge node[rotate=26.57, anchor=south] {$w_{0}\mid w_{0} w_{3}$} (v2.201.57);
            \path[->] (v0) edge[loop below] node {$w_{0}\mid w_{0}$} ();
            \end{tikzpicture}
            sage: view(T) # not tested

        TESTS::

            sage: T.latex_options(format_state_label='Nothing')
            Traceback (most recent call last):
            ...
            TypeError: format_state_label must be callable.
            sage: T.latex_options(format_letter='')
            Traceback (most recent call last):
            ...
            TypeError: format_letter must be callable.
            sage: T.latex_options(format_transition_label='')
            Traceback (most recent call last):
            ...
            TypeError: format_transition_label must be callable.
            sage: T.latex_options(loop_where=37)
            Traceback (most recent call last):
            ...
            TypeError: loop_where must be a callable or a
            dictionary.
            sage: T.latex_options(loop_where=lambda x: 'top')
            Traceback (most recent call last):
            ...
            ValueError: loop_where for 4 must be in ['above',
            'below', 'left', 'right'].
            sage: T.latex_options(initial_where=90)
            Traceback (most recent call last):
            ...
            TypeError: initial_where must be a callable or a
            dictionary.
            sage: T.latex_options(initial_where=lambda x: 'top')
            Traceback (most recent call last):
            ...
            ValueError: initial_where for 4 must be in ['above',
            'below', 'left', 'right'].
            sage: T.latex_options(accepting_style='fancy')
            Traceback (most recent call last):
            ...
            ValueError: accepting_style must be in ['accepting by
            arrow', 'accepting by double'].
            sage: T.latex_options(accepting_where=90)
            Traceback (most recent call last):
            ...
            TypeError: accepting_where must be a callable or a
            dictionary.
            sage: T.latex_options(accepting_where=lambda x: 'top')
            Traceback (most recent call last):
            ...
            ValueError: accepting_where for 0 must be in ['above',
            'below', 'left', 'right'].
            sage: T.latex_options(accepting_where={0: 'above', 3: 'top'})
            Traceback (most recent call last):
            ...
            ValueError: accepting_where for 3 must be a real number or
            be in ['above', 'below', 'left', 'right'].
        """
        if coordinates is not None:
            self.set_coordinates(coordinates)

        if format_state_label is not None:
            if not callable(format_state_label):
                raise TypeError('format_state_label must be callable.')
            self.format_state_label = format_state_label

        if format_letter is not None:
            if not callable(format_letter):
                raise TypeError('format_letter must be callable.')
            self.format_letter = format_letter

        if format_transition_label is not None:
            if not callable(format_transition_label):
                raise TypeError('format_transition_label must be callable.')
            self.format_transition_label = format_transition_label

        if loop_where is not None:
            permissible = list(tikz_automata_where)
            for state in self.states():
                if callable(loop_where):
                    where = loop_where(state.label())
                else:
                    try:
                        where = loop_where[state.label()]
                    except TypeError:
                        raise TypeError("loop_where must be a "
                                        "callable or a dictionary.")
                    except KeyError:
                        continue
                if where in permissible:
                    state.loop_where = where
                else:
                    raise ValueError('loop_where for %s must be in %s.' %
                                     (state.label(), sorted(permissible)))

        if initial_where is not None:
            permissible = list(tikz_automata_where)
            for state in self.iter_initial_states():
                if callable(initial_where):
                    where = initial_where(state.label())
                else:
                    try:
                        where = initial_where[state.label()]
                    except TypeError:
                        raise TypeError("initial_where must be a "
                                        "callable or a dictionary.")
                    except KeyError:
                        continue
                if where in permissible:
                    state.initial_where = where
                else:
                    raise ValueError('initial_where for %s must be in %s.' %
                                     (state.label(), sorted(permissible)))

        if accepting_style is not None:
            permissible = ['accepting by double',
                           'accepting by arrow']
            if accepting_style in permissible:
                self.accepting_style = accepting_style
            else:
                raise ValueError('accepting_style must be in %s.' %
                    sorted(permissible))

        if accepting_distance is not None:
            self.accepting_distance = accepting_distance

        if accepting_where is not None:
            permissible = list(tikz_automata_where)
            for state in self.iter_final_states():
                if callable(accepting_where):
                    where = accepting_where(state.label())
                else:
                    try:
                        where = accepting_where[state.label()]
                    except TypeError:
                        raise TypeError("accepting_where must be a "
                                        "callable or a dictionary.")
                    except KeyError:
                        continue
                if where in permissible:
                    state.accepting_where = where
                elif hasattr(state, 'final_word_out') \
                        and state.final_word_out:
                    if where in RR:
                        state.accepting_where = where
                    else:
                        raise ValueError('accepting_where for %s must '
                                         'be a real number or be in %s.' %
                                         (state.label(), sorted(permissible)))

                else:
                    raise ValueError('accepting_where for %s must be in %s.' %
                                     (state.label(), sorted(permissible)))

        if accepting_show_empty is not None:
            self.accepting_show_empty = accepting_show_empty

    def _latex_(self):
        r"""
        Return a LaTeX code for the graph of the finite state machine.

        OUTPUT: string

        EXAMPLES::

            sage: F = FiniteStateMachine([('A', 'B', 1, 2)],
            ....:                        initial_states=['A'],
            ....:                        final_states=['B'])
            sage: F.state('A').initial_where='below'
            sage: print(latex(F))  # indirect doctest
            \begin{tikzpicture}[auto, initial text=, >=latex]
            \node[state, initial, initial where=below] (v0) at (3.000000, 0.000000) {$\text{\texttt{A}}$};
            \node[state, accepting] (v1) at (-3.000000, 0.000000) {$\text{\texttt{B}}$};
            \path[->] (v0) edge node[rotate=360.00, anchor=south] {$ $} (v1);
            \end{tikzpicture}

        TESTS:

            Check that :issue:`16943` is fixed::

                sage: latex(Transducer(
                ....:     [(0, 1), (1, 1), (2, 2), (3, 3), (4, 4)]))
                \begin{tikzpicture}[auto, initial text=, >=latex]
                \node[state] (v0) at (3.000000, 0.000000) {$0$};
                \node[state] (v1) at (0.927051, 2.853170) {$1$};
                \node[state] (v2) at (-2.427051, 1.763356) {$2$};
                \node[state] (v3) at (-2.427051, -1.763356) {$3$};
                \node[state] (v4) at (0.927051, -2.853170) {$4$};
                \path[->] (v0) edge node[rotate=306.00, anchor=south] {$\varepsilon\mid \varepsilon$} (v1);
                \path[->] (v1) edge[loop above] node {$\varepsilon\mid \varepsilon$} ();
                \path[->] (v2) edge[loop above] node {$\varepsilon\mid \varepsilon$} ();
                \path[->] (v3) edge[loop above] node {$\varepsilon\mid \varepsilon$} ();
                \path[->] (v4) edge[loop above] node {$\varepsilon\mid \varepsilon$} ();
                \end{tikzpicture}
        """
        from math import sin, cos, pi, atan2

        def label_rotation(angle, both_directions):
            """
            Given an angle of a transition, compute the TikZ string to
            rotate the label.
            """
            angle_label = angle
            anchor_label = "south"
            if angle > 90 or angle <= -90:
                angle_label = angle + 180
                if both_directions:
                    # if transitions in both directions, the transition to the
                    # left has its label below the transition, otherwise above
                    anchor_label = "north"
            if hasattr(angle_label, 'n'):
                # we may need to convert symbolic expressions to floats,
                # but int does not have .n()
                angle_label = angle_label.n()
            return "rotate=%.2f, anchor=%s" % (angle_label, anchor_label)

        setup_latex_preamble()

        options = ["auto", "initial text=", ">=latex"]

        nonempty_final_word_out = False
        for state in self.iter_final_states():
            if state.final_word_out:
                nonempty_final_word_out = True
                break

        if hasattr(self, "accepting_style"):
            accepting_style = self.accepting_style
        elif nonempty_final_word_out:
            accepting_style = "accepting by arrow"
        else:
            accepting_style = "accepting by double"

        if accepting_style == "accepting by arrow":
            options.append("accepting text=")
            options.append("accepting/.style=%s" % accepting_style)

        if hasattr(self, "accepting_distance"):
            accepting_distance = self.accepting_distance
        elif nonempty_final_word_out:
            accepting_distance = "7ex"
        else:
            accepting_distance = None
        if accepting_style == "accepting by arrow" and accepting_distance:
            options.append("accepting distance=%s"
                           % accepting_distance)

        if hasattr(self, "accepting_show_empty"):
            accepting_show_empty = self.accepting_show_empty
        else:
            accepting_show_empty = False

        result = "\\begin{tikzpicture}[%s]\n" % ", ".join(options)
        j = 0
        for vertex in self.iter_states():
            if not hasattr(vertex, "coordinates"):
                vertex.coordinates = (3*cos(2*pi*j/len(self.states())),
                                      3*sin(2*pi*j/len(self.states())))
            options = ""
            if vertex.is_final:
                if not (vertex.final_word_out
                        and accepting_style == "accepting by arrow") \
                        and not accepting_show_empty:
                    # otherwise, we draw a custom made accepting path
                    # with label below
                    options += ", accepting"
                    if hasattr(vertex, "accepting_where"):
                        options += ", accepting where=%s" % (
                            vertex.accepting_where,)
            if vertex.is_initial:
                options += ", initial"
            if hasattr(vertex, "initial_where"):
                options += ", initial where=%s" % vertex.initial_where
            if hasattr(vertex, "format_label"):
                label = vertex.format_label()
            elif hasattr(self, "format_state_label"):
                label = self.format_state_label(vertex)
            else:
                label = latex(vertex.label())
            result += "\\node[state%s] (v%d) at (%f, %f) {$%s$};\n" % (
                options, j, vertex.coordinates[0],
                vertex.coordinates[1], label)
            vertex._number_ = j
            if vertex.is_final and (vertex.final_word_out or accepting_show_empty):
                angle = 0
                if hasattr(vertex, "accepting_where"):
                    angle = tikz_automata_where.get(vertex.accepting_where,
                                                    vertex.accepting_where)
                result += "\\path[->] (v%d.%.2f) edge node[%s] {$%s \\mid %s$} ++(%.2f:%s);\n" % (
                    j, angle,
                    label_rotation(angle, False),
                    EndOfWordLaTeX,
                    self.format_transition_label(vertex.final_word_out),
                    angle, accepting_distance)

            j += 1

        def key_function(s):
            return (s.from_state, s.to_state)
        # We use an OrderedDict instead of a dict in order to have a
        # defined ordering of the transitions in the output. See
        # https://github.com/sagemath/sage/issues/16580#comment:3 . As the
        # transitions have to be sorted anyway, the performance
        # penalty should be bearable; nevertheless, this is only
        # required for doctests.
        adjacent = OrderedDict(
            (pair, list(transitions))
            for pair, transitions in
            itertools.groupby(
                sorted(self.iter_transitions(),
                       key=key_function),
                key=key_function
                ))

        for ((source, target), transitions) in adjacent.items():
            if transitions:
                labels = []
                for transition in transitions:
                    if hasattr(transition, "format_label"):
                        labels.append(transition.format_label())
                    else:
                        labels.append(self._latex_transition_label_(
                            transition, self.format_transition_label))
                label = ", ".join(labels)
                if source != target:
                    angle = atan2(
                        target.coordinates[1] - source.coordinates[1],
                        target.coordinates[0] - source.coordinates[0]) * 180/pi
                    both_directions = (target, source) in adjacent
                    if both_directions:
                        angle_source = ".%.2f" % (angle + 5)
                        angle_target = ".%.2f" % (angle + 175)
                    else:
                        angle_source = ""
                        angle_target = ""
                    result += "\\path[->] (v%d%s) edge node[%s] {$%s$} (v%d%s);\n" % (
                        source._number_, angle_source,
                        label_rotation(angle, both_directions),
                        label,
                        target._number_, angle_target)
                else:
                    loop_where = "above"
                    if hasattr(source, "loop_where"):
                        loop_where = source.loop_where
                    rotation = {'left': '[rotate=90, anchor=south]',
                                'right': '[rotate=90, anchor=north]'}
                    result += "\\path[->] (v%d) edge[loop %s] node%s {$%s$} ();\n" % (
                        source._number_,
                        loop_where, rotation.get(loop_where, ''),
                        label)

        result += "\\end{tikzpicture}"
        return result

    def _latex_transition_label_(self, transition,
                                 format_function=None):
        r"""
        Return the proper transition label.

        INPUT:

        - ``transition`` -- a transition

        - ``format_function`` -- a function formatting the labels

        OUTPUT: string

        TESTS::

            sage: F = FiniteStateMachine([('A', 'B', 0, 1)])
            sage: t = F.transitions()[0]
            sage: F._latex_transition_label_(t)
            ' '
        """
        return ' '

    def set_coordinates(self, coordinates, default=True):
        """
        Set coordinates of the states for the LaTeX representation by
        a dictionary or a function mapping labels to coordinates.

        INPUT:

        - ``coordinates`` -- dictionary or a function mapping labels
          of states to pairs interpreted as coordinates

        - ``default`` -- if ``True``, then states not given by
          ``coordinates`` get a default position on a circle of
          radius 3.

        OUTPUT: nothing

        EXAMPLES::

            sage: F = Automaton([[0, 1, 1], [1, 2, 2], [2, 0, 0]])
            sage: F.set_coordinates({0: (0, 0), 1: (2, 0), 2: (1, 1)})
            sage: F.state(0).coordinates
            (0, 0)

        We can also use a function to determine the coordinates::

            sage: F = Automaton([[0, 1, 1], [1, 2, 2], [2, 0, 0]])
            sage: F.set_coordinates(lambda l: (l, 3/(l+1)))
            sage: F.state(2).coordinates
            (2, 1)
        """
        from math import sin, cos, pi

        states_without_coordinates = []
        for state in self.iter_states():
            try:
                state.coordinates = coordinates[state.label()]
                continue
            except (KeyError, TypeError):
                pass

            try:
                state.coordinates = coordinates(state.label())
                continue
            except TypeError:
                pass

            states_without_coordinates.append(state)

        if default:
            n = len(states_without_coordinates)
            for j, state in enumerate(states_without_coordinates):
                state.coordinates = (3*cos(2*pi*j/n),
                                     3*sin(2*pi*j/n))

    # ************************************************************************
    # other
    # ************************************************************************

    def _matrix_(self, R=None):
        """
        Return the adjacency matrix of the finite state machine.

        See :meth:`.adjacency_matrix` for more information.

        EXAMPLES::

            sage: B = FiniteStateMachine({0: {0: (0, 0), 'a': (1, 0)},
            ....:                         'a': {2: (0, 0), 3: (1, 0)},
            ....:                         2:{0:(1, 1), 4:(0, 0)},
            ....:                         3:{'a':(0, 1), 2:(1, 1)},
            ....:                         4:{4:(1, 1), 3:(0, 1)}},
            ....:                        initial_states=[0])
            sage: B._matrix_()                                                          # needs sage.symbolic
            [1 1 0 0 0]
            [0 0 1 1 0]
            [x 0 0 0 1]
            [0 x x 0 0]
            [0 0 0 x x]
        """
        return self.adjacency_matrix()

    def adjacency_matrix(self, input=None,
                         entry=None):
        """
        Return the adjacency matrix of the underlying graph.

        INPUT:

        - ``input`` -- only transitions with input label ``input`` are
          respected

        - ``entry`` -- the function ``entry`` takes a transition and the
          return value is written in the matrix as the entry
          ``(transition.from_state, transition.to_state)``. The default
          value (``None``) of entry takes the variable ``x`` to the
          power of the sum of the output word of the transition.

        OUTPUT: a matrix

        If any label of a state is not an integer, the finite state
        machine is relabeled at the beginning.  If there are more than
        one transitions between two states, then the different return
        values of ``entry`` are added up.

        EXAMPLES::

            sage: B = FiniteStateMachine({0:{0:(0, 0), 'a':(1, 0)},
            ....:                         'a':{2:(0, 0), 3:(1, 0)},
            ....:                         2:{0:(1, 1), 4:(0, 0)},
            ....:                         3:{'a':(0, 1), 2:(1, 1)},
            ....:                         4:{4:(1, 1), 3:(0, 1)}},
            ....:                        initial_states=[0])
            sage: B.adjacency_matrix()                                                  # needs sage.symbolic
            [1 1 0 0 0]
            [0 0 1 1 0]
            [x 0 0 0 1]
            [0 x x 0 0]
            [0 0 0 x x]

        This is equivalent to::

            sage: matrix(B)                                                             # needs sage.symbolic
            [1 1 0 0 0]
            [0 0 1 1 0]
            [x 0 0 0 1]
            [0 x x 0 0]
            [0 0 0 x x]

        It is also possible to use other entries in the adjacency matrix::

            sage: B.adjacency_matrix(entry=(lambda transition: 1))
            [1 1 0 0 0]
            [0 0 1 1 0]
            [1 0 0 0 1]
            [0 1 1 0 0]
            [0 0 0 1 1]
            sage: var('t')                                                              # needs sage.symbolic
            t
            sage: B.adjacency_matrix(1, entry=(lambda transition:                       # needs sage.symbolic
            ....:     exp(I*transition.word_out[0]*t)))
            [      0       1       0       0       0]
            [      0       0       0       1       0]
            [e^(I*t)       0       0       0       0]
            [      0       0 e^(I*t)       0       0]
            [      0       0       0       0 e^(I*t)]
            sage: a = Automaton([(0, 1, 0),
            ....:                (1, 2, 0),
            ....:                (2, 0, 1),
            ....:                (2, 1, 0)],
            ....:               initial_states=[0],
            ....:               final_states=[0])
            sage: a.adjacency_matrix()                                                  # needs sage.symbolic
            [0 1 0]
            [0 0 1]
            [1 1 0]
        """

        if entry is None:
            from sage.symbolic.ring import SR
            x = SR.var('x')

            def default_function(transition):
                return x**sum(transition.word_out)

            entry = default_function

        relabeledFSM = self
        l = len(relabeledFSM.states())
        for state in self.iter_states():
            if state.label() not in ZZ or state.label() >= l \
                     or state.label() < 0:
                relabeledFSM = self.relabeled()
                break
        dictionary = {}
        for transition in relabeledFSM.iter_transitions():
            if input is None or transition.word_in == [input]:
                if (transition.from_state.label(),
                    transition.to_state.label()) in dictionary:
                    dictionary[(transition.from_state.label(),
                                transition.to_state.label())] \
                                += entry(transition)
                else:
                    dictionary[(transition.from_state.label(),
                                transition.to_state.label())] \
                                = entry(transition)
        return matrix(
            len(relabeledFSM.states()), dictionary)

    def determine_input_alphabet(self, reset=True):
        """
        Determine the input alphabet according to the transitions
        of this finite state machine.

        INPUT:

        - ``reset`` -- boolean (default: ``True``); if ``True``, then
          the existing input alphabet is erased, otherwise new letters are
          appended to the existing alphabet.

        OUTPUT: nothing

        After this operation the input alphabet of this finite state machine
        is a list of letters.

        .. TODO::

            At the moment, the letters of the alphabet need to be hashable.

        EXAMPLES::

            sage: T = Transducer([(1, 1, 1, 0), (1, 2, 2, 1),
            ....:                 (2, 2, 1, 1), (2, 2, 0, 0)],
            ....:                final_states=[1],
            ....:                determine_alphabets=False)
            sage: (T.input_alphabet, T.output_alphabet)
            (None, None)
            sage: T.determine_input_alphabet()
            sage: (T.input_alphabet, T.output_alphabet)
            ([0, 1, 2], None)

        .. SEEALSO::

           :meth:`determine_output_alphabet`,
           :meth:`determine_alphabets`.
        """
        if reset:
            ain = set()
        else:
            ain = set(self.input_alphabet)

        for t in self.iter_transitions():
            for letter in t.word_in:
                ain.add(letter)
        self.input_alphabet = list(ain)

    def determine_output_alphabet(self, reset=True):
        """
        Determine the output alphabet according to the transitions
        of this finite state machine.

        INPUT:

        - ``reset`` -- boolean (default: ``True``); if ``True``, then
          the existing output alphabet is erased, otherwise new letters are
          appended to the existing alphabet.

        OUTPUT: nothing

        After this operation the output alphabet of this finite state machine
        is a list of letters.

        .. TODO::

            At the moment, the letters of the alphabet need to be hashable.

        EXAMPLES::

            sage: T = Transducer([(1, 1, 1, 0), (1, 2, 2, 1),
            ....:                 (2, 2, 1, 1), (2, 2, 0, 0)],
            ....:                final_states=[1],
            ....:                determine_alphabets=False)
            sage: T.state(1).final_word_out = [1, 4]
            sage: (T.input_alphabet, T.output_alphabet)
            (None, None)
            sage: T.determine_output_alphabet()
            sage: (T.input_alphabet, T.output_alphabet)
            (None, [0, 1, 4])

        .. SEEALSO::

           :meth:`determine_input_alphabet`,
           :meth:`determine_alphabets`.
        """
        if reset:
            aout = set()
        else:
            aout = set(self.output_alphabet)

        for t in self.iter_transitions():
            for letter in t.word_out:
                aout.add(letter)
        for s in self.iter_final_states():
            for letter in s.final_word_out:
                aout.add(letter)
        self.output_alphabet = list(aout)

    def determine_alphabets(self, reset=True):
        """
        Determine the input and output alphabet according to the
        transitions in this finite state machine.

        INPUT:

        - ``reset`` -- if reset is ``True``, then the existing input
          and output alphabets are erased, otherwise new letters are
          appended to the existing alphabets.

        OUTPUT: nothing

        After this operation the input alphabet and the output
        alphabet of this finite state machine are a list of letters.

        .. TODO::

            At the moment, the letters of the alphabets need to be hashable.

        EXAMPLES::

            sage: T = Transducer([(1, 1, 1, 0), (1, 2, 2, 1),
            ....:                 (2, 2, 1, 1), (2, 2, 0, 0)],
            ....:                final_states=[1],
            ....:                determine_alphabets=False)
            sage: T.state(1).final_word_out = [1, 4]
            sage: (T.input_alphabet, T.output_alphabet)
            (None, None)
            sage: T.determine_alphabets()
            sage: (T.input_alphabet, T.output_alphabet)
            ([0, 1, 2], [0, 1, 4])

        .. SEEALSO::

           :meth:`determine_input_alphabet`,
           :meth:`determine_output_alphabet`.
        """
        self.determine_input_alphabet(reset)
        self.determine_output_alphabet(reset)

    # ************************************************************************
    # get states and transitions
    # ************************************************************************

    def states(self):
        """
        Return the states of the finite state machine.

        OUTPUT: the states of the finite state machine as list

        EXAMPLES::

            sage: FSM = Automaton([('1', '2', 1), ('2', '2', 0)])
            sage: FSM.states()
            ['1', '2']
        """
        return copy(self._states_)

    def iter_states(self):
        """
        Return an iterator of the states.

        OUTPUT: an iterator of the states of the finite state machine

        EXAMPLES::

            sage: FSM = Automaton([('1', '2', 1), ('2', '2', 0)])
            sage: [s.label() for s in FSM.iter_states()]
            ['1', '2']
        """
        return iter(self._states_)

    def transitions(self, from_state=None):
        """
        Return a list of all transitions.

        INPUT:

        - ``from_state`` -- (default: ``None``) if ``from_state`` is
          given, then a list of transitions starting there is given

        OUTPUT: list of all transitions

        EXAMPLES::

            sage: FSM = Automaton([('1', '2', 1), ('2', '2', 0)])
            sage: FSM.transitions()
            [Transition from '1' to '2': 1|-,
             Transition from '2' to '2': 0|-]
        """
        return list(self.iter_transitions(from_state))

    def iter_transitions(self, from_state=None):
        """
        Return an iterator of all transitions.

        INPUT:

        - ``from_state`` -- (default: ``None``) if ``from_state`` is
          given, then a list of transitions starting there is given

        OUTPUT: an iterator of all transitions

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
        Return an iterator over all transitions.

        OUTPUT: an iterator over all transitions

        EXAMPLES::

            sage: FSM = Automaton([('1', '2', 1), ('2', '2', 0)])
            sage: [(t.from_state.label(), t.to_state.label())
            ....:     for t in FSM._iter_transitions_all_()]
            [('1', '2'), ('2', '2')]
        """
        for state in self.iter_states():
            yield from state.transitions

    def initial_states(self):
        """
        Return a list of all initial states.

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
        Return an iterator of the initial states.

        OUTPUT: an iterator over all initial states

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A', is_initial=True)
            sage: B = FSMState('B')
            sage: F = FiniteStateMachine([(A, B, 1, 0)])
            sage: [s.label() for s in F.iter_initial_states()]
            ['A']
        """
        return (s for s in self.iter_states() if s.is_initial)

    def final_states(self):
        """
        Return a list of all final states.

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
        Return an iterator of the final states.

        OUTPUT: an iterator over all initial states

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A', is_final=True)
            sage: B = FSMState('B', is_initial=True)
            sage: C = FSMState('C', is_final=True)
            sage: F = FiniteStateMachine([(A, B), (A, C)])
            sage: [s.label() for s in F.iter_final_states()]
            ['A', 'C']
        """
        return (s for s in self.iter_states() if s.is_final)

    def state(self, state):
        """
        Return the state of the finite state machine.

        INPUT:

        - ``state`` -- if ``state`` is not an instance of
          :class:`FSMState`, then it is assumed that it is the label
          of a state.

        OUTPUT: the state of the finite state machine corresponding to ``state``

        If no state is found, then a :exc:`LookupError` is thrown.

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
        switch = isinstance(state, FSMState)

        try:
            return self._states_dict_[what(state, switch)]
        except AttributeError:
            for s in self.iter_states():
                if what(s, not switch) == state:
                    return s
        except KeyError:
            pass
        raise LookupError("No state with label %s found." % (what(state, switch),))

    def transition(self, transition):
        """
        Return the transition of the finite state machine.

        INPUT:

        - ``transition`` -- if ``transition`` is not an instance of
          :class:`FSMTransition`, then it is assumed that it is a
          tuple ``(from_state, to_state, word_in, word_out)``

        OUTPUT: the transition of the finite state machine corresponding
        to ``transition``

        If no transition is found, then a :exc:`LookupError` is thrown.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: t = FSMTransition('A', 'B', 0)
            sage: F = FiniteStateMachine([t])
            sage: F.transition(('A', 'B', 0))
            Transition from 'A' to 'B': 0|-
            sage: id(t) == id(F.transition(('A', 'B', 0)))
            True
        """
        if not isinstance(transition, FSMTransition):
            transition = FSMTransition(*transition)
        for s in self.iter_transitions(transition.from_state):
            if s == transition:
                return s
        raise LookupError("No transition found.")

    # ************************************************************************
    # properties (state and transitions)
    # ************************************************************************

    def has_state(self, state):
        """
        Return whether ``state`` is one of the states of the finite
        state machine.

        INPUT:

        - ``state`` can be a :class:`FSMState` or a label of a state.

        OUTPUT: boolean

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
        Return whether ``transition`` is one of the transitions of
        the finite state machine.

        INPUT:

        - ``transition`` has to be a :class:`FSMTransition`.

        OUTPUT: boolean

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
        if isinstance(transition, FSMTransition):
            return transition in self.iter_transitions()
        raise TypeError("Transition is not an instance of FSMTransition.")

    def has_initial_state(self, state):
        """
        Return whether ``state`` is one of the initial states of the
        finite state machine.

        INPUT:

        - ``state`` can be a :class:`FSMState` or a label.

        OUTPUT: boolean

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
        Return whether the finite state machine has an initial state.

        OUTPUT: boolean

        EXAMPLES::

            sage: FiniteStateMachine().has_initial_states()
            False
        """
        return bool(self.initial_states())

    def has_final_state(self, state):
        """
        Return whether ``state`` is one of the final states of the
        finite state machine.

        INPUT:

        - ``state`` can be a :class:`FSMState` or a label.

        OUTPUT: boolean

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
        Return whether the finite state machine has a final state.

        OUTPUT: boolean

        EXAMPLES::

            sage: FiniteStateMachine().has_final_states()
            False
        """
        return bool(self.final_states())

    # ************************************************************************
    # properties
    # ************************************************************************

    def is_deterministic(self):
        """
        Return whether the finite finite state machine is deterministic.

        OUTPUT: boolean

        A finite state machine is considered to be deterministic if
        each transition has input label of length one and for each
        pair `(q,a)` where `q` is a state and `a` is an element of the
        input alphabet, there is at most one transition from `q` with
        input label `a`. Furthermore, the finite state may not have
        more than one initial state.

        EXAMPLES::

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

        Check that :issue:`18556` is fixed::

            sage: Automaton().is_deterministic()
            True
            sage: Automaton(initial_states=[0]).is_deterministic()
            True
            sage: Automaton(initial_states=[0, 1]).is_deterministic()
            False
        """
        if len(self.initial_states()) > 1:
            return False
        for state in self.iter_states():
            for transition in state.transitions:
                if len(transition.word_in) != 1:
                    return False

            transition_classes_by_word_in = full_group_by(
                state.transitions,
                key=lambda t: t.word_in)

            for _, transition_class in transition_classes_by_word_in:
                if len(transition_class) > 1:
                    return False
        return True

    def is_complete(self):
        """
        Return whether the finite state machine is complete.

        OUTPUT: boolean

        A finite state machine is considered to be complete if
        each transition has an input label of length one and for each
        pair `(q, a)` where `q` is a state and `a` is an element of the
        input alphabet, there is exactly one transition from `q` with
        input label `a`.

        EXAMPLES::

            sage: fsm = FiniteStateMachine([(0, 0, 0, 0),
            ....:                           (0, 1, 1, 1),
            ....:                           (1, 1, 0, 0)],
            ....:                          determine_alphabets=False)
            sage: fsm.is_complete()
            Traceback (most recent call last):
            ...
            ValueError: No input alphabet is given. Try calling determine_alphabets().
            sage: fsm.input_alphabet = [0, 1]
            sage: fsm.is_complete()
            False
            sage: fsm.add_transition((1, 1, 1, 1))
            Transition from 1 to 1: 1|1
            sage: fsm.is_complete()
            True
            sage: fsm.add_transition((0, 0, 1, 0))
            Transition from 0 to 0: 1|0
            sage: fsm.is_complete()
            False
        """
        if self.input_alphabet is None:
            raise ValueError("No input alphabet is given. "
                             "Try calling determine_alphabets().")

        for state in self.iter_states():
            for transition in state.transitions:
                if len(transition.word_in) != 1:
                    return False

            transition_classes_by_word_in = full_group_by(
                state.transitions,
                key=lambda t: t.word_in)

            for key, transition_class in transition_classes_by_word_in:
                if len(transition_class) > 1:
                    return False

            # all input labels are lists, extract the only element
            outgoing_alphabet = [key[0] for key, transition_class in
                                 transition_classes_by_word_in]
            if not sorted(self.input_alphabet) == sorted(outgoing_alphabet):
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

    # ************************************************************************
    # let the finite state machine work
    # ************************************************************************

    _process_default_options_ = {'full_output': True,
                                 'list_of_outputs': None,
                                 'only_accepted': False,
                                 'always_include_output': False,
                                 'automatic_output_type': False}

    def process(self, *args, **kwargs):
        """
        Return whether the finite state machine accepts the input, the state
        where the computation stops and which output is generated.

        INPUT:

        - ``input_tape`` -- the input tape can be a list or an
          iterable with entries from the input alphabet. If we are
          working with a multi-tape machine (see parameter
          ``use_multitape_input`` and notes below), then the tape is a
          list or tuple of tracks, each of which can be a list or an
          iterable with entries from the input alphabet.

        - ``initial_state`` or ``initial_states`` -- the initial
          state(s) in which the machine starts. Either specify a
          single one with ``initial_state`` or a list of them with
          ``initial_states``. If both are given, ``initial_state``
          will be appended to ``initial_states``. If neither is
          specified, the initial states of the finite state machine
          are taken.

        - ``list_of_outputs`` -- (default: ``None``) a boolean or
          ``None``. If ``True``, then the outputs are given in list form
          (even if we have no or only one single output). If
          ``False``, then the result is never a list (an exception is
          raised if the result cannot be returned). If
          ``list_of_outputs=None``, the method determines automatically
          what to do (e.g. if a non-deterministic machine returns more
          than one path, then the output is returned in list form).

        - ``only_accepted`` -- boolean (default: ``False``); if set, then the
          first argument in the output is guaranteed to be ``True`` (if the
          output is a list, then the first argument of each element will be
          ``True``)

        - ``always_include_output`` -- if set (not by default), always
          include the output. This is inconsequential for a
          :class:`FiniteStateMachine`, but can be used in derived
          classes where the output is suppressed by default,
          cf. :meth:`Automaton.process`.

        - ``format_output`` -- a function that translates the written
          output (which is in form of a list) to something more
          readable. By default (``None``) identity is used here.

        - ``check_epsilon_transitions`` -- boolean (default: ``True``); if
          ``False``, then epsilon transitions are not taken into consideration
          during process

        - ``write_final_word_out`` -- boolean (default: ``True``); whether the
          final output words should be written or not

        - ``use_multitape_input`` -- boolean (default: ``False``); if ``True``,
          then the multi-tape mode of the process iterator is activated. See
          also the notes below for multi-tape machines.

        - ``process_all_prefixes_of_input`` -- boolean (default: ``False``); if
          ``True``, then each prefix of the input word is processed (instead of
          processing the whole input word at once). Consequently, there is an
          output generated for each of these prefixes.

        - ``process_iterator_class`` -- (default: ``None``) a class
          inherited from :class:`FSMProcessIterator`. If ``None``,
          then :class:`FSMProcessIterator` is taken. An instance of this
          class is created and is used during the processing.

        - ``automatic_output_type`` -- boolean (default: ``False``); if set and
          the input has a parent, then the output will have the same parent. If
          the input does not have a parent, then the output will be of the same
          type as the input.

        OUTPUT:

        A triple (or a list of triples,
        cf. parameter ``list_of_outputs``), where

        - the first entry is ``True`` if the input string is accepted,

        - the second gives the reached state after processing the
          input tape (This is a state with label ``None`` if the input
          could not be processed, i.e., if at one point no
          transition to go on could be found.), and

        - the third gives a list of the output labels written during
          processing (in the case the finite state machine runs as
          transducer).

        Note that in the case the finite state machine is not
        deterministic, all possible paths are taken into account.

        This function uses an iterator which, in its simplest form, goes
        from one state to another in each step. To decide which way to
        go, it uses the input words of the outgoing transitions and
        compares them to the input tape. More precisely, in each step,
        the iterator takes an outgoing transition of the current state,
        whose input label equals the input letter of the tape. The
        output label of the transition, if present, is written on the
        output tape.

        If the choice of the outgoing transition is not unique (i.e.,
        we have a non-deterministic finite state machine), all
        possibilities are followed. This is done by splitting the
        process into several branches, one for each of the possible
        outgoing transitions.

        The process (iteration) stops if all branches are finished,
        i.e., for no branch, there is any transition whose input word
        coincides with the processed input tape. This can simply
        happen when the entire tape was read.

        Also see :meth:`~FiniteStateMachine.__call__` for a version of
        :meth:`.process` with shortened output.

        Internally this function creates and works with an instance of
        :class:`FSMProcessIterator`. This iterator can also be obtained
        with :meth:`iter_process`.

        If working with multi-tape finite state machines, all input
        words of transitions are words of `k`-tuples of letters.
        Moreover, the input tape has to consist of `k` tracks, i.e.,
        be a list or tuple of `k` iterators, one for each track.

        .. WARNING::

            Working with multi-tape finite state machines is still
            experimental and can lead to wrong outputs.

        EXAMPLES::

            sage: binary_inverter = FiniteStateMachine({'A': [('A', 0, 1), ('A', 1, 0)]},
            ....:                                      initial_states=['A'], final_states=['A'])
            sage: binary_inverter.process([0, 1, 0, 0, 1, 1])
            (True, 'A', [1, 0, 1, 1, 0, 0])

        Alternatively, we can invoke this function by::

            sage: binary_inverter([0, 1, 0, 0, 1, 1])
            (True, 'A', [1, 0, 1, 1, 0, 0])

        Below we construct a finite state machine which tests if an input
        is a non-adjacent form, i.e., no two neighboring letters are
        both nonzero (see also the example on
        :ref:`non-adjacent forms <finite_state_machine_recognizing_NAFs_example>`
        in the documentation of the module
        :doc:`finite_state_machine`)::

            sage: NAF = FiniteStateMachine(
            ....:     {'_': [('_', 0), (1, 1)], 1: [('_', 0)]},
            ....:     initial_states=['_'], final_states=['_', 1])
            sage: [NAF.process(w)[0] for w in [[0], [0, 1], [1, 1], [0, 1, 0, 1],
            ....:                           [0, 1, 1, 1, 0], [1, 0, 0, 1, 1]]]
            [True, True, False, True, False, False]

        Working only with the first component (i.e., returning whether
        accepted or not) usually corresponds to using the more
        specialized class :class:`Automaton`.

        Non-deterministic finite state machines can be handled as well.

        ::

            sage: T = Transducer([(0, 1, 0, 0), (0, 2, 0, 0)],
            ....:     initial_states=[0])
            sage: T.process([0])
            [(False, 1, [0]), (False, 2, [0])]

        Here is another non-deterministic finite state machine. Note
        that we use ``format_output`` (see
        :class:`FSMProcessIterator`) to convert the written outputs
        (all characters) to strings.

        ::

            sage: T = Transducer([(0, 1, [0, 0], 'a'), (0, 2, [0, 0, 1], 'b'),
            ....:                 (0, 1, 1, 'c'), (1, 0, [], 'd'),
            ....:                 (1, 1, 1, 'e')],
            ....:                initial_states=[0], final_states=[0, 1])
            sage: T.process([0], format_output=lambda o: ''.join(o))
            (False, None, None)
            sage: T.process([0, 0], format_output=lambda o: ''.join(o))
            [(True, 0, 'ad'), (True, 1, 'a')]
            sage: T.process([1], format_output=lambda o: ''.join(o))
            [(True, 0, 'cd'), (True, 1, 'c')]
            sage: T.process([1, 1], format_output=lambda o: ''.join(o))
            [(True, 0, 'cdcd'), (True, 0, 'ced'),
             (True, 1, 'cdc'), (True, 1, 'ce')]
            sage: T.process([0, 0, 1], format_output=lambda o: ''.join(o))
            [(False, 2, 'b'),
             (True, 0, 'adcd'),
             (True, 0, 'aed'),
             (True, 1, 'adc'),
             (True, 1, 'ae')]
            sage: T.process([0, 0, 1], format_output=lambda o: ''.join(o),
            ....:           only_accepted=True)
            [(True, 0, 'adcd'), (True, 0, 'aed'),
             (True, 1, 'adc'), (True, 1, 'ae')]

        A simple example of a multi-tape finite state machine is the
        following: It writes the length of the first tape many letters
        ``a`` and then the length of the second tape many letters
        ``b``::

            sage: M = FiniteStateMachine([(0, 0, (1, None), 'a'),
            ....:                         (0, 1, [], []),
            ....:                         (1, 1, (None, 1), 'b')],
            ....:                        initial_states=[0],
            ....:                        final_states=[1])
            sage: M.process(([1, 1], [1]), use_multitape_input=True)
            (True, 1, ['a', 'a', 'b'])

        .. SEEALSO::

            :meth:`Automaton.process`,
            :meth:`Transducer.process`,
            :meth:`~FiniteStateMachine.iter_process`,
            :meth:`~FiniteStateMachine.__call__`,
            :class:`FSMProcessIterator`.

        TESTS::

            sage: T = Transducer([(0, 1, [0, 0], 0), (0, 2, [0, 0, 1], 0),
            ....:                 (0, 1, 1, 2), (1, 0, [], 1), (1, 1, 1, 3)],
            ....:     initial_states=[0], final_states=[0, 1])
            sage: T.process([0])
            (False, None, None)
            sage: T.process([0, 0])
            [(True, 0, [0, 1]), (True, 1, [0])]
            sage: T.process([1])
            [(True, 0, [2, 1]), (True, 1, [2])]
            sage: T.process([1, 1])
            [(True, 0, [2, 1, 2, 1]), (True, 0, [2, 3, 1]),
             (True, 1, [2, 1, 2]), (True, 1, [2, 3])]

        ::

            sage: F = FiniteStateMachine([(0, 0, 0, 0)],
            ....:                        initial_states=[0])
            sage: F.process([0], only_accepted=True)
            []
            sage: F.process([0], only_accepted=True, list_of_outputs=False)
            Traceback (most recent call last):
            ...
            ValueError: No accepting output was found but according to the
            given options, an accepting output should be returned. Change
            only_accepted and/or list_of_outputs options.
            sage: F.process([0], only_accepted=True, list_of_outputs=True)
            []
            sage: F.process([0], only_accepted=False)
            (False, 0, [0])
            sage: F.process([0], only_accepted=False, list_of_outputs=False)
            (False, 0, [0])
            sage: F.process([0], only_accepted=False, list_of_outputs=True)
            [(False, 0, [0])]
            sage: F.process([1], only_accepted=True)
            []
            sage: F.process([1], only_accepted=True, list_of_outputs=False)
            Traceback (most recent call last):
            ...
            ValueError: No accepting output was found but according to the
            given options, an accepting output should be returned. Change
            only_accepted and/or list_of_outputs options.
            sage: F.process([1], only_accepted=True, list_of_outputs=True)
            []
            sage: F.process([1], only_accepted=False)
            (False, None, None)
            sage: F.process([1], only_accepted=False, list_of_outputs=False)
            (False, None, None)
            sage: F.process([1], only_accepted=False, list_of_outputs=True)
            []

        ::

            sage: F = FiniteStateMachine([(0, 1, 1, 'a'), (0, 2, 2, 'b')],
            ....:                        initial_states=[0],
            ....:                        final_states=[1])
            sage: A = Automaton([(0, 1, 1), (0, 2, 2)],
            ....:               initial_states=[0],
            ....:               final_states=[1])
            sage: T = Transducer([(0, 1, 1, 'a'), (0, 2, 2, 'b')],
            ....:                initial_states=[0],
            ....:                final_states=[1])
            sage: F.process([1])
            (True, 1, ['a'])
            sage: A.process([1])
            (True, 1)
            sage: T.process([1])
            (True, 1, ['a'])
            sage: F.process([2])
            (False, 2, ['b'])
            sage: A.process([2])
            (False, 2)
            sage: T.process([2])
            (False, 2, ['b'])
            sage: F.process([3])
            (False, None, None)
            sage: A.process([3])
            (False, None)
            sage: T.process([3])
            (False, None, None)
        """
        # set default values
        options = copy(self._process_default_options_)
        options.update(kwargs)

        # perform iteration
        it = self.iter_process(*args, **options)
        for _ in it:
            pass

        # process output: filtering accepting results
        only_accepted = options['only_accepted']
        it_output = [result for result in it.result()
                     if not only_accepted or result[0]]

        # process output: returning a list output
        if (len(it_output) > 1 and options['list_of_outputs'] is None or
                options['list_of_outputs']):
            return [self._process_convert_output_(out, **options)
                    for out in sorted(it_output)]

        # process output: cannot return output to due input parameters
        if options['list_of_outputs'] is False:
            if not it_output and only_accepted:
                raise ValueError('No accepting output was found but according '
                                 'to the given options, an accepting output '
                                 'should be returned. Change only_accepted '
                                 'and/or list_of_outputs options.')
            elif len(it_output) > 1:
                raise ValueError('Got more than one output, but only allowed '
                                 'to show one. Change list_of_outputs option.')
        # At this point it_output has length 0 or 1.

        # process output: create non-accepting output if needed
        if not it_output:
            if only_accepted:
                return []
            NoneState = FSMState(None, allow_label_None=True)
            it_output = [(False, NoneState, None)]

        return self._process_convert_output_(it_output[0], **options)

    def _process_convert_output_(self, output_data, **kwargs):
        """
        Helper function which converts the output of
        :meth:`FiniteStateMachine.process`. This is the identity.

        INPUT:

        - ``output_data`` -- a triple

        - ``full_output`` -- boolean

        OUTPUT: the converted output

        This function is overridden in :class:`Automaton` and
        :class:`Transducer`.

        TESTS::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: F = FiniteStateMachine()
            sage: F._process_convert_output_((True, FSMState('a'), [1, 0, 1]),
            ....:                            full_output=False)
            (True, 'a', [1, 0, 1])
            sage: F._process_convert_output_((True, FSMState('a'), [1, 0, 1]),
            ....:                            full_output=True)
            (True, 'a', [1, 0, 1])
        """
        accept_input, current_state, output = output_data
        return (accept_input, current_state, output)

    def iter_process(self, input_tape=None, initial_state=None,
                     process_iterator_class=None,
                     iterator_type=None,
                     automatic_output_type=False, **kwargs):
        r"""
        This function returns an iterator for processing the input.
        See :meth:`.process` (which runs this iterator until the end)
        for more information.

        INPUT:

        - ``iterator_type`` -- if ``None`` (default), then
          an instance of :class:`FSMProcessIterator` is returned. If
          this is ``'simple'`` only an iterator over one output is
          returned (an exception is raised if this is not the case, i.e.,
          if the process has branched).

        See :meth:`process` for a description of the other parameters.

        OUTPUT: an iterator

        EXAMPLES:

        We can use :meth:`iter_process` to deal with infinite words::

            sage: inverter = Transducer({'A': [('A', 0, 1), ('A', 1, 0)]},
            ....:     initial_states=['A'], final_states=['A'])

            sage: # needs sage.combinat
            sage: words.FibonacciWord()
            word: 0100101001001010010100100101001001010010...
            sage: it = inverter.iter_process(
            ....:     words.FibonacciWord(), iterator_type='simple')
            sage: Words([0,1])(it)
            word: 1011010110110101101011011010110110101101...

        This can also be done by::

            sage: inverter.iter_process(words.FibonacciWord(),                          # needs sage.combinat
            ....:                       iterator_type='simple',
            ....:                       automatic_output_type=True)
            word: 1011010110110101101011011010110110101101...

        or even simpler by::

            sage: inverter(words.FibonacciWord())                                       # needs sage.combinat
            word: 1011010110110101101011011010110110101101...

        To see what is going on, we use :meth:`iter_process` without
        arguments::

            sage: # needs sage.combinat
            sage: from itertools import islice
            sage: it = inverter.iter_process(words.FibonacciWord())
            sage: for current in islice(it, 4r):
            ....:     print(current)
            process (1 branch)
            + at state 'A'
            +-- tape at 1, [[1]]
            process (1 branch)
            + at state 'A'
            +-- tape at 2, [[1, 0]]
            process (1 branch)
            + at state 'A'
            +-- tape at 3, [[1, 0, 1]]
            process (1 branch)
            + at state 'A'
            +-- tape at 4, [[1, 0, 1, 1]]

        The following show the difference between using the ``'simple'``-option
        and not using it. With this option, we have
        ::

            sage: it = inverter.iter_process(input_tape=[0, 1, 1],
            ....:                            iterator_type='simple')
            sage: for i, o in enumerate(it):
            ....:     print('step %s: output %s' % (i, o))
            step 0: output 1
            step 1: output 0
            step 2: output 0

        So :meth:`iter_process` is a generator expression which gives
        a new output letter in each step (and not more). In many cases
        this is sufficient.

        Doing the same without the ``'simple'``-option does not give
        the output directly; it has to be extracted first. On the
        other hand, additional information is presented::

            sage: it = inverter.iter_process(input_tape=[0, 1, 1])
            sage: for current in it:
            ....:     print(current)
            process (1 branch)
            + at state 'A'
            +-- tape at 1, [[1]]
            process (1 branch)
            + at state 'A'
            +-- tape at 2, [[1, 0]]
            process (1 branch)
            + at state 'A'
            +-- tape at 3, [[1, 0, 0]]
            process (0 branches)
            sage: it.result()
            [Branch(accept=True, state='A', output=[1, 0, 0])]

        One can see the growing of the output (the list of lists at
        the end of each entry).

        Even if the transducer has transitions with empty or multiletter
        output, the simple iterator returns one new output letter in
        each step::

            sage: T = Transducer([(0, 0, 0, []),
            ....:                 (0, 0, 1, [1]),
            ....:                 (0, 0, 2, [2, 2])],
            ....:                initial_states=[0])
            sage: it = T.iter_process(input_tape=[0, 1, 2, 0, 1, 2],
            ....:                     iterator_type='simple')
            sage: for i, o in enumerate(it):
            ....:     print('step %s: output %s' % (i, o))
            step 0: output 1
            step 1: output 2
            step 2: output 2
            step 3: output 1
            step 4: output 2
            step 5: output 2

        .. SEEALSO::

            :meth:`FiniteStateMachine.process`,
            :meth:`Automaton.process`,
            :meth:`Transducer.process`,
            :meth:`~FiniteStateMachine.__call__`,
            :class:`FSMProcessIterator`.
        """
        if automatic_output_type and 'format_output' in kwargs:
            raise ValueError("Parameter 'automatic_output_type' set, but "
                             "'format_output' specified as well.")
        if automatic_output_type:
            try:
                kwargs['format_output'] = input_tape.parent()
            except AttributeError:
                kwargs['format_output'] = type(input_tape)

        if process_iterator_class is None:
            process_iterator_class = FSMProcessIterator
        it = process_iterator_class(self,
                                    input_tape=input_tape,
                                    initial_state=initial_state,
                                    **kwargs)
        if iterator_type is None:
            return it
        elif iterator_type == 'simple':
            simple_it = self._iter_process_simple_(it)
            try:
                return kwargs['format_output'](simple_it)
            except KeyError:
                return simple_it
        else:
            raise ValueError('Iterator type %s unknown.' % (iterator_type,))

    def _iter_process_simple_(self, iterator):
        r"""
        Convert a :class:`process iterator <FSMProcessIterator>` to a simpler
        iterator, which only outputs the written letters.

        INPUT:

        - ``iterator`` -- in instance of :class:`FSMProcessIterator`

        OUTPUT: a generator

        An exception is raised if the process branches.

        EXAMPLES::

            sage: inverter = Transducer({'A': [('A', 0, 1), ('A', 1, 0)]},
            ....:     initial_states=['A'], final_states=['A'])
            sage: it = inverter.iter_process(words.FibonacciWord()[:10])                # needs sage.combinat
            sage: it_simple = inverter._iter_process_simple_(it)                        # needs sage.combinat
            sage: list(it_simple)                                                       # needs sage.combinat
            [1, 0, 1, 1, 0, 1, 0, 1, 1, 0]

        .. SEEALSO::

            :meth:`iter_process`,
            :meth:`FiniteStateMachine.process`,
            :meth:`Automaton.process`,
            :meth:`Transducer.process`,
            :meth:`~FiniteStateMachine.__call__`,
            :class:`FSMProcessIterator`.

        TESTS::

            sage: T = Transducer([(0, 0, [0, 0], 0), (0, 1, 0, 0)],
            ....:                initial_states=[0], final_states=[0])
            sage: list(T.iter_process([0, 0], iterator_type='simple'))
            Traceback (most recent call last):
            ...
            RuntimeError: Process has branched (2 branches exist).
            The 'simple' iterator cannot be used here.
            sage: T = Transducer([(0, 0, 0, 0), (0, 1, 0, 0)],
            ....:                initial_states=[0], final_states=[0])
            sage: list(T.iter_process([0], iterator_type='simple'))
            Traceback (most recent call last):
            ...
            RuntimeError: Process has branched (visiting 2 states in branch).
            The 'simple' iterator cannot be used here.
            sage: T = Transducer([(0, 1, 0, 1), (0, 1, 0, 2)],
            ....:                initial_states=[0], final_states=[0])
            sage: list(T.iter_process([0], iterator_type='simple'))
            Traceback (most recent call last):
            ...
            RuntimeError: Process has branched. (2 different outputs in branch).
            The 'simple' iterator cannot be used here.
        """
        for current in iterator:
            if not current:
                return

            if len(current) > 1:
                raise RuntimeError("Process has branched "
                                   "(%s branches exist). The "
                                   "'simple' iterator cannot be used "
                                   "here." %
                                   (len(current),))
            _, states = next(iter(current.items()))
            if len(states) > 1:
                raise RuntimeError("Process has branched "
                                   "(visiting %s states in branch). The "
                                   "'simple' iterator cannot be used "
                                   "here." %
                                   (len(states),))
            _, branch = next(iter(states.items()))
            if len(branch.outputs) > 1:
                raise RuntimeError("Process has branched. "
                                   "(%s different outputs in branch). The "
                                   "'simple' iterator cannot be used "
                                   "here." %
                                   (len(branch.outputs),))

            yield from branch.outputs[0]
            branch.outputs[0] = []
            # Reset output so that in the next round
            # (of "for current in iterator") only new
            # output is returned (by the yield).

    # ************************************************************************
    # change finite state machine (add/remove state/transitions)
    # ************************************************************************

    def add_state(self, state):
        """
        Add a state to the finite state machine and returns the new
        state. If the state already exists, that existing state is
        returned.

        INPUT:

        - ``state`` is either an instance of
          :class:`FSMState` or,
          otherwise, a label of a state.

        OUTPUT: the new or existing state

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
        if isinstance(state, FSMState):
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
        Add several states. See add_state for more information.

        INPUT:

        - ``states`` -- list of states or iterator over states

        OUTPUT: nothing

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
        Add a transition to the finite state machine and returns the
        new transition.

        If the transition already exists, the return value of
        ``self.on_duplicate_transition`` is returned. See the
        documentation of :class:`FiniteStateMachine`.

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

        If the states ``A`` and ``B`` are not instances of
        :class:`FSMState`, then it is assumed that they are labels of
        states.

        OUTPUT:

        The new transition.
        """
        if len(args) + len(kwargs) == 0:
            return
        if len(args) + len(kwargs) == 1:
            if len(args) == 1:
                d = args[0]
                if isinstance(d, FSMTransition):
                    return self._add_fsm_transition_(d)
            else:
                d = next(iter(kwargs.values()))
            if isinstance(d, Mapping):
                args = []
                kwargs = d
            elif isinstance(d, Iterable):
                args = d
                kwargs = {}
            else:
                raise TypeError("Cannot decide what to do with input.")

        data = dict(zip(
                ('from_state', 'to_state', 'word_in', 'word_out', 'hook'),
                args))
        data.update(kwargs)

        data['from_state'] = self.add_state(data['from_state'])
        data['to_state'] = self.add_state(data['to_state'])

        return self._add_fsm_transition_(FSMTransition(**data))

    def _add_fsm_transition_(self, t):
        """
        Add a transition.

        INPUT:

        - ``t`` -- an instance of :class:`FSMTransition`

        OUTPUT: the new transition

        TESTS::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: F = FiniteStateMachine()
            sage: F._add_fsm_transition_(FSMTransition('A', 'B'))
            Transition from 'A' to 'B': -|-
        """
        try:
            existing_transition = self.transition(t)
        except LookupError:
            pass
        else:
            return self.on_duplicate_transition(existing_transition, t)
        from_state = self.add_state(t.from_state)
        self.add_state(t.to_state)
        from_state.transitions.append(t)
        return t

    def add_from_transition_function(self, function, initial_states=None,
                                     explore_existing_states=True):
        """
        Construct a finite state machine from a transition function.

        INPUT:

        - ``function`` may return a tuple (new_state, output_word) or a
          list of such tuples

        - ``initial_states`` -- if no initial states are given, the
          already existing initial states of ``self`` are taken

        - ``explore_existing_states`` -- boolean (default: ``True``); if
          ``True`` (default), then already existing states in ``self`` (e.g.
          already given final states) will also be processed if they are
          reachable from the initial statess

        OUTPUT: nothing

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
            ....:     return 1 - state, []
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
            ....:     return 1 - state, []
            sage: F.add_from_transition_function(transition_function,
            ....:                                explore_existing_states=False)
            sage: F.transitions()
            [Transition from 0 to 1: 0|-]

        TESTS::

            sage: F = FiniteStateMachine(initial_states=['A'])
            sage: def f(state, input):
            ....:     return [('A', input), ('B', 1-input)]
            sage: F.add_from_transition_function(f)
            Traceback (most recent call last):
            ...
            ValueError: No input alphabet is given.
            Try calling determine_alphabets().

        ::

            sage: def transition(state, where):
            ....:     return (vector([0, 0]), 1)
            sage: Transducer(transition, input_alphabet=[0], initial_states=[0])
            Traceback (most recent call last):
            ...
            TypeError: mutable vectors are unhashable
        """
        if self.input_alphabet is None:
            raise ValueError("No input alphabet is given. "
                               "Try calling determine_alphabets().")

        if initial_states is None:
            not_done = self.initial_states()
        elif isinstance(initial_states, Iterable):
            not_done = []
            for s in initial_states:
                state = self.add_state(s)
                state.is_initial = True
                not_done.append(state)
        else:
            raise TypeError('Initial states must be iterable '
                            '(e.g. a list of states).')
        if not not_done:
            raise ValueError("No state is initial.")
        if explore_existing_states:
            ignore_done = self.states()
            for s in not_done:
                try:
                    ignore_done.remove(s)
                except ValueError:
                    pass
        else:
            ignore_done = []
        while not_done:
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
                        pass
                except TypeError:
                    raise ValueError("The callback function for "
                                     "add_from_transition is expected "
                                     "to return a pair (new_state, "
                                     "output_label) or a list of such pairs. "
                                     "For the state %s and the input "
                                     "letter %s, it however returned %s, "
                                     "which is not acceptable."
                                     % (s.label(), letter, return_value))
                for (st_label, word) in return_value:
                    if not self.has_state(st_label):
                        not_done.append(self.add_state(st_label))
                    elif ignore_done:
                        u = self.state(st_label)
                        if u in ignore_done:
                            not_done.append(u)
                            ignore_done.remove(u)
                    self.add_transition(s, st_label,
                                        word_in=letter, word_out=word)

    def add_transitions_from_function(self, function, labels_as_input=True):
        """
        Add one or more transitions if ``function(state, state)``
        says that there are some.

        INPUT:

        - ``function`` -- a transition function. Given two states
          ``from_state`` and ``to_state`` (or their labels if
          ``label_as_input`` is true), this function shall return a
          tuple ``(word_in, word_out)`` to add a transition from
          ``from_state`` to ``to_state`` with input and output labels
          ``word_in`` and ``word_out``, respectively. If no such
          addition is to be added, the transition function shall
          return ``None``. The transition function may also return
          a list of such tuples in order to add multiple transitions
          between the pair of states.

        - ``label_as_input`` -- (default: ``True``)

        OUTPUT: nothing

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

        Multiple transitions are also possible::

            sage: F = FiniteStateMachine()
            sage: F.add_states([0, 1])
            sage: def f(state1, state2):
            ....:     if state1 != state2:
            ....:          return [(0, 1), (1, 0)]
            ....:     else:
            ....:          return None
            sage: F.add_transitions_from_function(f)
            sage: F.transitions()
            [Transition from 0 to 1: 0|1,
             Transition from 0 to 1: 1|0,
             Transition from 1 to 0: 0|1,
             Transition from 1 to 0: 1|0]

        TESTS::

            sage: F = FiniteStateMachine()
            sage: F.add_state(0)
            0
            sage: def f(state1, state2):
            ....:     return 1
            sage: F.add_transitions_from_function(f)
            Traceback (most recent call last):
            ...
            ValueError: The callback function for add_transitions_from_function
            is expected to return a pair (word_in, word_out) or a list of such
            pairs. For states 0 and 0 however, it returned 1,
            which is not acceptable.
        """
        for s_from in self.iter_states():
            for s_to in self.iter_states():
                try:
                    if labels_as_input:
                        return_value = function(s_from.label(), s_to.label())
                    else:
                        return_value = function(s_from, s_to)
                except LookupError:
                    continue
                if return_value is None:
                    continue
                if not hasattr(return_value, "pop"):
                    transitions = [return_value]
                else:
                    transitions = return_value
                for t in transitions:
                    if not hasattr(t, '__getitem__'):
                        raise ValueError("The callback function for "
                                         "add_transitions_from_function "
                                         "is expected to return a "
                                         "pair (word_in, word_out) or a "
                                         "list of such pairs. For "
                                         "states %s and %s however, it "
                                         "returned %s, which is not "
                                         "acceptable." % (s_from, s_to, return_value))
                    label_in = t[0]
                    try:
                        label_out = t[1]
                    except LookupError:
                        label_out = None
                    self.add_transition(s_from, s_to, label_in, label_out)

    def delete_transition(self, t):
        """
        Delete a transition by removing it from the list of transitions of
        the state, where the transition starts.

        INPUT:

        - ``t`` -- a transition

        OUTPUT: nothing

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
        Delete a state and all transitions coming or going to this state.

        INPUT:

        - ``s`` -- a label of a state or an :class:`FSMState`

        OUTPUT: nothing

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMTransition
            sage: t1 = FSMTransition('A', 'B', 0)
            sage: t2 = FSMTransition('B', 'B', 1)
            sage: F = FiniteStateMachine([t1, t2])
            sage: F.delete_state('A')
            sage: F.transitions()
            [Transition from 'B' to 'B': 1|-]

        TESTS:

        This shows that :issue:`16024` is fixed. ::

            sage: F._states_
            ['B']
            sage: F._states_dict_
            {'B': 'B'}
        """
        state = self.state(s)
        for transition in self.transitions():
            if transition.to_state == state:
                self.delete_transition(transition)
        self._states_.remove(state)
        try:
            del self._states_dict_[state.label()]
        except AttributeError:
            pass

    def remove_epsilon_transitions(self):
        """
        TESTS::

            sage: FiniteStateMachine().remove_epsilon_transitions()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def epsilon_successors(self, state):
        """
        Return the dictionary with states reachable from ``state``
        without reading anything from an input tape as keys. The
        values are lists of outputs.

        INPUT:

        - ``state`` -- the state whose epsilon successors should be
          determined

        OUTPUT: a dictionary mapping states to a list of output words

        The states in the output are the epsilon successors of
        ``state``. Each word of the list of output words is a word
        written when taking a path from ``state`` to the corresponding
        state.

        EXAMPLES::

            sage: T = Transducer([(0, 1, None, 'a'), (1, 2, None, 'b')])
            sage: T.epsilon_successors(0)
            {1: [['a']], 2: [['a', 'b']]}
            sage: T.epsilon_successors(1)
            {2: [['b']]}
            sage: T.epsilon_successors(2)
            {}

        If there is a cycle with only epsilon transitions, then this
        cycle is only processed once and there is no infinite loop::

            sage: S = Transducer([(0, 1, None, 'a'), (1, 0, None, 'b')])
            sage: S.epsilon_successors(0)
            {0: [['a', 'b']], 1: [['a']]}
            sage: S.epsilon_successors(1)
            {0: [['b']], 1: [['b', 'a']]}
        """
        return self.state(state)._epsilon_successors_(self)

    def accessible_components(self):
        """
        Return a new finite state machine with the accessible states
        of ``self`` and all transitions between those states.

        OUTPUT:

        A finite state machine with the accessible states of ``self`` and
        all transitions between those states.

        A state is accessible if there is a directed path from an
        initial state to the state. If ``self`` has no initial states then
        a copy of the finite state machine ``self`` is returned.

        EXAMPLES::

            sage: F = Automaton([(0, 0, 0), (0, 1, 1), (1, 1, 0), (1, 0, 1)],
            ....:               initial_states=[0])
            sage: F.accessible_components()
            Automaton with 2 states

        ::

            sage: F = Automaton([(0, 0, 1), (0, 0, 1), (1, 1, 0), (1, 0, 1)],
            ....:               initial_states=[0])
            sage: F.accessible_components()
            Automaton with 1 state

        .. SEEALSO::
            :meth:`coaccessible_components`

        TESTS:

        Check whether input of length > 1 works::

            sage: F = Automaton([(0, 1, [0, 1]), (0, 2, 0)],
            ....:               initial_states=[0])
            sage: F.accessible_components()
            Automaton with 3 states
        """
        if not self.initial_states():
            return deepcopy(self)
        memo = {}

        def accessible(from_state, read):
            return [(deepcopy(x.to_state, memo), x.word_out)
                    for x in self.iter_transitions(from_state)
                    if x.word_in[0] == read]

        new_initial_states = [deepcopy(x, memo) for x in self.initial_states()]
        result = self.empty_copy()
        result.add_from_transition_function(accessible,
                                            initial_states=new_initial_states)
        for final_state in self.iter_final_states():
            try:
                new_final_state = result.state(final_state.label)
                new_final_state.is_final = True
            except LookupError:
                pass
        return result

    def coaccessible_components(self):
        r"""
        Return the sub-machine induced by the coaccessible states of this
        finite state machine.

        OUTPUT:

        A finite state machine of the same type as this finite state
        machine.

        EXAMPLES::

            sage: A = automata.ContainsWord([1, 1],
            ....:     input_alphabet=[0, 1]).complement().minimization().relabeled()
            sage: A.transitions()
            [Transition from 0 to 0: 0|-,
             Transition from 0 to 0: 1|-,
             Transition from 1 to 2: 0|-,
             Transition from 1 to 0: 1|-,
             Transition from 2 to 2: 0|-,
             Transition from 2 to 1: 1|-]
            sage: A.initial_states()
            [2]
            sage: A.final_states()
            [1, 2]
            sage: C = A.coaccessible_components()
            sage: C.transitions()
            [Transition from 1 to 2: 0|-,
             Transition from 2 to 2: 0|-,
             Transition from 2 to 1: 1|-]

        .. SEEALSO::
            :meth:`accessible_components`,
            :meth:`induced_sub_finite_state_machine`
        """
        DG = self.digraph().reverse()
        coaccessible_states = DG.breadth_first_search(
            [_.label() for _ in self.iter_final_states()])
        return self.induced_sub_finite_state_machine(
            [self.state(_) for _ in coaccessible_states])

    # *************************************************************************
    # creating new finite state machines
    # *************************************************************************

    def disjoint_union(self, other):
        """
        Return the disjoint union of this and another finite state
        machine.

        INPUT:

        -   ``other`` -- a :class:`FiniteStateMachine`

        OUTPUT:

        A finite state machine of the same type as this finite state
        machine.

        In general, the disjoint union of two finite state machines is
        non-deterministic. In the case of a automata, the language
        accepted by the disjoint union is the union of the languages
        accepted by the constituent automata. In the case of
        transducer, for each successful path in one of the constituent
        transducers, there will be one successful path with the same input
        and output labels in the disjoint union.

        The labels of the states of the disjoint union are pairs ``(i,
        s)``: for each state ``s`` of this finite state machine, there
        is a state ``(0, s)`` in the disjoint union; for each state
        ``s`` of the other finite state machine, there is a state ``(1,
        s)`` in the disjoint union.

        The input alphabet is the union of the input alphabets (if
        possible) and ``None`` otherwise. In the latter case, try
        calling :meth:`.determine_alphabets`.

        The disjoint union can also be written as ``A + B`` or ``A | B``.

        EXAMPLES::

            sage: A = Automaton([(0, 1, 0), (1, 0, 1)],
            ....:               initial_states=[0],
            ....:               final_states=[0])
            sage: A([0, 1, 0, 1])
            True
            sage: B = Automaton([(0, 1, 0), (1, 2, 0), (2, 0, 1)],
            ....:               initial_states=[0],
            ....:               final_states=[0])
            sage: B([0, 0, 1])
            True
            sage: C = A.disjoint_union(B)
            sage: C
            Automaton with 5 states
            sage: C.transitions()
            [Transition from (0, 0) to (0, 1): 0|-,
             Transition from (0, 1) to (0, 0): 1|-,
             Transition from (1, 0) to (1, 1): 0|-,
             Transition from (1, 1) to (1, 2): 0|-,
             Transition from (1, 2) to (1, 0): 1|-]
            sage: C([0, 0, 1])
            True
            sage: C([0, 1, 0, 1])
            True
            sage: C([1])
            False
            sage: C.initial_states()
            [(0, 0), (1, 0)]

        Instead of ``.disjoint_union``, alternative notations are
        available::

            sage: C1 = A + B
            sage: C1 == C
            True
            sage: C2 = A | B
            sage: C2 == C
            True

        In general, the disjoint union is not deterministic.::

            sage: C.is_deterministic()
            False
            sage: D = C.determinisation().minimization()
            sage: D.is_equivalent(Automaton([(0, 0, 0), (0, 0, 1),
            ....:    (1, 7, 0), (1, 0, 1), (2, 6, 0), (2, 0, 1),
            ....:    (3, 5, 0), (3, 0, 1), (4, 0, 0), (4, 2, 1),
            ....:    (5, 0, 0), (5, 3, 1), (6, 4, 0), (6, 0, 1),
            ....:    (7, 4, 0), (7, 3, 1)],
            ....:    initial_states=[1],
            ....:    final_states=[1, 2, 3]))
            True

        Disjoint union of transducers::

            sage: T1 = Transducer([(0, 0, 0, 1)],
            ....:                 initial_states=[0],
            ....:                 final_states=[0])
            sage: T2 = Transducer([(0, 0, 0, 2)],
            ....:                 initial_states=[0],
            ....:                 final_states=[0])
            sage: T1([0])
            [1]
            sage: T2([0])
            [2]
            sage: T = T1.disjoint_union(T2)
            sage: T([0])
            Traceback (most recent call last):
            ...
            ValueError: Found more than one accepting path.
            sage: T.process([0])
            [(True, (0, 0), [1]), (True, (1, 0), [2])]

        Handling of the input alphabet (see :issue:`18989`)::

            sage: A = Automaton([(0, 0, 0)])
            sage: B = Automaton([(0, 0, 1)], input_alphabet=[1, 2])
            sage: C = Automaton([(0, 0, 2)], determine_alphabets=False)
            sage: D = Automaton([(0, 0, [[0, 0]])], input_alphabet=[[0, 0]])
            sage: A.input_alphabet
            [0]
            sage: B.input_alphabet
            [1, 2]
            sage: C.input_alphabet is None
            True
            sage: D.input_alphabet
            [[0, 0]]
            sage: (A + B).input_alphabet
            [0, 1, 2]
            sage: (A + C).input_alphabet is None
            True
            sage: (A + D).input_alphabet is None
            True

        .. SEEALSO::

            :meth:`Automaton.intersection`,
            :meth:`Transducer.intersection`,
            :meth:`.determine_alphabets`.
        """
        result = self.empty_copy()
        for s in self.iter_states():
            result.add_state(s.relabeled((0, s)))
        for s in other.iter_states():
            result.add_state(s.relabeled((1, s)))
        for t in self.iter_transitions():
            result.add_transition((0, t.from_state),
                                  (0, t.to_state),
                                  t.word_in,
                                  t.word_out)
        for t in other.iter_transitions():
            result.add_transition((1, t.from_state),
                                  (1, t.to_state),
                                  t.word_in,
                                  t.word_out)
        try:
            result.input_alphabet = list(set(self.input_alphabet)
                                         | set(other.input_alphabet))
        except TypeError:
            # e.g. None or unhashable letters
            result.input_alphabet = None

        return result

    def concatenation(self, other):
        r"""
        Concatenate this finite state machine with another finite
        state machine.

        INPUT:

        - ``other`` -- a :class:`FiniteStateMachine`

        OUTPUT:

        A :class:`FiniteStateMachine` of the same type as this finite
        state machine.

        Assume that both finite state machines are automata. If
        `\mathcal{L}_1` is the language accepted by this automaton and
        `\mathcal{L}_2` is the language accepted by the other automaton,
        then the language accepted by the concatenated automaton is
        `\{ w_1w_2 \mid w_1\in\mathcal{L}_1, w_2\in\mathcal{L}_2\}` where
        `w_1w_2` denotes the concatenation of the words `w_1` and `w_2`.

        Assume that both finite state machines are transducers and that
        this transducer maps words `w_1\in\mathcal{L}_1` to words
        `f_1(w_1)` and that the other transducer maps words
        `w_2\in\mathcal{L}_2` to words `f_2(w_2)`. Then the concatenated
        transducer maps words `w_1w_2` with `w_1\in\mathcal{L}_1` and
        `w_2\in\mathcal{L}_2` to `f_1(w_1)f_2(w_2)`. Here, `w_1w_2` and
        `f_1(w_1)f_2(w_2)` again denote concatenation of words.

        The input alphabet is the union of the input alphabets (if
        possible) and ``None`` otherwise. In the latter case, try
        calling :meth:`.determine_alphabets`.

        Instead of ``A.concatenation(B)``, the notation ``A * B`` can be
        used.

        EXAMPLES:

        Concatenation of two automata::

            sage: A = automata.Word([0])
            sage: B = automata.Word([1])
            sage: C = A.concatenation(B)
            sage: C.transitions()
            [Transition from (0, 0) to (0, 1): 0|-,
             Transition from (0, 1) to (1, 0): -|-,
             Transition from (1, 0) to (1, 1): 1|-]
            sage: [w
            ....:  for w in ([0, 0], [0, 1], [1, 0], [1, 1])
            ....:  if C(w)]
            [[0, 1]]
            sage: from sage.combinat.finite_state_machine import (
            ....:     Automaton, Transducer)
            sage: isinstance(C, Automaton)
            True

        Concatenation of two transducers::

            sage: A = Transducer([(0, 1, 0, 1), (0, 1, 1, 2)],
            ....:                initial_states=[0],
            ....:                final_states=[1])
            sage: B = Transducer([(0, 1, 0, 1), (0, 1, 1, 0)],
            ....:                initial_states=[0],
            ....:                final_states=[1])
            sage: C = A.concatenation(B)
            sage: C.transitions()
            [Transition from (0, 0) to (0, 1): 0|1,
             Transition from (0, 0) to (0, 1): 1|2,
             Transition from (0, 1) to (1, 0): -|-,
             Transition from (1, 0) to (1, 1): 0|1,
             Transition from (1, 0) to (1, 1): 1|0]
            sage: [(w, C(w)) for w in ([0, 0], [0, 1], [1, 0], [1, 1])]
            [([0, 0], [1, 1]),
             ([0, 1], [1, 0]),
             ([1, 0], [2, 1]),
             ([1, 1], [2, 0])]
            sage: isinstance(C, Transducer)
            True


        Alternative notation as multiplication::

            sage: C == A * B
            True

        Final output words are taken into account::

            sage: A = Transducer([(0, 1, 0, 1)],
            ....:                initial_states=[0],
            ....:                final_states=[1])
            sage: A.state(1).final_word_out = 2
            sage: B = Transducer([(0, 1, 0, 3)],
            ....:                initial_states=[0],
            ....:                final_states=[1])
            sage: B.state(1).final_word_out = 4
            sage: C = A * B
            sage: C([0, 0])
            [1, 2, 3, 4]

        Handling of the input alphabet::

            sage: A = Automaton([(0, 0, 0)])
            sage: B = Automaton([(0, 0, 1)], input_alphabet=[1, 2])
            sage: C = Automaton([(0, 0, 2)], determine_alphabets=False)
            sage: D = Automaton([(0, 0, [[0, 0]])], input_alphabet=[[0, 0]])
            sage: A.input_alphabet
            [0]
            sage: B.input_alphabet
            [1, 2]
            sage: C.input_alphabet is None
            True
            sage: D.input_alphabet
            [[0, 0]]
            sage: (A * B).input_alphabet
            [0, 1, 2]
            sage: (A * C).input_alphabet is None
            True
            sage: (A * D).input_alphabet is None
            True

        .. SEEALSO::

            :meth:`~.disjoint_union`,
            :meth:`.determine_alphabets`.

        TESTS::

            sage: A = Automaton()
            sage: F = FiniteStateMachine()
            sage: A * F
            Traceback (most recent call last):
            ...
            TypeError: Cannot concatenate finite state machines of
            different types.
            sage: F * A
            Traceback (most recent call last):
            ...
            TypeError: Cannot concatenate finite state machines of
            different types.
            sage: F * 5
            Traceback (most recent call last):
            ...
            TypeError: A finite state machine can only be concatenated
            with a another finite state machine.
        """
        if not isinstance(other, FiniteStateMachine):
            raise TypeError('A finite state machine can only be concatenated '
                            'with a another finite state machine.')
        if isinstance(other, Automaton) != isinstance(self, Automaton):
            raise TypeError('Cannot concatenate finite state machines of '
                            'different types.')

        result = self.empty_copy()
        first_states = {}
        second_states = {}
        for s in self.iter_states():
            new_state = s.relabeled((0, s.label()))
            new_state.final_word_out = None
            new_state.is_final = False
            first_states[s] = new_state
            result.add_state(new_state)

        for s in other.iter_states():
            new_state = s.relabeled((1, s.label()))
            new_state.is_initial = False
            second_states[s] = new_state
            result.add_state(new_state)

        for t in self.iter_transitions():
            result.add_transition(first_states[t.from_state],
                                  first_states[t.to_state],
                                  t.word_in,
                                  t.word_out)

        for t in other.iter_transitions():
            result.add_transition(second_states[t.from_state],
                                  second_states[t.to_state],
                                  t.word_in,
                                  t.word_out)

        for s in self.iter_final_states():
            first_state = first_states[s]
            for t in other.iter_initial_states():
                second_state = second_states[t]
                result.add_transition(first_state,
                                      second_state,
                                      [],
                                      s.final_word_out)

        try:
            result.input_alphabet = list(set(self.input_alphabet)
                                         | set(other.input_alphabet))
        except TypeError:
            # e.g. None or unhashable letters
            result.input_alphabet = None

        return result

    __mul__ = concatenation

    def kleene_star(self):
        r"""
        Compute the Kleene closure of this finite state machine.

        OUTPUT:

        A :class:`FiniteStateMachine` of the same type as this finite
        state machine.

        Assume that this finite state machine is an automaton
        recognizing the language `\mathcal{L}`.  Then the Kleene star
        recognizes the language `\mathcal{L}^*=\{ w_1\ldots w_n \mid
        n\ge 0, w_j\in\mathcal{L} \text{ for all } j\}`.

        Assume that this finite state machine is a transducer realizing
        a function `f` on some alphabet `\mathcal{L}`. Then the Kleene
        star realizes a function `g` on `\mathcal{L}^*` with
        `g(w_1\ldots w_n)=f(w_1)\ldots f(w_n)`.

        EXAMPLES:

        Kleene star of an automaton::

            sage: A = automata.Word([0, 1])
            sage: B = A.kleene_star()
            sage: B.transitions()
            [Transition from 0 to 1: 0|-,
             Transition from 2 to 0: -|-,
             Transition from 1 to 2: 1|-]
            sage: from sage.combinat.finite_state_machine import (
            ....:     is_Automaton, is_Transducer)
            sage: isinstance(B, Automaton)
            True
            sage: [w for w in ([], [0, 1], [0, 1, 0], [0, 1, 0, 1], [0, 1, 1, 1])
            ....:  if B(w)]
            [[],
             [0, 1],
             [0, 1, 0, 1]]

        Kleene star of a transducer::

            sage: T = Transducer([(0, 1, 0, 1), (0, 1, 1, 0)],
            ....:                initial_states=[0],
            ....:                final_states=[1])
            sage: S = T.kleene_star()
            sage: S.transitions()
            [Transition from 0 to 1: 0|1,
             Transition from 0 to 1: 1|0,
             Transition from 1 to 0: -|-]
            sage: isinstance(S, Transducer)
            True
            sage: for w in ([], [0], [1], [0, 0], [0, 1]):
            ....:     print("{} {}".format(w, S.process(w)))
            []     (True, 0, [])
            [0]    [(True, 0, [1]), (True, 1, [1])]
            [1]    [(True, 0, [0]), (True, 1, [0])]
            [0, 0] [(True, 0, [1, 1]), (True, 1, [1, 1])]
            [0, 1] [(True, 0, [1, 0]), (True, 1, [1, 0])]

        Final output words are taken into account::

            sage: T = Transducer([(0, 1, 0, 1)],
            ....:                initial_states=[0],
            ....:                final_states=[1])
            sage: T.state(1).final_word_out = 2
            sage: S = T.kleene_star()
            sage: sorted(S.process([0, 0]))
            [(True, 0, [1, 2, 1, 2]), (True, 1, [1, 2, 1, 2])]

        Final output words may lead to undesirable situations if initial
        states and final states coincide::

            sage: T = Transducer(initial_states=[0], final_states=[0])
            sage: T.state(0).final_word_out = 1
            sage: T([])
            [1]
            sage: S = T.kleene_star()
            sage: S([])
            Traceback (most recent call last):
            ...
            RuntimeError: State 0 is in an epsilon cycle (no input), but
            output is written.
        """
        result = deepcopy(self)
        for initial in result.iter_initial_states():
            for final in result.iter_final_states():
                result.add_transition(final, initial, [], final.final_word_out)

        for initial in result.iter_initial_states():
            initial.is_final = True

        return result

    def intersection(self, other):
        """
        TESTS::

            sage: FiniteStateMachine().intersection(FiniteStateMachine())
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def product_FiniteStateMachine(self, other, function,
                                   new_input_alphabet=None,
                                   only_accessible_components=True,
                                   final_function=None,
                                   new_class=None):
        r"""
        Return a new finite state machine whose states are
        `d`-tuples of states of the original finite state machines.

        INPUT:

        - ``other`` -- a finite state machine (for `d=2`) or a list
          (or iterable) of `d-1` finite state machines

        - ``function`` has to accept `d` transitions from `A_j` to `B_j`
          for `j\in\{1, \ldots, d\}` and returns a pair ``(word_in, word_out)``
          which is the label of the transition `A=(A_1, \ldots, A_d)` to `B=(B_1,
          \ldots, B_d)`. If there is no transition from `A` to `B`,
          then ``function`` should raise a :exc:`LookupError`.

        - ``new_input_alphabet`` -- (optional) the new input alphabet
          as a list

        - ``only_accessible_components`` -- if ``True`` (default), then
          the result is piped through :meth:`.accessible_components`. If no
          ``new_input_alphabet`` is given, it is determined by
          :meth:`.determine_alphabets`.

        - ``final_function`` -- a function mapping `d` final states of
          the original finite state machines to the final output of
          the corresponding state in the new finite state machine. By
          default, the final output is the empty word if both final
          outputs of the constituent states are empty; otherwise, a
          :exc:`ValueError` is raised.

        - ``new_class`` -- class of the new finite state machine. By
          default (``None``), the class of ``self`` is used

        OUTPUT:

        A finite state machine whose states are `d`-tuples of states of the
        original finite state machines. A state is initial or
        final if all constituent states are initial or final,
        respectively.

        The labels of the transitions are defined by ``function``.

        The final output of a final state is determined by calling
        ``final_function`` on the constituent states.

        The color of a new state is the tuple of colors of the
        constituent states of ``self`` and ``other``. However,
        if all constituent states have color ``None``, then
        the state has color ``None``, too.

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
            sage: [s.color for s in H.iter_states()]
            [None, None]
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

        Also final output words are considered according to the function
        ``final_function``::

            sage: F = Transducer([(0, 1, 0, 1), (1, 1, 1, 1), (1, 1, 0, 1)],
            ....:                final_states=[1])
            sage: F.state(1).final_word_out = 1
            sage: G = Transducer([(0, 0, 0, 1), (0, 0, 1, 0)], final_states=[0])
            sage: G.state(0).final_word_out = 1
            sage: def minus(t1, t2):
            ....:     return (t1.word_in[0] - t2.word_in[0],
            ....:                t1.word_out[0] - t2.word_out[0])
            sage: H = F.product_FiniteStateMachine(G, minus)
            Traceback (most recent call last):
            ...
            ValueError: A final function must be given.
            sage: def plus(s1, s2):
            ....:     return s1.final_word_out[0] + s2.final_word_out[0]
            sage: H = F.product_FiniteStateMachine(G, minus,
            ....:                                  final_function=plus)
            sage: H.final_states()
            [(1, 0)]
            sage: H.final_states()[0].final_word_out
            [2]

        Products of more than two finite state machines are also possible::

            sage: def plus(s1, s2, s3):
            ....:     if s1.word_in == s2.word_in == s3.word_in:
            ....:          return (s1.word_in,
            ....:                  sum(s.word_out[0] for s in (s1, s2, s3)))
            ....:     else:
            ....:         raise LookupError
            sage: T0 = transducers.CountSubblockOccurrences([0, 0], [0, 1, 2])
            sage: T1 = transducers.CountSubblockOccurrences([1, 1], [0, 1, 2])
            sage: T2 = transducers.CountSubblockOccurrences([2, 2], [0, 1, 2])
            sage: T = T0.product_FiniteStateMachine([T1, T2], plus)
            sage: T.transitions()
            [Transition from ((), (), ()) to ((0,), (), ()): 0|0,
             Transition from ((), (), ()) to ((), (1,), ()): 1|0,
             Transition from ((), (), ()) to ((), (), (2,)): 2|0,
             Transition from ((0,), (), ()) to ((0,), (), ()): 0|1,
             Transition from ((0,), (), ()) to ((), (1,), ()): 1|0,
             Transition from ((0,), (), ()) to ((), (), (2,)): 2|0,
             Transition from ((), (1,), ()) to ((0,), (), ()): 0|0,
             Transition from ((), (1,), ()) to ((), (1,), ()): 1|1,
             Transition from ((), (1,), ()) to ((), (), (2,)): 2|0,
             Transition from ((), (), (2,)) to ((0,), (), ()): 0|0,
             Transition from ((), (), (2,)) to ((), (1,), ()): 1|0,
             Transition from ((), (), (2,)) to ((), (), (2,)): 2|1]
            sage: T([0, 0, 1, 1, 2, 2, 0, 1, 2, 2])
            [0, 1, 0, 1, 0, 1, 0, 0, 0, 1]

        ``other`` can also be an iterable::

            sage: T == T0.product_FiniteStateMachine(iter([T1, T2]), plus)
            True

        TESTS:

        Check that colors are correctly dealt with. In particular, the
        new colors have to be hashable such that
        :meth:`Automaton.determinisation` does not fail::

            sage: A = Automaton([[0, 0, 0]], initial_states=[0])
            sage: B = A.product_FiniteStateMachine(A,
            ....:                                  lambda t1, t2: (0, None))
            sage: B.states()[0].color is None
            True
            sage: B.determinisation()
            Automaton with 1 state

        Check handling of the parameter ``other``::

            sage: A.product_FiniteStateMachine(None, plus)
            Traceback (most recent call last):
            ...
            ValueError: other must be a finite state machine or a list
            of finite state machines.
            sage: A.product_FiniteStateMachine([None], plus)
            Traceback (most recent call last):
            ...
            ValueError: other must be a finite state machine or a list
            of finite state machines.

        Test whether ``new_class`` works::

            sage: T = Transducer()
            sage: type(T.product_FiniteStateMachine(T, None))
            <class 'sage.combinat.finite_state_machine.Transducer'>
            sage: type(T.product_FiniteStateMachine(T, None,
            ....:      new_class=Automaton))
            <class 'sage.combinat.finite_state_machine.Automaton'>

        Check that isolated vertices are kept (:issue:`16762`)::

            sage: F = Transducer(initial_states=[0])
            sage: F.add_state(1)
            1
            sage: G = Transducer(initial_states=['A'])
            sage: F.product_FiniteStateMachine(G, None).states()
            [(0, 'A')]
            sage: F.product_FiniteStateMachine(
            ....:     G, None, only_accessible_components=False).states()
            [(0, 'A'), (1, 'A')]
        """
        def default_final_function(*args):
            if any(s.final_word_out for s in args):
                raise ValueError("A final function must be given.")
            return []

        if final_function is None:
            final_function = default_final_function

        result = self.empty_copy(new_class=new_class)
        if new_input_alphabet is not None:
            result.input_alphabet = new_input_alphabet
        else:
            result.input_alphabet = None

        if isinstance(other, Iterable):
            machines = [self]
            machines.extend(other)
            if not all(isinstance(m, FiniteStateMachine) for m in machines):
                raise ValueError("other must be a finite state machine "
                                 "or a list of finite state machines.")
        elif isinstance(other, FiniteStateMachine):
            machines = [self, other]
        else:
            raise ValueError("other must be a finite state machine or "
                             "a list of finite state machines.")

        for transitions in itertools.product(
            *(m.iter_transitions() for m in machines)):
            try:
                word = function(*transitions)
            except LookupError:
                continue
            result.add_transition(tuple(t.from_state for t in transitions),
                                  tuple(t.to_state for t in transitions),
                                  word[0], word[1])

        if only_accessible_components:
            state_iterator = itertools.product(
                *(m.iter_initial_states() for m in machines))
        else:
            state_iterator = itertools.product(
                *(m.iter_states() for m in machines))

        for state in state_iterator:
            result.add_state(state)

        for state in result.states():
            if all(s.is_initial for s in state.label()):
                state.is_initial = True
            if all(s.is_final for s in state.label()):
                state.is_final = True
                state.final_word_out = final_function(*state.label())
            if all(s.color is None for s in state.label()):
                state.color = None
            else:
                state.color = tuple(s.color for s in state.label())

        if only_accessible_components:
            if result.input_alphabet is None:
                result.determine_alphabets()
            return result.accessible_components()
        else:
            return result

    def composition(self, other, algorithm=None,
                    only_accessible_components=True):
        """
        Return a new transducer which is the composition of ``self``
        and ``other``.

        INPUT:

        - ``other`` -- a transducer

        - ``algorithm`` -- can be one of the following

          - ``direct`` -- the composition is calculated directly

            There can be arbitrarily many initial and final states,
            but the input and output labels must have length `1`.

            .. WARNING::

                The output of ``other`` is fed into ``self``.

          - ``explorative`` -- an explorative algorithm is used

            The input alphabet of ``self`` has to be specified.

            .. WARNING::

                The output of ``other`` is fed into ``self``.

          If algorithm is ``None``, then the algorithm is chosen
          automatically (at the moment always ``direct``, except when
          there are output words of ``other`` or input words of ``self``
          of length greater than `1`).

        OUTPUT: a new transducer

        The labels of the new finite state machine are pairs of states
        of the original finite state machines. The color of a new
        state is the tuple of colors of the constituent states.

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
            sage: He = F.composition(G, algorithm='explorative')
            sage: He.initial_states()
            [(1, 'A'), (1, 'B')]
            sage: He.transitions()
            [Transition from (1, 'A') to (2, 'B'): 0|0,
             Transition from (1, 'B') to (1, 'A'): 1|1,
             Transition from (2, 'B') to (2, 'A'): 0|1,
             Transition from (2, 'A') to (2, 'B'): 1|0]
            sage: Hd == He
            True

        The following example has output of length `> 1`, so the
        explorative algorithm has to be used (and is selected
        automatically).

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
            sage: Ha = G.composition(F)
            sage: Ha == He
            True

        Final output words are also considered::

            sage: F = Transducer([('A', 'B', 1, 0), ('B', 'A', 0, 1)],
            ....:                initial_states=['A', 'B'],
            ....:                final_states=['A', 'B'])
            sage: F.state('A').final_word_out = 0
            sage: F.state('B').final_word_out = 1
            sage: G = Transducer([(1, 1, 1, 0), (1, 2, 0, 1),
            ....:                 (2, 2, 1, 1), (2, 2, 0, 0)],
            ....:                initial_states=[1], final_states=[2])
            sage: G.state(2).final_word_out = 0
            sage: Hd = F.composition(G, algorithm='direct')
            sage: Hd.final_states()
            [(2, 'B')]
            sage: He = F.composition(G, algorithm='explorative')
            sage: He.final_states()
            [(2, 'B')]

        Note that ``(2, 'A')`` is not final, as the final output `0`
        of state `2` of `G` cannot be processed in state ``'A'`` of
        `F`.

        ::

            sage: [s.final_word_out for s in Hd.final_states()]
            [[1, 0]]
            sage: [s.final_word_out for s in He.final_states()]
            [[1, 0]]
            sage: Hd == He
            True

        Here is a non-deterministic example with intermediate output
        length `>1`.

        ::

            sage: F = Transducer([(1, 1, 1, ['a', 'a']), (1, 2, 1, 'b'),
            ....:                 (2, 1, 2, 'a'), (2, 2, 2, 'b')],
            ....:                initial_states=[1, 2])
            sage: G = Transducer([('A', 'A', 'a', 'i'),
            ....:                 ('A', 'B', 'a', 'l'),
            ....:                 ('B', 'B', 'b', 'e')],
            ....:                initial_states=['A', 'B'])
            sage: G(F).transitions()
            [Transition from (1, 'A') to (1, 'A'): 1|'i','i',
             Transition from (1, 'A') to (1, 'B'): 1|'i','l',
             Transition from (1, 'B') to (2, 'B'): 1|'e',
             Transition from (2, 'A') to (1, 'A'): 2|'i',
             Transition from (2, 'A') to (1, 'B'): 2|'l',
             Transition from (2, 'B') to (2, 'B'): 2|'e']

        Be aware that after composition, different transitions may
        share the same output label (same python object)::

            sage: F = Transducer([ ('A','B',0,0), ('B','A',0,0)],
            ....:                initial_states=['A'],
            ....:                final_states=['A'])
            sage: F.transitions()[0].word_out is F.transitions()[1].word_out
            False
            sage: G = Transducer([('C','C',0,1)],
            ....:                initial_states=['C'],
            ....:                final_states=['C'])
            sage: H = G.composition(F)
            sage: H.transitions()[0].word_out is H.transitions()[1].word_out
            True

        TESTS:

        In the explorative algorithm, transducers with non-empty final
        output words are implemented in :issue:`16548`::

            sage: A = transducers.GrayCode()
            sage: B = transducers.abs([0, 1])
            sage: A.composition(B, algorithm='explorative').transitions()
            [Transition from (0, 0) to (0, 1): 0|-,
             Transition from (0, 0) to (0, 2): 1|-,
             Transition from (0, 1) to (0, 1): 0|0,
             Transition from (0, 1) to (0, 2): 1|1,
             Transition from (0, 2) to (0, 1): 0|1,
             Transition from (0, 2) to (0, 2): 1|0]

        Similarly, the explorative algorithm can handle
        non-deterministic finite state machines as of :issue:`16548`::

            sage: A = Transducer([(0, 0, 0, 0), (0, 1, 0, 0)],
            ....:                initial_states=[0])
            sage: B = transducers.Identity([0])
            sage: A.composition(B, algorithm='explorative').transitions()
            [Transition from (0, 0) to (0, 0): 0|0,
             Transition from (0, 0) to (0, 1): 0|0]
            sage: B.composition(A, algorithm='explorative').transitions()
            [Transition from (0, 0) to (0, 0): 0|0,
             Transition from (0, 0) to (1, 0): 0|0]

        In the following example, ``algorithm='direct'`` is inappropriate
        as there are edges with output labels of length greater than 1::

            sage: F = Transducer([('A', 'B', 1, [1, 0]), ('B', 'B', 1, 1),
            ....:                 ('B', 'B', 0, 0)],
            ....:                initial_states=['A'], final_states=['B'])
            sage: G = Transducer([(1, 1, 0, 0), (1, 2, 1, 0),
            ....:                 (2, 2, 0, 1), (2, 1, 1, 1)],
            ....:                initial_states=[1], final_states=[1])
            sage: Hd = G.composition(F, algorithm='direct')

        In the following examples, we compose transducers and automata
        and check whether the types are correct.

        ::

            sage: from sage.combinat.finite_state_machine import (
            ....:     is_Automaton, is_Transducer)
            sage: T = Transducer([(0, 0, 0, 0)], initial_states=[0])
            sage: A = Automaton([(0, 0, 0)], initial_states=[0])
            sage: isinstance(T.composition(T, algorithm='direct'), Transducer)
            True
            sage: isinstance(T.composition(T, algorithm='explorative'), Transducer)
            True
            sage: T.composition(A, algorithm='direct')
            Traceback (most recent call last):
            ...
            TypeError: Composition with automaton is not possible.
            sage: T.composition(A, algorithm='explorative')
            Traceback (most recent call last):
            ...
            TypeError: Composition with automaton is not possible.
            sage: A.composition(A, algorithm='direct')
            Traceback (most recent call last):
            ...
            TypeError: Composition with automaton is not possible.
            sage: A.composition(A, algorithm='explorative')
            Traceback (most recent call last):
            ...
            TypeError: Composition with automaton is not possible.
            sage: isinstance(A.composition(T, algorithm='direct'), Automaton)
            True
            sage: isinstance(A.composition(T, algorithm='explorative'), Automaton)
            True

        Non-deterministic final output cannot be handled::

            sage: F = Transducer([('I', 'A', 0, 42), ('I', 'B', 0, 42)],
            ....:                initial_states=['I'],
            ....:                final_states=['A', 'B'])
            sage: G = Transducer(initial_states=[0],
            ....:                final_states=[0],
            ....:                input_alphabet=[0])
            sage: G.state(0).final_word_out = 0
            sage: H = F.composition(G, algorithm='explorative')
            sage: for s in H.final_states():
            ....:     print("{} {}".format(s, s.final_word_out))
            (0, 'I') [42]
            sage: F.state('A').final_word_out = 'a'
            sage: F.state('B').final_word_out = 'b'
            sage: F.composition(G, algorithm='explorative')
            Traceback (most recent call last):
            ...
            NotImplementedError: Stopping in state (0, 'I') leads to
            non-deterministic final output.

        Check that the output and input alphabets are set correctly::

            sage: F = Transducer([(0, 0, 1, 'A')],
            ....:                initial_states=[0],
            ....:                determine_alphabets=False)
            sage: G = Transducer([(2, 2, 'A', 'a')],
            ....:                initial_states=[2],
            ....:                determine_alphabets=False)
            sage: Hd = G(F, algorithm='direct')
            sage: Hd.input_alphabet, Hd.output_alphabet
            ([1], ['a'])
            sage: He = G(F, algorithm='explorative')
            Traceback (most recent call last):
            ...
            ValueError: No input alphabet is given. Try calling
            determine_alphabets().
            sage: F.input_alphabet = [1]
            sage: Hd = G(F, algorithm='direct')
            sage: Hd.input_alphabet, Hd.output_alphabet
            ([1], ['a'])
            sage: He = G(F, algorithm='explorative')
            sage: He.input_alphabet, He.output_alphabet
            ([1], None)
            sage: G.output_alphabet = ['a']
            sage: Hd = G(F, algorithm='direct')
            sage: Hd.input_alphabet, Hd.output_alphabet
            ([1], ['a'])
            sage: He = G(F, algorithm='explorative')
            sage: He.input_alphabet, He.output_alphabet
            ([1], ['a'])
            sage: Hd == He
            True
            sage: F.input_alphabet = None
            sage: Hd = G(F, algorithm='direct')
            sage: Hd.input_alphabet, Hd.output_alphabet
            ([1], ['a'])
            sage: He = G(F, algorithm='explorative')
            Traceback (most recent call last):
            ...
            ValueError: No input alphabet is given. Try calling
            determine_alphabets().
        """
        if not other._allow_composition_:
            raise TypeError("Composition with automaton is not "
                            "possible.")

        if algorithm is None:
            if (any(len(t.word_out) > 1 for t in other.iter_transitions())
                or
                any(len(t.word_in) != 1 for t in self.iter_transitions())):
                algorithm = 'explorative'
            else:
                algorithm = 'direct'
        if algorithm == 'direct':
            return self._composition_direct_(other, only_accessible_components)
        elif algorithm == 'explorative':
            return self._composition_explorative_(other)
        else:
            raise ValueError("Unknown algorithm %s." % (algorithm,))

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

        result = other.product_FiniteStateMachine(
            self, function,
            only_accessible_components=only_accessible_components,
            final_function=lambda s1, s2: [],
            new_class=self.__class__)

        for state_result in result.iter_states():
            state = state_result.label()[0]
            if state.is_final:
                accept, _, output = self.process(
                    state.final_word_out,
                    initial_state=self.state(state_result.label()[1]))
                if not accept:
                    state_result.is_final = False
                else:
                    state_result.is_final = True
                    state_result.final_word_out = output

        return result

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

        Check that colors are correctly dealt with, cf. :issue:`19199`.
        In particular, the new colors have to be hashable such that
        :meth:`Automaton.determinisation` does not fail::

            sage: T = Transducer([[0, 0, 0, 0]], initial_states=[0])
            sage: A = T.input_projection()
            sage: B = A.composition(T, algorithm='explorative')
            sage: B.states()[0].color is None
            True
            sage: A.state(0).color = 0
            sage: B = A.composition(T, algorithm='explorative')
            sage: B.states()[0].color
            (None, 0)
            sage: B.determinisation()
            Automaton with 1 state
        """
        def composition_transition(states, input):
            (state1, state2) = states
            return [((new_state1, new_state2), output_second)
                    for _, new_state1, output_first in
                    first.process([input],
                                  list_of_outputs=True,
                                  initial_state=state1,
                                  write_final_word_out=False)
                    for _, new_state2, output_second in
                    second.process(output_first,
                                   list_of_outputs=True,
                                   initial_state=state2,
                                   write_final_word_out=False,
                                   always_include_output=True)]

        first = other
        if any(len(t.word_in) > 1
               for t in first.iter_transitions()):
            first = first.split_transitions()

        second = self
        if any(len(t.word_in) > 1
               for t in second.iter_transitions()):
            second = second.split_transitions()

        F = first.empty_copy(new_class=second.__class__)
        new_initial_states = itertools.product(
            first.iter_initial_states(),
            second.iter_initial_states())
        F.add_from_transition_function(composition_transition,
                                       initial_states=new_initial_states)

        for state in F.iter_states():
            (state1, state2) = state.label()
            if state1.is_final:
                final_output_second = second.process(
                    state1.final_word_out,
                    list_of_outputs=True,
                    initial_state=state2,
                    only_accepted=True,
                    always_include_output=True)
                if (len(final_output_second) > 1 and
                    not equal(r[2] for r in final_output_second)):
                    raise NotImplementedError("Stopping in state %s "
                                              "leads to "
                                              "non-deterministic final "
                                              "output." % state)
                if final_output_second:
                    state.is_final = True
                    state.final_word_out = final_output_second[0][2]

            if all(s.color is None for s in state.label()):
                state.color = None
            else:
                state.color = tuple(s.color for s in state.label())

        F.output_alphabet = second.output_alphabet
        return F

    def input_projection(self):
        """
        Return an automaton where the output of each transition of
        ``self`` is deleted.

        OUTPUT: an automaton

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
        Return a automaton where the input of each transition of self
        is deleted and the new input is the original output.

        OUTPUT: an automaton

        EXAMPLES::

            sage: F = FiniteStateMachine([('A', 'B', 0, 1), ('A', 'A', 1, 1),
            ....:                         ('B', 'B', 1, 0)])
            sage: G = F.output_projection()
            sage: G.transitions()
            [Transition from 'A' to 'B': 1|-,
             Transition from 'A' to 'A': 1|-,
             Transition from 'B' to 'B': 0|-]

        Final output words are also considered correctly::

            sage: H = Transducer([('A', 'B', 0, 1), ('A', 'A', 1, 1),
            ....:                 ('B', 'B', 1, 0), ('A', ('final', 0), 0, 0)],
            ....:                final_states=['A', 'B'])
            sage: H.state('B').final_word_out = 2
            sage: J = H.output_projection()
            sage: J.states()
            ['A', 'B', ('final', 0), ('final', 1)]
            sage: J.transitions()
            [Transition from 'A' to 'B': 1|-,
             Transition from 'A' to 'A': 1|-,
             Transition from 'A' to ('final', 0): 0|-,
             Transition from 'B' to 'B': 0|-,
             Transition from 'B' to ('final', 1): 2|-]
            sage: J.final_states()
            ['A', ('final', 1)]
        """
        return self.projection(what='output')

    def projection(self, what='input'):
        """
        Return an Automaton which transition labels are the projection
        of the transition labels of the input.

        INPUT:

        - ``what`` -- (default: ``input``) either ``input`` or ``output``

        OUTPUT: an automaton

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
        # TODO: use empty_copy() in order to
        # preserve on_duplicate_transition and future extensions.
        # for this, empty_copy would need a new optional argument
        # use_class=None ?

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

        if what == 'output':
            states = [s for s in self.iter_final_states() if s.final_word_out]
            if not states:
                return new
            number = 0
            while new.has_state(('final', number)):
                number += 1
            final = new.add_state(('final', number))
            final.is_final = True
            for state in states:
                output = state.final_word_out
                new.state(state_mapping[state]).final_word_out = []
                new.state(state_mapping[state]).is_final = False
                new.add_transition((state_mapping[state], final, output, None))

        return new

    def transposition(self, reverse_output_labels=True):
        """
        Return a new finite state machine, where all transitions of the
        input finite state machine are reversed.

        INPUT:

        - ``reverse_output_labels`` -- boolean (default: ``True``); whether to
          reverse output labels

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

        ::

            sage: A = Automaton([(0, 1, [1, 0])],
            ....:     initial_states=[0],
            ....:     final_states=[1])
            sage: A([1, 0])
            True
            sage: A.transposition()([0, 1])
            True

        ::

            sage: T = Transducer([(0, 1, [1, 0], [1, 0])],
            ....:     initial_states=[0],
            ....:     final_states=[1])
            sage: T([1, 0])
            [1, 0]
            sage: T.transposition()([0, 1])
            [0, 1]
            sage: T.transposition(reverse_output_labels=False)([0, 1])
            [1, 0]


        TESTS:

        If a final state of ``self`` has a non-empty final output word,
        transposition is not implemented::

            sage: T = Transducer([('1', '1', 1, 0), ('1', '2', 0, 1),
            ....:                 ('2', '2', 0, 2)],
            ....:                 initial_states=['1'],
            ....:                 final_states=['1', '2'])
            sage: T.state('1').final_word_out = [2, 5]
            sage: T.transposition()
            Traceback (most recent call last):
            ...
            NotImplementedError: Transposition for transducers with
            final output words is not implemented.
        """
        if reverse_output_labels:
            rewrite_output = lambda word: list(reversed(word))
        else:
            rewrite_output = lambda word: word

        transposition = self.empty_copy()

        for state in self.iter_states():
            transposition.add_state(deepcopy(state))

        for transition in self.iter_transitions():
            transposition.add_transition(
                transition.to_state.label(), transition.from_state.label(),
                list(reversed(transition.word_in)),
                rewrite_output(transition.word_out))

        for initial in self.iter_initial_states():
            state = transposition.state(initial.label())
            if not initial.is_final:
                state.is_final = True
                state.is_initial = False

        for final in self.iter_final_states():
            state = transposition.state(final.label())
            if final.final_word_out:
                raise NotImplementedError("Transposition for transducers "
                                          "with final output words is not "
                                          "implemented.")
            if not final.is_initial:
                state.is_final = False
                state.is_initial = True

        return transposition

    def split_transitions(self):
        """
        Return a new transducer, where all transitions in ``self`` with input
        labels consisting of more than one letter
        are replaced by a path of the corresponding length.

        OUTPUT: a new transducer

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

    def final_components(self):
        """
        Return the final components of a finite state machine as finite
        state machines.

        OUTPUT:

        A list of finite state machines, each representing a final
        component of ``self``.

        A final component of a transducer ``T`` is a strongly connected
        component ``C`` such that there are no transitions of ``T``
        leaving ``C``.

        The final components are the only parts of a transducer which
        influence the main terms of the asymptotic behaviour of the sum
        of output labels of a transducer, see [HKP2015]_ and [HKW2015]_.

        EXAMPLES::

            sage: T = Transducer([['A', 'B', 0, 0], ['B', 'C', 0, 1],
            ....:                 ['C', 'B', 0, 1], ['A', 'D', 1, 0],
            ....:                 ['D', 'D', 0, 0], ['D', 'B', 1, 0],
            ....:                 ['A', 'E', 2, 0], ['E', 'E', 0, 0]])
            sage: FC = T.final_components()
            sage: sorted(FC[0].transitions())
            [Transition from 'B' to 'C': 0|1,
             Transition from 'C' to 'B': 0|1]
            sage: FC[1].transitions()
            [Transition from 'E' to 'E': 0|0]

        Another example (cycle of length 2)::

            sage: T = Automaton([[0, 1, 0], [1, 0, 0]])
            sage: len(T.final_components()) == 1
            True
            sage: T.final_components()[0].transitions()
            [Transition from 0 to 1: 0|-,
             Transition from 1 to 0: 0|-]
        """
        DG = self.digraph()
        condensation = DG.strongly_connected_components_digraph()
        return [self.induced_sub_finite_state_machine([self.state(_) for _ in component])
                for component in condensation.vertices(sort=True)
                if condensation.out_degree(component) == 0]

    def completion(self, sink=None):
        """
        Return a completion of this finite state machine.

        INPUT:

        - ``sink`` -- either an instance of :class:`FSMState` or a label
          for the sink (default: ``None``); if ``None``, the least
          available nonzero integer is used

        OUTPUT:

        A :class:`FiniteStateMachine` of the same type as this finite
        state machine.

        The resulting finite state machine is a complete version of this
        finite state machine.  A finite state machine is considered to
        be complete if each transition has an input label of length one
        and for each pair `(q, a)` where `q` is a state and `a` is an
        element of the input alphabet, there is exactly one transition
        from `q` with input label `a`.

        If this finite state machine is already complete, a deep copy is
        returned. Otherwise, a new non-final state (usually called a
        sink) is created and transitions to this sink are introduced as
        appropriate.

        EXAMPLES::

            sage: F = FiniteStateMachine([(0, 0, 0, 0),
            ....:                         (0, 1, 1, 1),
            ....:                         (1, 1, 0, 0)])
            sage: F.is_complete()
            False
            sage: G1 = F.completion()
            sage: G1.is_complete()
            True
            sage: G1.transitions()
            [Transition from 0 to 0: 0|0,
             Transition from 0 to 1: 1|1,
             Transition from 1 to 1: 0|0,
             Transition from 1 to 2: 1|-,
             Transition from 2 to 2: 0|-,
             Transition from 2 to 2: 1|-]
            sage: G2 = F.completion('Sink')
            sage: G2.is_complete()
            True
            sage: G2.transitions()
            [Transition from 0 to 0: 0|0,
             Transition from 0 to 1: 1|1,
             Transition from 1 to 1: 0|0,
             Transition from 1 to 'Sink': 1|-,
             Transition from 'Sink' to 'Sink': 0|-,
             Transition from 'Sink' to 'Sink': 1|-]
            sage: F.completion(1)
            Traceback (most recent call last):
            ...
            ValueError: The finite state machine already contains a state
            '1'.

        An input alphabet must be given::

            sage: F = FiniteStateMachine([(0, 0, 0, 0),
            ....:                         (0, 1, 1, 1),
            ....:                         (1, 1, 0, 0)],
            ....:                        determine_alphabets=False)
            sage: F.is_complete()
            Traceback (most recent call last):
            ...
            ValueError: No input alphabet is given. Try calling
            determine_alphabets().

        Non-deterministic machines are not allowed. ::

            sage: F = FiniteStateMachine([(0, 0, 0, 0), (0, 1, 0, 0)])
            sage: F.is_complete()
            False
            sage: F.completion()
            Traceback (most recent call last):
            ...
            ValueError: The finite state machine must be deterministic.
            sage: F = FiniteStateMachine([(0, 0, [0, 0], 0)])
            sage: F.is_complete()
            False
            sage: F.completion()
            Traceback (most recent call last):
            ...
            ValueError: The finite state machine must be deterministic.

        .. SEEALSO::

            :meth:`is_complete`,
            :meth:`split_transitions`,
            :meth:`determine_alphabets`,
            :meth:`is_deterministic`.

        TESTS:

        Test the use of an :class:`FSMState` as sink::

            sage: F = FiniteStateMachine([(0, 0, 0, 0),
            ....:                         (0, 1, 1, 1),
            ....:                         (1, 1, 0, 0)])
            sage: from sage.combinat.finite_state_machine import FSMState
            sage: F.completion(FSMState(1))
            Traceback (most recent call last):
            ...
            ValueError: The finite state machine already contains a state
            '1'.
            sage: s = FSMState(2)
            sage: G = F.completion(s)
            sage: G.state(2) is s
            True
        """
        result = deepcopy(self)
        if result.is_complete():
            return result
        if not result.is_deterministic():
            raise ValueError(
                "The finite state machine must be deterministic.")

        if sink is not None:
            try:
                s = result.state(sink)
                raise ValueError("The finite state machine already "
                                 "contains a state '%s'." % s.label())
            except LookupError:
                pass
        else:
            sink = 1 + max(itertools.chain(
                    [-1],
                    (s.label() for s in result.iter_states()
                     if s.label() in ZZ)))

        sink_state = result.add_state(sink)

        for state in result.iter_states():
            for transition in state.transitions:
                if len(transition.word_in) != 1:
                    raise ValueError(
                        "Transitions with input labels of length greater "
                        "than one are not allowed. Try calling "
                        "split_transitions().")

            existing = set(transition.word_in[0]
                           for transition in state.transitions)
            for missing in set(result.input_alphabet) - existing:
                result.add_transition(state, sink_state, missing)

        return result

    # *************************************************************************
    # simplifications
    # *************************************************************************

    def prepone_output(self):
        """
        For all paths, shift the output of the path from one
        transition to the earliest possible preceding transition of
        the path.

        OUTPUT: nothing

        Apply the following to each state `s` (except initial states) of the
        finite state machine as often as possible:

        If the letter `a` is a prefix of the output label of all transitions from
        `s` (including the final output of `s`), then remove it from all these
        labels and append it to all output labels of all transitions leading
        to `s`.

        We assume that the states have no output labels, but final outputs are
        allowed.

        EXAMPLES::

            sage: A = Transducer([('A', 'B', 1, 1),
            ....:                 ('B', 'B', 0, 0),
            ....:                 ('B', 'C', 1, 0)],
            ....:                initial_states=['A'],
            ....:                final_states=['C'])
            sage: A.prepone_output()
            sage: A.transitions()
            [Transition from 'A' to 'B': 1|1,0,
             Transition from 'B' to 'B': 0|0,
             Transition from 'B' to 'C': 1|-]

        ::

            sage: B = Transducer([('A', 'B', 0, 1),
            ....:                 ('B', 'C', 1, [1, 1]),
            ....:                 ('B', 'C', 0, 1)],
            ....:                initial_states=['A'],
            ....:                final_states=['C'])
            sage: B.prepone_output()
            sage: B.transitions()
            [Transition from 'A' to 'B': 0|1,1,
             Transition from 'B' to 'C': 1|1,
             Transition from 'B' to 'C': 0|-]

        If initial states are not labeled as such, unexpected results may be
        obtained::

            sage: C = Transducer([(0,1,0,0)])
            sage: C.prepone_output()
            verbose 0 (...: finite_state_machine.py, prepone_output)
            All transitions leaving state 0 have an output label with
            prefix 0.  However, there is no inbound transition and it
            is not an initial state. This routine (possibly called by
            simplification) therefore erased this prefix from all
            outbound transitions.
            sage: C.transitions()
            [Transition from 0 to 1: 0|-]

        Also the final output of final states can be changed::

            sage: T = Transducer([('A', 'B', 0, 1),
            ....:                 ('B', 'C', 1, [1, 1]),
            ....:                 ('B', 'C', 0, 1)],
            ....:                initial_states=['A'],
            ....:                final_states=['B'])
            sage: T.state('B').final_word_out = [1]
            sage: T.prepone_output()
            sage: T.transitions()
            [Transition from 'A' to 'B': 0|1,1,
             Transition from 'B' to 'C': 1|1,
             Transition from 'B' to 'C': 0|-]
            sage: T.state('B').final_word_out
            []

        ::

            sage: S = Transducer([('A', 'B', 0, 1),
            ....:                 ('B', 'C', 1, [1, 1]),
            ....:                 ('B', 'C', 0, 1)],
            ....:                initial_states=['A'],
            ....:                final_states=['B'])
            sage: S.state('B').final_word_out = [0]
            sage: S.prepone_output()
            sage: S.transitions()
            [Transition from 'A' to 'B': 0|1,
             Transition from 'B' to 'C': 1|1,1,
             Transition from 'B' to 'C': 0|1]
            sage: S.state('B').final_word_out
            [0]

        Output labels do not have to be hashable::

            sage: C = Transducer([(0, 1, 0, []),
            ....:                 (1, 0, 0, [vector([0, 0]), 0]),
            ....:                 (1, 1, 1, [vector([0, 0]), 1]),
            ....:                 (0, 0, 1, 0)],
            ....:                 determine_alphabets=False,
            ....:                 initial_states=[0])
            sage: C.prepone_output()
            sage: sorted(C.transitions())
            [Transition from 0 to 1: 0|(0, 0),
             Transition from 0 to 0: 1|0,
             Transition from 1 to 0: 0|0,
             Transition from 1 to 1: 1|1,(0, 0)]
        """
        def find_common_output(state):
            if (any(transition for transition in self.transitions(state)
                    if not transition.word_out)
                        or state.is_final and not state.final_word_out):
                return tuple()
            first_letters = [transition.word_out[0]
                             for transition in self.transitions(state)]
            if state.is_final:
                first_letters = first_letters + [state.final_word_out[0]]
            if not first_letters:
                return tuple()
            first_item = first_letters.pop()
            if all(item == first_item for item in first_letters):
                return (first_item,)
            return tuple()

        changed = 1
        iteration = 0
        while changed > 0:
            changed = 0
            iteration += 1
            for state in self.iter_states():
                if state.is_initial:
                    continue
                if state.word_out:
                    raise NotImplementedError(
                        "prepone_output assumes that all states have "
                        "empty output word, but state %s has output "
                        "word %s" % (state, state.word_out))
                common_output = find_common_output(state)
                if common_output:
                    changed += 1
                    if state.is_final:
                        assert state.final_word_out[0] == common_output[0]
                        state.final_word_out = state.final_word_out[1:]
                    for transition in self.transitions(state):
                        assert transition.word_out[0] == common_output[0]
                        transition.word_out = transition.word_out[1:]
                    found_inbound_transition = False
                    for transition in self.iter_transitions():
                        if transition.to_state == state:
                            transition.word_out = transition.word_out \
                                + [common_output[0]]
                            found_inbound_transition = True
                    if not found_inbound_transition:
                        verbose(
                            "All transitions leaving state %s have an "
                            "output label with prefix %s. However, "
                            "there is no inbound transition and it is "
                            "not an initial state. This routine "
                            "(possibly called by simplification) "
                            "therefore erased this prefix from all "
                            "outbound transitions." %
                            (state, common_output[0]),
                            level=0)

    def equivalence_classes(self):
        r"""
        Return a list of equivalence classes of states.

        Two states `a` and `b` are equivalent if and only if there is
        a bijection `\varphi` between paths starting at `a` and paths
        starting at `b` with the following properties: Let `p_a` be a
        path from `a` to `a'` and `p_b` a path from `b` to `b'` such
        that `\varphi(p_a)=p_b`, then

        - `p_a.\mathit{word}_\mathit{in}=p_b.\mathit{word}_\mathit{in}`,
        - `p_a.\mathit{word}_\mathit{out}=p_b.\mathit{word}_\mathit{out}`,
        - `a'` and `b'` have the same output label, and
        - `a'` and `b'` are both final or both non-final and have the
          same final output word.

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
            sage: sorted(fsm.equivalence_classes())
            [['A', 'C'], ['B', 'D']]
            sage: fsm.state("A").is_final = True
            sage: sorted(fsm.equivalence_classes())
            [['A'], ['B'], ['C'], ['D']]
            sage: fsm.state("C").is_final = True
            sage: sorted(fsm.equivalence_classes())
            [['A', 'C'], ['B', 'D']]
            sage: fsm.state("A").final_word_out = 1
            sage: sorted(fsm.equivalence_classes())
            [['A'], ['B'], ['C'], ['D']]
            sage: fsm.state("C").final_word_out = 1
            sage: sorted(fsm.equivalence_classes())
            [['A', 'C'], ['B', 'D']]
        """

        # Two states `a` and `b` are j-equivalent if and only if there
        # is a bijection `\varphi` between paths of length <= j
        # starting at `a` and paths starting at `b` with the following
        # properties: Let `p_a` be a path from `a` to `a'` and `p_b` a
        # path from `b` to `b'` such that `\varphi(p_a)=p_b`, then
        #
        # - `p_a.\mathit{word}_{in}=p_b.\mathit{word}_{in}`,
        # - `p_a.\mathit{word}_{out}=p_b.\mathit{word}_{out}`,
        # - `a'` and `b'` have the same output label, and
        # - `a'` and `b'` are both final or both non-final.

        # If for some j the relations j-1 equivalent and j-equivalent
        # coincide, then they are equal to the equivalence relation
        # described in the docstring.

        # classes_current holds the equivalence classes of
        # j-equivalence, classes_previous holds the equivalence
        # classes of j-1 equivalence.

        # initialize with 0-equivalence
        classes_previous = []
        key_0 = lambda state: (state.is_final, state.color, state.word_out,
                               state.final_word_out)
        states_grouped = full_group_by(self.states(), key=key_0)
        classes_current = [equivalence_class for
                           (key, equivalence_class) in states_grouped]

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
                                       (key, equivalence_class) in states_grouped])

        return classes_current

    def quotient(self, classes):
        r"""
        Construct the quotient with respect to the equivalence classes.

        INPUT:

        - ``classes`` is a list of equivalence classes of states.

        OUTPUT: a finite state machine

        The labels of the new states are tuples of states of the
        ``self``, corresponding to ``classes``.

        Assume that `c` is a class, and `a` and `b` are states in
        `c`. Then there is a bijection `\varphi` between the
        transitions from `a` and the transitions from `b` with the
        following properties: if `\varphi(t_a)=t_b`, then

        - `t_a.\mathit{word}_\mathit{in}=t_b.\mathit{word}_\mathit{in}`,
        - `t_a.\mathit{word}_\mathit{out}=t_b.\mathit{word}_\mathit{out}`, and
        - `t_a` and `t_b` lead to some equivalent states `a'` and `b'`.

        Non-initial states may be merged with initial states, the
        resulting state is an initial state.

        All states in a class must have the same ``is_final``,
        ``final_word_out`` and ``word_out`` values.

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
            AssertionError: Transitions of state 'A' and 'B' are incompatible.

        TESTS::

            sage: fsm = FiniteStateMachine([("A", "B", 0, 1), ("A", "B", 1, 0),
            ....:                           ("B", "C", 0, 0), ("B", "C", 1, 1),
            ....:                           ("C", "D", 0, 1), ("C", "D", 1, 0),
            ....:                           ("D", "A", 0, 0), ("D", "A", 1, 1)],
            ....:                           final_states=["A", "C"])
            sage: fsm.state("A").final_word_out = 1
            sage: fsm.state("C").final_word_out = 2
            sage: fsmq = fsm.quotient([[fsm.state("A"), fsm.state("C")],
            ....:                      [fsm.state("B"), fsm.state("D")]])
            Traceback (most recent call last):
                ...
            AssertionError: Class ['A', 'C'] mixes
            final states with different final output words.
        """
        new = self.empty_copy()
        state_mapping = {}

        # Create new states and build state_mapping
        for c in classes:
            new_label = tuple(c)
            new_state = c[0].relabeled(new_label)
            new.add_state(new_state)
            for state in c:
                state_mapping[state] = new_state

        # Copy data from old transducer
        for c in classes:
            new_state = state_mapping[c[0]]
            sorted_transitions = sorted(
                [(state_mapping[t.to_state], t.word_in, t.word_out)
                 for t in c[0].transitions])
            for transition in self.iter_transitions(c[0]):
                new.add_transition(
                    from_state=new_state,
                    to_state=state_mapping[transition.to_state],
                    word_in=transition.word_in,
                    word_out=transition.word_out)

            # check that all class members have the same information (modulo classes)
            for state in c:
                new_state.is_initial = new_state.is_initial or state.is_initial
                assert new_state.is_final == state.is_final, \
                    "Class %s mixes final and non-final states" % (c,)
                assert new_state.word_out == state.word_out, \
                    "Class %s mixes different word_out" % (c,)
                assert new_state.color == state.color, \
                    "Class %s mixes different colors" % (c,)
                assert sorted_transitions == sorted(
                    [(state_mapping[t.to_state], t.word_in, t.word_out)
                     for t in state.transitions]), \
                    "Transitions of state %s and %s are incompatible." % (c[0], state)
                assert new_state.final_word_out == state.final_word_out, \
                    "Class %s mixes final states with different " \
                    "final output words." % (c,)
        return new

    def merged_transitions(self):
        """
        Merges transitions which have the same ``from_state``,
        ``to_state`` and ``word_out`` while adding their ``word_in``.

        OUTPUT:

        A finite state machine with merged transitions. If no mergers occur,
        return ``self``.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import duplicate_transition_add_input
            sage: T = Transducer([[1, 2, 1/4, 1], [1, -2, 1/4, 1], [1, -2, 1/2, 1],
            ....:                 [2, 2, 1/4, 1], [2, -2, 1/4, 1], [-2, -2, 1/4, 1],
            ....:                 [-2, 2, 1/4, 1], [2, 3, 1/2, 1], [-2, 3, 1/2, 1]],
            ....:                on_duplicate_transition=duplicate_transition_add_input)
            sage: T1 = T.merged_transitions()
            sage: T1 is T
            False
            sage: sorted(T1.transitions())
            [Transition from -2 to -2: 1/4|1,
             Transition from -2 to 2: 1/4|1,
             Transition from -2 to 3: 1/2|1,
             Transition from 1 to 2: 1/4|1,
             Transition from 1 to -2: 3/4|1,
             Transition from 2 to -2: 1/4|1,
             Transition from 2 to 2: 1/4|1,
             Transition from 2 to 3: 1/2|1]

        Applying the function again does not change the result::

            sage: T2 = T1.merged_transitions()
            sage: T2 is T1
            True
        """
        def key(transition):
            return (transition.to_state, transition.word_out)

        new = self.empty_copy()
        changed = False
        state_dict = {}
        memo = {}

        for state in self.states():
            new_state = deepcopy(state, memo)
            state_dict[state] = new_state
            new.add_state(new_state)

        for state in self.states():
            grouped_transitions = itertools.groupby(sorted(state.transitions, key=key), key=key)
            for (to_state, word_out), transitions in grouped_transitions:
                transition_list = list(transitions)
                changed = changed or len(transition_list) > 1
                word_in = 0
                for transition in transition_list:
                    if isinstance(transition.word_in, Iterable) and len(transition.word_in) == 1:
                        word_in += transition.word_in[0]
                    else:
                        raise TypeError('%s does not have a list of length 1 as word_in' % transition)
                new.add_transition((state, to_state, word_in, word_out))

        if changed:
            return new
        else:
            return self

    def markov_chain_simplification(self):
        """
        Consider ``self`` as Markov chain with probabilities as input labels
        and simplify it.

        OUTPUT: simplified version of ``self``

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import duplicate_transition_add_input
            sage: T = Transducer([[1, 2, 1/4, 0], [1, -2, 1/4, 0], [1, -2, 1/2, 0],
            ....:                 [2, 2, 1/4, 1], [2, -2, 1/4, 1], [-2, -2, 1/4, 1],
            ....:                 [-2, 2, 1/4, 1], [2, 3, 1/2, 2], [-2, 3, 1/2, 2]],
            ....:                initial_states=[1],
            ....:                final_states=[3],
            ....:                on_duplicate_transition=duplicate_transition_add_input)
            sage: T1 = T.markov_chain_simplification()
            sage: sorted(T1.transitions())
            [Transition from ((1,),) to ((2, -2),): 1|0,
             Transition from ((2, -2),) to ((2, -2),): 1/2|1,
             Transition from ((2, -2),) to ((3,),): 1/2|2]
        """
        current = self.merged_transitions()
        number_states = len(current.states())

        while True:
            current = current.simplification()
            new_number_states = len(current.states())
            new = current.merged_transitions()
            if new is current and number_states == new_number_states:
                return new
            current = new
            number_states = new_number_states

    def with_final_word_out(self, letters, allow_non_final=True):
        """
        Construct a new finite state machine with final output words
        for all states by implicitly reading trailing letters until a
        final state is reached.

        INPUT:

        - ``letters`` -- either an element of the input alphabet or a
          list of such elements. This is repeated cyclically when
          needed.

        - ``allow_non_final`` -- boolean (default: ``True``); whether we allow
          that some states may be non-final in the resulting finite state
          machine. I.e., if ``False`` then each state has to have a path to a
          final state with input label matching ``letters``.

        OUTPUT: a finite state machine

        The inplace version of this function is
        :meth:`.construct_final_word_out`.

        Suppose for the moment a single element ``letter`` as input
        for ``letters``. This is equivalent to ``letters = [letter]``.
        We will discuss the general case below.

        Let ``word_in`` be a word over the input alphabet and assume
        that the original finite state machine transforms ``word_in`` to
        ``word_out`` reaching a possibly non-final state ``s``. Let
        further `k` be the minimum number of letters ``letter`` such
        that there is a path from ``s`` to some final state ``f`` whose
        input label consists of `k` copies of ``letter`` and whose
        output label is ``path_word_out``. Then the state ``s`` of the
        resulting finite state machine is a final state with final
        output ``path_word_out + f.final_word_out``. Therefore, the new
        finite state machine transforms ``word_in`` to ``word_out +
        path_word_out + f.final_word_out``.

        This is e.g. useful for finite state machines operating on digit
        expansions: there, it is sometimes required to read a sufficient
        number of trailing zeros (at the most significant positions) in
        order to reach a final state and to flush all carries. In this
        case, this method constructs an essentially equivalent finite
        state machine in the sense that it not longer requires adding
        sufficiently many trailing zeros. However, it is the
        responsibility of the user to make sure that if adding trailing
        zeros to the input anyway, the output is equivalent.

        If ``letters`` consists of more than one letter, then it is
        assumed that (not necessarily complete) cycles of ``letters``
        are appended as trailing input.

        .. SEEALSO::

            :ref:`example on Gray code <finite_state_machine_gray_code_example>`

        EXAMPLES:

            #.  A simple transducer transforming `00` blocks to `01`
                blocks::

                    sage: T = Transducer([(0, 1, 0, 0), (1, 0, 0, 1)],
                    ....:                initial_states=[0],
                    ....:                final_states=[0])
                    sage: T.process([0, 0, 0])
                    (False, 1, [0, 1, 0])
                    sage: T.process([0, 0, 0, 0])
                    (True, 0, [0, 1, 0, 1])
                    sage: F = T.with_final_word_out(0)
                    sage: for f in F.iter_final_states():
                    ....:     print("{} {}".format(f, f.final_word_out))
                    0 []
                    1 [1]
                    sage: F.process([0, 0, 0])
                    (True, 1, [0, 1, 0, 1])
                    sage: F.process([0, 0, 0, 0])
                    (True, 0, [0, 1, 0, 1])

            #.  A more realistic example: Addition of `1` in binary. We
                construct a transition function transforming the input
                to its binary expansion::

                    sage: def binary_transition(carry, input):
                    ....:     value = carry + input
                    ....:     if value.mod(2) == 0:
                    ....:         return (value/2, 0)
                    ....:     else:
                    ....:         return ((value-1)/2, 1)

                Now, we only have to start with a carry of `1` to
                get the required transducer::

                    sage: T = Transducer(binary_transition,
                    ....:                input_alphabet=[0, 1],
                    ....:                initial_states=[1],
                    ....:                final_states=[0])

                We test this for the binary expansion of `7`::

                    sage: T.process([1, 1, 1])
                    (False, 1, [0, 0, 0])

                The final carry `1` has not be flushed yet, we have to add a
                trailing zero::

                    sage: T.process([1, 1, 1, 0])
                    (True, 0, [0, 0, 0, 1])

                We check that with this trailing zero, the transducer
                performs as advertised::

                    sage: all(ZZ(T(k.bits()+[0]), base=2) == k + 1
                    ....:     for k in srange(16))
                    True

                However, most of the time, we produce superfluous trailing
                zeros::

                    sage: T(11.bits()+[0])
                    [0, 0, 1, 1, 0]

                We now use this method::

                    sage: F = T.with_final_word_out(0)
                    sage: for f in F.iter_final_states():
                    ....:     print("{} {}".format(f, f.final_word_out))
                    1 [1]
                    0 []

                The same tests as above, but we do not have to pad with
                trailing zeros anymore::

                    sage: F.process([1, 1, 1])
                    (True, 1, [0, 0, 0, 1])
                    sage: all(ZZ(F(k.bits()), base=2) == k + 1
                    ....:     for k in srange(16))
                    True

                No more trailing zero in the output::

                    sage: F(11.bits())
                    [0, 0, 1, 1]
                    sage: all(F(k.bits())[-1] == 1
                    ....:     for k in srange(16))
                    True

            #.  Here is an example, where we allow trailing repeated `10`::

                    sage: T = Transducer([(0, 1, 0, 'a'),
                    ....:                 (1, 2, 1, 'b'),
                    ....:                 (2, 0, 0, 'c')],
                    ....:                initial_states=[0],
                    ....:                final_states=[0])
                    sage: F = T.with_final_word_out([1, 0])
                    sage: for f in F.iter_final_states():
                    ....:     print(str(f) + ' ' + ''.join(f.final_word_out))
                    0
                    1 bc

                Trying this with trailing repeated `01` does not produce
                a ``final_word_out`` for state ``1``, but for state ``2``::

                    sage: F = T.with_final_word_out([0, 1])
                    sage: for f in F.iter_final_states():
                    ....:     print(str(f) + ' ' + ''.join(f.final_word_out))
                    0
                    2 c

            #.  Here another example with a more-letter trailing input::

                    sage: T = Transducer([(0, 1, 0, 'a'),
                    ....:                 (1, 2, 0, 'b'), (1, 2, 1, 'b'),
                    ....:                 (2, 3, 0, 'c'), (2, 0, 1, 'e'),
                    ....:                 (3, 1, 0, 'd'), (3, 1, 1, 'd')],
                    ....:                initial_states=[0],
                    ....:                final_states=[0],
                    ....:                with_final_word_out=[0, 0, 1, 1])
                    sage: for f in T.iter_final_states():
                    ....:     print(str(f) + ' ' + ''.join(f.final_word_out))
                    0
                    1 bcdbcdbe
                    2 cdbe
                    3 dbe

        TESTS:

            #.  Reading copies of ``letter`` may result in a cycle. In
                this simple example, we have no final state at all::

                    sage: T = Transducer([(0, 1, 0, 0), (1, 0, 0, 0)],
                    ....:                initial_states=[0])
                    sage: T.with_final_word_out(0)
                    Traceback (most recent call last):
                    ...
                    ValueError: The finite state machine contains
                    a cycle starting at state 0 with input label 0
                    and no final state.

            #.  A unique transition with input word ``letter`` is
                required::

                    sage: T = Transducer([(0, 1, 0, 0), (0, 2, 0, 0)])
                    sage: T.with_final_word_out(0)
                    Traceback (most recent call last):
                    ...
                    ValueError: No unique transition leaving state 0
                    with input label 0.

                It is not a problem if there is no transition starting
                at state ``1`` with input word ``letter``::

                    sage: T = Transducer([(0, 1, 0, 0)])
                    sage: F = T.with_final_word_out(0)
                    sage: for f in F.iter_final_states():
                    ....:     print(f, f.final_word_out)

                Anyhow, you can override this by::

                    sage: T = Transducer([(0, 1, 0, 0)])
                    sage: T.with_final_word_out(0, allow_non_final=False)
                    Traceback (most recent call last):
                    ...
                    ValueError: No unique transition leaving state 1
                    with input label 0.

            #.  All transitions must have input labels of length `1`::

                    sage: T = Transducer([(0, 0, [], 0)])
                    sage: T.with_final_word_out(0)
                    Traceback (most recent call last):
                    ...
                    NotImplementedError: All transitions must have input
                    labels of length 1. Consider calling split_transitions().
                    sage: T = Transducer([(0, 0, [0, 1], 0)])
                    sage: T.with_final_word_out(0)
                    Traceback (most recent call last):
                    ...
                    NotImplementedError: All transitions must have input
                    labels of length 1. Consider calling split_transitions().

            #.  An empty list as input is not allowed::

                    sage: T = Transducer([(0, 0, [], 0)])
                    sage: T.with_final_word_out([])
                    Traceback (most recent call last):
                    ...
                    ValueError: letters is not allowed to be an empty list.
        """
        new = deepcopy(self)
        new.construct_final_word_out(letters, allow_non_final)
        return new

    def construct_final_word_out(self, letters, allow_non_final=True):
        """
        This is an inplace version of :meth:`.with_final_word_out`. See
        :meth:`.with_final_word_out` for documentation and examples.

        TESTS::

            sage: T = Transducer([(0, 1, 0, 0), (1, 0, 0, 1)],
            ....:                initial_states=[0],
            ....:                final_states=[0])
            sage: F = T.with_final_word_out(0)
            sage: T.construct_final_word_out(0)
            sage: T == F  # indirect doctest
            True
            sage: T = Transducer([(0, 1, 0, None)],
            ....:                final_states=[1])
            sage: F = T.with_final_word_out(0)
            sage: F.state(0).final_word_out
            []
        """

        if not isinstance(letters, list):
            letters = [letters]
        elif not letters:
            raise ValueError(
                "letters is not allowed to be an empty list.")

        in_progress = set()
        cache = {}

        def find_final_word_out(state):
            # The return value is the output which is produced when
            # reading the given letters until a final state is reached.
            # If no final state can be reached, then None is returned.
            # For final states, the final word out is returned.
            # For final states with empty final output, that is [].
            position, letter = next(trailing_letters)
            if state.is_final:
                return state.final_word_out

            if (state, position) in cache:
                return cache[state, position]

            if (state, position) in in_progress:
                raise ValueError(
                    "The finite state machine contains a cycle "
                    "starting at state %s with input label %s "
                    "and no final state." % (state, letter))

            if any(len(t.word_in) != 1 for t in state.transitions):
                raise NotImplementedError(
                    "All transitions must have input labels of length "
                    "1. Consider calling split_transitions().")

            transitions = [t for t in state.transitions
                           if t.word_in == [letter]]
            if allow_non_final and not transitions:
                final_word_out = None
            elif len(transitions) != 1:
                raise ValueError(
                    "No unique transition leaving state %s with input "
                    "label %s." % (state, letter))
            else:
                in_progress.add((state, position))
                next_word = find_final_word_out(transitions[0].to_state)
                if next_word is not None:
                    final_word_out = transitions[0].word_out + next_word
                else:
                    final_word_out = None
                in_progress.remove((state, position))

            cache[state, position] = final_word_out
            return final_word_out

        for state in self.iter_states():
            assert not in_progress
            # trailing_letters is an infinite iterator additionally
            # marking positions
            trailing_letters = itertools.cycle(enumerate(letters))
            find_final_word_out(state)

        # actual modifications can only be carried out after all final words
        # have been computed as it may not be permissible to stop at a
        # formerly non-final state unless a cycle has been completed.

        for (state, position), final_word_out in cache.items():
            if position == 0 and final_word_out is not None:
                state.is_final = True
                state.final_word_out = final_word_out

    # *************************************************************************
    # other
    # *************************************************************************

    def graph(self, edge_labels='words_in_out'):
        """
        Return the graph of the finite state machine with labeled
        vertices and labeled edges.

        INPUT:

        - ``edge_label`` -- (default: ``'words_in_out'``) can be
            - ``'words_in_out'`` (labels will be strings ``'i|o'``)
            - a function with which takes as input a transition
              and outputs (returns) the label

        OUTPUT: a :class:`directed graph <DiGraph>`

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = FSMState('A')
            sage: T = Transducer()
            sage: T.graph()
            Looped multi-digraph on 0 vertices
            sage: T.add_state(A)
            'A'
            sage: T.graph()
            Looped multi-digraph on 1 vertex
            sage: T.add_transition(('A', 'A', 0, 1))
            Transition from 'A' to 'A': 0|1
            sage: T.graph()
            Looped multi-digraph on 1 vertex

        .. SEEALSO::

            :class:`DiGraph`
        """
        if edge_labels == 'words_in_out':
            label_fct = lambda t: t._in_out_label_()
        elif callable(edge_labels):
            label_fct = edge_labels
        else:
            raise TypeError('Wrong argument for edge_labels.')

        graph_data = []
        isolated_vertices = []
        for state in self.iter_states():
            transitions = state.transitions
            if not transitions:
                isolated_vertices.append(state.label())
            graph_data.extend((t.from_state.label(), t.to_state.label(),
                               label_fct(t))
                              for t in transitions)

        G = DiGraph(graph_data, multiedges=True, loops=True)
        G.add_vertices(isolated_vertices)
        return G

    digraph = graph

    def plot(self):
        """
        Plots a graph of the finite state machine with labeled
        vertices and labeled edges.

        OUTPUT: a plot of the graph of the finite state machine

        TESTS::

            sage: FiniteStateMachine([('A', 'A', 0)]).plot()                            # needs sage.plot
            Graphics object consisting of 3 graphics primitives
        """
        return self.graph(edge_labels='words_in_out').plot()

    def predecessors(self, state, valid_input=None):
        """
        List all predecessors of a state.

        INPUT:

        - ``state`` -- the state from which the predecessors should be
          listed

        - ``valid_input`` -- if ``valid_input`` is a list, then we
          only consider transitions whose input labels are contained
          in ``valid_input``. ``state`` has to be a :class:`FSMState`
          (not a label of a state). If input labels of length larger
          than `1` are used, then ``valid_input`` has to be a list of
          lists.

        OUTPUT: list of states

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
        if valid_input is not None:
            valid_list = list()
            for input in valid_input:
                input_list = input
                if not isinstance(input_list, list):
                    input_list = [input]
                valid_list.append(input_list)
            valid_input = valid_list

        unhandeled_direct_predecessors = {s: [] for s in self.states()}
        for t in self.transitions():
            if valid_input is None or t.word_in in valid_input:
                unhandeled_direct_predecessors[t.to_state].append(t.from_state)
        done = []
        open = [state]
        while open:
            s = open.pop()
            candidates = unhandeled_direct_predecessors[s]
            if candidates is not None:
                open.extend(candidates)
                unhandeled_direct_predecessors[s] = None
                done.append(s)
        return done

    def number_of_words(self, variable=None,
                        base_ring=None):
        r"""
        Return the number of successful input words of given length.

        INPUT:

        - ``variable`` -- a symbol denoting the length of the words,
          by default `n`

        - ``base_ring`` -- ring (default: ``QQbar``) in which to
          compute the eigenvalues

        OUTPUT: a symbolic expression

        EXAMPLES::

            sage: NAFpm = Automaton([(0, 0, 0), (0, 1, 1),
            ....:                    (0, 1, -1), (1, 0, 0)],
            ....:                   initial_states=[0],
            ....:                   final_states=[0, 1])
            sage: N = NAFpm.number_of_words(); N                                        # needs sage.symbolic
            4/3*2^n - 1/3*(-1)^n
            sage: all(len(list(NAFpm.language(s)))                                      # needs sage.symbolic
            ....:     - len(list(NAFpm.language(s-1))) == N.subs(n=s)
            ....:     for s in srange(1, 6))
            True

        An example with non-rational eigenvalues. By default,
        eigenvalues are elements of the
        :mod:`field of algebraic numbers <sage.rings.qqbar>`. ::

            sage: NAFp = Automaton([(0, 0, 0), (0, 1, 1),  (1, 0, 0)],
            ....:                 initial_states=[0],
            ....:                 final_states=[0, 1])
            sage: N = NAFp.number_of_words(); N                                         # needs sage.rings.number_field sage.symbolic
            1.170820393249937?*1.618033988749895?^n
            - 0.1708203932499369?*(-0.618033988749895?)^n
            sage: all(len(list(NAFp.language(s)))                                       # needs sage.rings.number_field sage.symbolic
            ....:     - len(list(NAFp.language(s-1))) == N.subs(n=s)
            ....:     for s in srange(1, 6))
            True

        We specify a suitable ``base_ring`` to obtain a radical
        expression. To do so, we first compute the characteristic
        polynomial and then construct a number field generated by its
        roots. ::

            sage: # needs sage.rings.number_field sage.symbolic
            sage: M = NAFp.adjacency_matrix(entry=lambda t: 1)
            sage: M.characteristic_polynomial()
            x^2 - x - 1
            sage: R.<phi> = NumberField(x^2 - x - 1, embedding=1.6)
            sage: N = NAFp.number_of_words(base_ring=R); N
            1/2*(1/2*sqrt(5) + 1/2)^n*(3*sqrt(1/5) + 1)
            - 1/2*(-1/2*sqrt(5) + 1/2)^n*(3*sqrt(1/5) - 1)
            sage: all(len(list(NAFp.language(s)))
            ....:     - len(list(NAFp.language(s-1))) == N.subs(n=s)
            ....:     for s in srange(1, 6))
            True

        In this special case, we might also use the constant
        :class:`golden_ratio <sage.symbolic.constants.GoldenRatio>`::

            sage: # needs sage.rings.number_field sage.symbolic
            sage: R.<phi> = NumberField(x^2-x-1, embedding=golden_ratio)
            sage: N = NAFp.number_of_words(base_ring=R); N
            1/5*(3*golden_ratio + 1)*golden_ratio^n
            - 1/5*(3*golden_ratio - 4)*(-golden_ratio + 1)^n
            sage: all(len(list(NAFp.language(s)))
            ....:     - len(list(NAFp.language(s-1))) == N.subs(n=s)
            ....:     for s in srange(1, 6))
            True

        The adjacency matrix of the following example is a Jordan
        matrix of size 3 to the eigenvalue 4::

            sage: J3 = Automaton([(0, 1, -1), (1, 2, -1)],
            ....:     initial_states=[0],
            ....:     final_states=[0, 1, 2])
            sage: for i in range(3):
            ....:     for j in range(4):
            ....:         new_transition = J3.add_transition(i, i, j)
            sage: J3.adjacency_matrix(entry=lambda t: 1)
            [4 1 0]
            [0 4 1]
            [0 0 4]
            sage: N = J3.number_of_words(); N                                           # needs sage.symbolic
            1/2*4^(n - 2)*(n - 1)*n + 4^(n - 1)*n + 4^n
            sage: all(len(list(J3.language(s)))                                         # needs sage.symbolic
            ....:     - len(list(J3.language(s-1))) == N.subs(n=s)
            ....:     for s in range(1, 6))
            True

        Here is an automaton without cycles, so with eigenvalue `0`. ::

            sage: A = Automaton([(j, j+1, 0) for j in range(3)],
            ....:               initial_states=[0],
            ....:               final_states=list(range(3)))
            sage: A.number_of_words()                                                   # needs sage.symbolic
            1/2*0^(n - 2)*(n - 1)*n + 0^(n - 1)*n + 0^n

        TESTS::

            sage: A = Automaton([(0, 0, 0), (0, 1, 0)],
            ....:               initial_states=[0])
            sage: A.number_of_words()                                                   # needs sage.symbolic
            Traceback (most recent call last):
            ...
            NotImplementedError: Finite State Machine must be deterministic.
        """
        from sage.modules.free_module_element import vector
        from sage.arith.misc import binomial
        from sage.symbolic.ring import SR
        if base_ring is None:
            from sage.rings.qqbar import QQbar as base_ring
        if variable is None:
            variable = SR.symbol('n')

        def jordan_block_power(block, exponent):
            eigenvalue = SR(block[0, 0])
            return matrix(block.nrows(),
                          block.nrows(),
                          lambda i, j: eigenvalue**(exponent-(j-i)) *
                          binomial(exponent, j - i)
                          if j >= i else 0)

        if not self.is_deterministic():
            raise NotImplementedError("Finite State Machine must be deterministic.")

        left = vector(ZZ(s.is_initial) for s in self.iter_states())
        right = vector(ZZ(s.is_final) for s in self.iter_states())
        A = self.adjacency_matrix(entry=lambda t: 1)
        J, T = A.jordan_form(base_ring, transformation=True)
        Jpower = matrix.block_diagonal(
            [jordan_block_power(J.subdivision(j, j), variable)
             for j in range(len(J.subdivisions()[0]) + 1)])
        T_inv_right = T.solve_right(right).change_ring(SR)
        left_T = (left * T).change_ring(SR)
        return left_T * Jpower * T_inv_right

    def asymptotic_moments(self, variable=None):
        r"""
        Return the main terms of expectation and variance of the sum
        of output labels and its covariance with the sum of input
        labels.

        INPUT:

        - ``variable`` -- a symbol denoting the length of the input,
          (default: `n`)

        OUTPUT: a dictionary consisting of

        - ``expectation`` -- `e n + \operatorname{Order}(1)`,
        - ``variance`` -- `v n + \operatorname{Order}(1)`,
        - ``covariance`` -- `c n + \operatorname{Order}(1)`

        for suitable constants `e`, `v` and `c`.

        Assume that all input and output labels are numbers and that
        ``self`` is complete and has only one final component. Assume
        further that this final component is aperiodic. Furthermore,
        assume that there is exactly one initial state and that all
        states are final.

        Denote by `X_n` the sum of output labels written by the
        finite state machine when reading a random input word of
        length `n` over the input alphabet (assuming
        equidistribution).

        Then the expectation of `X_n` is `en+O(1)`, the variance
        of `X_n` is `vn+O(1)` and the covariance of `X_n` and
        the sum of input labels is `cn+O(1)`, cf. [HKW2015]_,
        Theorem 3.9.

        In the case of non-integer input or output labels, performance
        degrades significantly. For rational input and output labels,
        consider rescaling to integers. This limitation comes from the
        fact that determinants over polynomial rings can be computed
        much more efficiently than over the symbolic ring. In fact, we
        compute (parts) of a trivariate generating function where the
        input and output labels are exponents of some indeterminates,
        see [HKW2015]_, Theorem 3.9 for details. If those exponents are
        integers, we can use a polynomial ring.

        EXAMPLES:

        #.  A trivial example: write the negative of the input::

                sage: T = Transducer([(0, 0, 0, 0), (0, 0, 1, -1)],
                ....:                initial_states=[0],
                ....:                final_states=[0])
                sage: T([0, 1, 1])
                [0, -1, -1]

                sage: # needs sage.symbolic
                sage: moments = T.asymptotic_moments()
                sage: moments['expectation']
                -1/2*n + Order(1)
                sage: moments['variance']
                1/4*n + Order(1)
                sage: moments['covariance']
                -1/4*n + Order(1)

        #.  For the case of the Hamming weight of the non-adjacent-form
            (NAF) of integers, cf. the :wikipedia:`Non-adjacent_form`
            and the :ref:`example on recognizing NAFs
            <finite_state_machine_recognizing_NAFs_example>`, the
            following agrees with the results in [HP2007]_.

            We first use the transducer to convert the standard binary
            expansion to the NAF given in [HP2007]_. We use the parameter
            ``with_final_word_out`` such that we do not have to add
            sufficiently many trailing zeros::

                sage: NAF = Transducer([(0, 0, 0, 0),
                ....:                   (0, '.1', 1, None),
                ....:                   ('.1', 0, 0, [1, 0]),
                ....:                   ('.1', 1, 1, [-1, 0]),
                ....:                   (1, 1, 1, 0),
                ....:                   (1, '.1', 0, None)],
                ....:                  initial_states=[0],
                ....:                  final_states=[0],
                ....:                  with_final_word_out=[0])

            As an example, we compute the NAF of `27` by this
            transducer.

            ::

                sage: binary_27 = 27.bits()
                sage: binary_27
                [1, 1, 0, 1, 1]
                sage: NAF_27 = NAF(binary_27)
                sage: NAF_27
                [-1, 0, -1, 0, 0, 1, 0]
                sage: ZZ(NAF_27, base=2)
                27

            Next, we are only interested in the Hamming weight::

                sage: def weight(state, input):
                ....:     if input is None:
                ....:         result = 0
                ....:     else:
                ....:         result = ZZ(input != 0)
                ....:     return (0, result)
                sage: weight_transducer = Transducer(weight,
                ....:                                input_alphabet=[-1, 0, 1],
                ....:                                initial_states=[0],
                ....:                                final_states=[0])
                sage: NAFweight = weight_transducer.composition(NAF)
                sage: NAFweight.transitions()
                [Transition from (0, 0) to (0, 0): 0|0,
                 Transition from (0, 0) to ('.1', 0): 1|-,
                 Transition from ('.1', 0) to (0, 0): 0|1,0,
                 Transition from ('.1', 0) to (1, 0): 1|1,0,
                 Transition from (1, 0) to ('.1', 0): 0|-,
                 Transition from (1, 0) to (1, 0): 1|0]
                sage: NAFweight(binary_27)
                [1, 0, 1, 0, 0, 1, 0]

            Now, we actually compute the asymptotic moments::

                sage: # needs sage.symbolic
                sage: moments = NAFweight.asymptotic_moments()
                sage: moments['expectation']
                1/3*n + Order(1)
                sage: moments['variance']
                2/27*n + Order(1)
                sage: moments['covariance']
                Order(1)

        #.  This is Example 3.16 in [HKW2015]_, where a transducer with
            variable output labels is given. There, the aim was to
            choose the output labels of this very simple transducer such
            that the input and output sum are asymptotically
            independent, i.e., the constant `c` vanishes.

            ::

                sage: # needs sage.symbolic
                sage: var('a_1, a_2, a_3, a_4')
                (a_1, a_2, a_3, a_4)
                sage: T = Transducer([[0, 0, 0, a_1], [0, 1, 1, a_3],
                ....:                 [1, 0, 0, a_4], [1, 1, 1, a_2]],
                ....:                initial_states=[0], final_states=[0, 1])
                sage: moments = T.asymptotic_moments()
                verbose 0 (...) Non-integer output weights lead to
                significant performance degradation.
                sage: moments['expectation']
                1/4*(a_1 + a_2 + a_3 + a_4)*n + Order(1)
                sage: moments['covariance']
                -1/4*(a_1 - a_2)*n + Order(1)

            Therefore, the asymptotic covariance vanishes if and only if
            `a_2=a_1`.

        #.  This is Example 4.3 in [HKW2015]_, dealing with the
            transducer converting the binary expansion of an integer
            into Gray code (cf. the :wikipedia:`Gray_code` and the
            :ref:`example on Gray code
            <finite_state_machine_gray_code_example>`)::

                sage: # needs sage.symbolic
                sage: moments = transducers.GrayCode().asymptotic_moments()
                sage: moments['expectation']
                1/2*n + Order(1)
                sage: moments['variance']
                1/4*n + Order(1)
                sage: moments['covariance']
                Order(1)

        #.  This is the first part of Example 4.4 in [HKW2015]_,
            counting the number of 10 blocks in the standard binary
            expansion. The least significant digit is at the left-most
            position::

                sage: block10 = transducers.CountSubblockOccurrences(
                ....:     [1, 0],
                ....:     input_alphabet=[0, 1])
                sage: sorted(block10.transitions())
                [Transition from () to (): 0|0,
                 Transition from () to (1,): 1|0,
                 Transition from (1,) to (): 0|1,
                 Transition from (1,) to (1,): 1|0]

                sage: # needs sage.symbolic
                sage: moments = block10.asymptotic_moments()
                sage: moments['expectation']
                1/4*n + Order(1)
                sage: moments['variance']
                1/16*n + Order(1)
                sage: moments['covariance']
                Order(1)

        #.  This is the second part of Example 4.4 in [HKW2015]_,
            counting the number of 11 blocks in the standard binary
            expansion. The least significant digit is at the left-most
            position::

                sage: block11 = transducers.CountSubblockOccurrences(
                ....:     [1, 1],
                ....:     input_alphabet=[0, 1])
                sage: sorted(block11.transitions())
                [Transition from () to (): 0|0,
                 Transition from () to (1,): 1|0,
                 Transition from (1,) to (): 0|0,
                 Transition from (1,) to (1,): 1|1]

                sage: # needs sage.symbolic
                sage: var('N')
                N
                sage: moments = block11.asymptotic_moments(N)
                sage: moments['expectation']
                1/4*N + Order(1)
                sage: moments['variance']
                5/16*N + Order(1)
                sage: correlation = (moments['covariance'].coefficient(N) /
                ....:                (1/2 * sqrt(moments['variance'].coefficient(N))))
                sage: correlation
                2/5*sqrt(5)

        #.  This is Example 4.5 in [HKW2015]_, counting the number of
            01 blocks minus the number of 10 blocks in the standard binary
            expansion. The least significant digit is at the left-most
            position::

                sage: block01 = transducers.CountSubblockOccurrences(
                ....:     [0, 1],
                ....:     input_alphabet=[0, 1])
                sage: product_01x10 = block01.cartesian_product(block10)
                sage: block_difference = transducers.sub([0, 1])(product_01x10)
                sage: T = block_difference.simplification().relabeled()
                sage: T.transitions()
                [Transition from 0 to 2: 0|-1,
                 Transition from 0 to 0: 1|0,
                 Transition from 1 to 2: 0|0,
                 Transition from 1 to 0: 1|0,
                 Transition from 2 to 2: 0|0,
                 Transition from 2 to 0: 1|1]

                sage: # needs sage.symbolic
                sage: moments = T.asymptotic_moments()
                sage: moments['expectation']
                Order(1)
                sage: moments['variance']
                Order(1)
                sage: moments['covariance']
                Order(1)

        #.  The finite state machine must have a unique final component::

                sage: T = Transducer([(0, -1, -1, -1), (0, 1, 1, 1),
                ....:                 (-1, -1, -1, -1), (-1, -1, 1, -1),
                ....:                 (1, 1, -1, 1), (1, 1, 1, 1)],
                ....:                initial_states=[0],
                ....:                final_states=[0, 1, -1])
                sage: T.asymptotic_moments()
                Traceback (most recent call last):
                ...
                NotImplementedError: asymptotic_moments is only
                implemented for finite state machines with one final
                component.

            In this particular example, the first letter of the input
            decides whether we reach the loop at `-1` or the loop at
            `1`. In the first case, we have `X_n = -n`, while we have
            `X_n = n` in the second case. Therefore, the expectation
            `E(X_n)` of `X_n` is `E(X_n) = 0`. We get `(X_n-E(X_n))^2 =
            n^2` in all cases, which results in a variance of `n^2`.

            So this example shows that the variance may be non-linear if
            there is more than one final component.

        TESTS:

        #.  An input alphabet must be given::

                sage: T = Transducer([[0, 0, 0, 0]],
                ....:                initial_states=[0], final_states=[0],
                ....:                determine_alphabets=False)
                sage: T.asymptotic_moments()
                Traceback (most recent call last):
                ...
                ValueError: No input alphabet is given.
                Try calling determine_alphabets().

        #.  The finite state machine must have a unique initial state::

                sage: T = Transducer([(0, 0, 0, 0)])
                sage: T.asymptotic_moments()
                Traceback (most recent call last):
                ...
                ValueError: A unique initial state is required.

        #.  The finite state machine must be complete::

                sage: T = Transducer([[0, 0, 0, 0]],
                ....:                initial_states=[0], final_states=[0],
                ....:                input_alphabet=[0, 1])
                sage: T.asymptotic_moments()
                Traceback (most recent call last):
                ...
                NotImplementedError: This finite state machine is
                not complete.

        #.  The final component of the finite state machine must be
            aperiodic::

                sage: T = Transducer([(0, 1, 0, 0), (1, 0, 0, 0)],
                ....:                initial_states=[0], final_states=[0, 1])
                sage: T.asymptotic_moments()
                Traceback (most recent call last):
                ...
                NotImplementedError: asymptotic_moments is only
                implemented for finite state machines whose unique final
                component is aperiodic.

        #.  Non-integer input or output labels lead to a warning::

                sage: # needs sage.symbolic
                sage: T = Transducer([[0, 0, 0, 0], [0, 0, 1, -1/2]],
                ....:                initial_states=[0], final_states=[0])
                sage: moments = T.asymptotic_moments()
                verbose 0 (...) Non-integer output weights lead to
                significant performance degradation.
                sage: moments['expectation']
                -1/4*n + Order(1)
                sage: moments['variance']
                1/16*n + Order(1)
                sage: moments['covariance']
                -1/8*n + Order(1)

            This warning can be silenced by :func:`~sage.misc.verbose.set_verbose`::

                sage: # needs sage.symbolic
                sage: from sage.misc.verbose import set_verbose
                sage: set_verbose(-1, "finite_state_machine.py")
                sage: moments = T.asymptotic_moments()
                sage: moments['expectation']
                -1/4*n + Order(1)
                sage: moments['variance']
                1/16*n + Order(1)
                sage: moments['covariance']
                -1/8*n + Order(1)
                sage: set_verbose(0, "finite_state_machine.py")

        #.  Check whether ``word_out`` of ``FSMState`` are correctly
            dealt with::

                sage: from sage.combinat.finite_state_machine import FSMState
                sage: s = FSMState(0, word_out=2,
                ....:              is_initial=True,
                ....:              is_final=True)
                sage: T = Transducer([(s, s, 0, 1)],
                ....:                initial_states=[s], final_states=[s])
                sage: T([0, 0])
                [2, 1, 2, 1, 2]
                sage: T.asymptotic_moments()['expectation']                             # needs sage.symbolic
                3*n + Order(1)

            The same test for non-integer output::

                sage: from sage.combinat.finite_state_machine import FSMState
                sage: s = FSMState(0, word_out=2/3)
                sage: T = Transducer([(s, s, 0, 1/2)],
                ....:                initial_states=[s], final_states=[s])
                sage: T.asymptotic_moments()['expectation']                             # needs sage.symbolic
                verbose 0 (...) Non-integer output weights lead to
                significant performance degradation.
                7/6*n + Order(1)

        #.  All states of ``self`` have to be final::

                sage: T = Transducer([(0, 1, 1, 4)], initial_states=[0])
                sage: T.asymptotic_moments()
                Traceback (most recent call last):
                ...
                ValueError: Not all states are final.

        ALGORITHM:

        See [HKW2015]_, Theorem 3.9.

        REFERENCES:

        .. [HP2007] Clemens Heuberger and Helmut Prodinger, *The Hamming
           Weight of the Non-Adjacent-Form under Various Input Statistics*,
           Periodica Mathematica Hungarica Vol. 55 (1), 2007, pp. 81--96,
           :doi:`10.1007/s10998-007-3081-z`.
        """
        if self.input_alphabet is None:
            raise ValueError("No input alphabet is given. "
                             "Try calling determine_alphabets().")

        if len(self.initial_states()) != 1:
            raise ValueError("A unique initial state is required.")

        if not all(state.is_final for state in self.iter_states()):
            raise ValueError("Not all states are final.")

        if not self.is_complete():
            raise NotImplementedError("This finite state machine is "
                                      "not complete.")

        final_components = self.final_components()
        if len(final_components) != 1:
            raise NotImplementedError("asymptotic_moments is only "
                                      "implemented for finite state machines "
                                      "with one final component.")
        final_component = final_components[0]

        if not final_component.digraph().is_aperiodic():
            raise NotImplementedError("asymptotic_moments is only "
                                      "implemented for finite state machines "
                                      "whose unique final component is "
                                      "aperiodic.")

        from sage.calculus.functional import derivative
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.rings.rational_field import QQ
        from sage.symbolic.ring import SR

        if variable is None:
            variable = SR.symbol('n')

        def get_matrix(fsm, x, y):
            return fsm.adjacency_matrix(
                entry=lambda transition: x**sum(transition.word_in) *
                                         y**(sum(transition.word_out) +
                                         sum(transition.from_state.word_out)))

        K = len(self.input_alphabet)
        R = PolynomialRing(QQ, ("x", "y", "z"))
        (x, y, z) = R.gens()
        try:
            M = get_matrix(self, x, y)
        except (TypeError, ValueError):
            verbose(
                "Non-integer output weights lead to "
                "significant performance degradation.", level=0)
            # fall back to symbolic ring
            R = SR
            x = R.symbol()
            y = R.symbol()
            z = R.symbol()
            M = get_matrix(self, x, y)

            def substitute_one(g):
                return g.subs({x: 1, y: 1, z: 1})
        else:

            def substitute_one(g):
                # the result of the substitution shall live in QQ,
                # not in the polynomial ring R, so the method
                # subs does not achieve the result.
                # Therefore, we need this helper function.
                return g(1, 1, 1)

        f = (M.parent().identity_matrix() - z/K*M).det()
        f_x = substitute_one(derivative(f, x))
        f_y = substitute_one(derivative(f, y))
        f_z = substitute_one(derivative(f, z))
        f_xy = substitute_one(derivative(f, x, y))
        f_xz = substitute_one(derivative(f, x, z))
        f_yz = substitute_one(derivative(f, y, z))
        f_yy = substitute_one(derivative(f, y, y))
        f_zz = substitute_one(derivative(f, z, z))

        e_2 = f_y / f_z
        v_2 = (f_y**2 * (f_zz+f_z) + f_z**2 * (f_yy+f_y)
               - 2*f_y*f_z*f_yz) / f_z**3
        c = (f_x * f_y * (f_zz+f_z) + f_z**2 * f_xy - f_y*f_z*f_xz
             - f_x*f_z*f_yz) / f_z**3

        return {'expectation': e_2*variable + SR(1).Order(),
                'variance': v_2*variable + SR(1).Order(),
                'covariance': c*variable + SR(1).Order()}

    def moments_waiting_time(self, test=bool, is_zero=None,
                             expectation_only=False):
        r"""
        If this finite state machine acts as a Markov chain, return
        the expectation and variance of the number of steps until
        first writing ``True``.

        INPUT:

        - ``test`` -- (default: ``bool``) a callable deciding whether
          an output label is to be considered ``True``. By default, the
          standard conversion to boolean is used.

        - ``is_zero`` -- (default: ``None``) a callable deciding
          whether an expression for a probability is zero. By default,
          checking for zero is simply done by
          :meth:`~sage.structure.element.Element.is_zero`.  This
          parameter can be used to provide a more sophisticated check
          for zero, e.g. in the case of symbolic probabilities, see
          the examples below. This parameter is passed on to
          :meth:`is_Markov_chain`. This parameter only affects the
          input of the Markov chain.

        - ``expectation_only`` -- boolean (default: ``False``); if set, the
          variance is not computed (in order to save time). By default,
          the variance is computed.

        OUTPUT:

        A dictionary (if ``expectation_only=False``) consisting of

        - ``expectation``,
        - ``variance``.

        Otherwise, just the expectation is returned (no dictionary for
        ``expectation_only=True``).

        Expectation and variance of the number of steps until first
        writing ``True`` (as determined by the parameter ``test``).

        ALGORITHM:

        Relies on a (classical and easy) probabilistic argument,
        cf. [FGT1992]_, Eqns. (6) and (7).

        For the variance, see [FHP2015]_, Section 2.

        EXAMPLES:

        #.  The simplest example is to wait for the first `1` in a
            `0`-`1`-string where both digits appear with probability
            `1/2`. In fact, the waiting time equals `k` if and only if
            the string starts with `0^{k-1}1`. This event occurs with
            probability `2^{-k}`. Therefore, the expected waiting time
            and the variance are `\sum_{k\ge 1} k2^{-k}=2` and
            `\sum_{k\ge 1} (k-2)^2 2^{-k}=2`::

                sage: var('k')                                                          # needs sage.symbolic
                k
                sage: sum(k * 2^(-k), k, 1, infinity)                                   # needs sage.symbolic
                2
                sage: sum((k-2)^2 * 2^(-k), k, 1, infinity)                             # needs sage.symbolic
                2

            We now compute the same expectation and variance by using a
            Markov chain::

                sage: from sage.combinat.finite_state_machine import (
                ....:     duplicate_transition_add_input)
                sage: T = Transducer(
                ....:     [(0, 0, 1/2, 0), (0, 0, 1/2, 1)],
                ....:     on_duplicate_transition=\
                ....:         duplicate_transition_add_input,
                ....:     initial_states=[0],
                ....:     final_states=[0])
                sage: T.moments_waiting_time()
                {'expectation': 2, 'variance': 2}
                sage: T.moments_waiting_time(expectation_only=True)
                2

            In the following, we replace the output ``0`` by ``-1`` and
            demonstrate the use of the parameter ``test``::

                sage: T.delete_transition((0, 0, 1/2, 0))
                sage: T.add_transition((0, 0, 1/2, -1))
                Transition from 0 to 0: 1/2|-1
                sage: T.moments_waiting_time(test=lambda x: x<0)
                {'expectation': 2, 'variance': 2}

        #.  Make sure that the transducer is actually a Markov
            chain. Although this is checked by the code, unexpected
            behaviour may still occur if the transducer looks like a
            Markov chain. In the following example, we 'forget' to
            assign probabilities, but due to a coincidence, all
            'probabilities' add up to one. Nevertheless, `0` is never
            written, so the expectation is `1`.

            ::

                sage: T = Transducer([(0, 0, 0, 0), (0, 0, 1, 1)],
                ....:                on_duplicate_transition=\
                ....:                    duplicate_transition_add_input,
                ....:                initial_states=[0],
                ....:                final_states=[0])
                sage: T.moments_waiting_time()
                {'expectation': 1, 'variance': 0}

        #.  If ``True`` is never written, the moments are
            ``+Infinity``::

                sage: T = Transducer([(0, 0, 1, 0)],
                ....:                on_duplicate_transition=\
                ....:                    duplicate_transition_add_input,
                ....:                initial_states=[0],
                ....:                final_states=[0])
                sage: T.moments_waiting_time()
                {'expectation': +Infinity, 'variance': +Infinity}

        #.  Let `h` and `r` be positive integers. We consider random
            strings of letters `1`, `\ldots`, `r` where the letter `j`
            occurs with probability `p_j`. Let `B` be the random
            variable giving the first position of a block of `h`
            consecutive identical letters. Then

            .. MATH::

                \begin{aligned}
                \mathbb{E}(B)&=\frac1{\displaystyle\sum_{i=1}^r
                \frac1{p_i^{-1}+\cdots+p_i^{-h}}},\\
                \mathbb{V}(B)&=\frac{\displaystyle\sum_{i=1}^r\biggl(
                \frac{p_i +p_i^h}{1-p_i^h}
                - 2h\frac{ p_i^h(1-p_i)}{(1-p_i^h)^2}\biggr)}
                {\displaystyle\biggl(\sum_{i=1}^r
                \frac1{p_i^{-1}+\cdots+p_i^{-h}}\biggr)^2}
                \end{aligned}

            cf. [S1986]_, p. 62, or [FHP2015]_, Theorem 1. We now
            verify this with a transducer approach.

            ::

                sage: # needs sage.libs.singular
                sage: def test(h, r):
                ....:     R = PolynomialRing(
                ....:             QQ,
                ....:             names=['p_%d' % j for j in range(r)])
                ....:     p = R.gens()
                ....:     def is_zero(polynomial):
                ....:         return polynomial in (sum(p) - 1) * R
                ....:     theory_expectation = 1/(sum(1/sum(p[j]^(-i)
                ....:                     for i in range(1, h+1))
                ....:                     for j in range(r)))
                ....:     theory_variance = sum(
                ....:         (p[i] + p[i]^h)/(1 - p[i]^h)
                ....:         - 2*h*p[i]^h * (1 - p[i])/(1 - p[i]^h)^2
                ....:         for i in range(r)
                ....:         ) * theory_expectation^2
                ....:     alphabet = list(range(r))
                ....:     counters = [
                ....:         transducers.CountSubblockOccurrences([j]*h,
                ....:                     alphabet)
                ....:         for j in alphabet]
                ....:     all_counter = counters[0].cartesian_product(
                ....:         counters[1:])
                ....:     adder = transducers.add(input_alphabet=[0, 1],
                ....:         number_of_operands=r)
                ....:     probabilities = Transducer(
                ....:        [(0, 0, p[j], j) for j in alphabet],
                ....:        initial_states=[0],
                ....:        final_states=[0],
                ....:        on_duplicate_transition=\
                ....:            duplicate_transition_add_input)
                ....:     chain = adder(all_counter(probabilities))
                ....:     result = chain.moments_waiting_time(
                ....:        is_zero=is_zero)
                ....:     return is_zero((result['expectation'] -
                ....:                theory_expectation).numerator()) \
                ....:            and \
                ....:            is_zero((result['variance'] -
                ....:                 theory_variance).numerator())
                sage: test(2, 2)
                True
                sage: test(2, 3)
                True
                sage: test(3, 3)
                True

        #.  Consider the alphabet `\{0, \ldots, r-1\}`, some `1\le j\le
            r` and some `h\ge 1`.  For some probabilities `p_0`,
            `\ldots`, `p_{r-1}`, we consider infinite words where the
            letters occur independently with the given probabilities.
            The random variable `B_j` is the first position `n` such
            that there exist `j` of the `r` letters having an `h`-run.
            The expectation of `B_j` is given in [FHP2015]_, Theorem 2.
            Here, we verify this result by using transducers::

                sage: # needs sage.libs.singular
                sage: def test(h, r, j):
                ....:     R = PolynomialRing(
                ....:             QQ,
                ....:             names=['p_%d' % i for i in range(r)])
                ....:     p = R.gens()
                ....:     def is_zero(polynomial):
                ....:         return polynomial in (sum(p) - 1) * R
                ....:     alphabet = list(range(r))
                ....:     counters = [
                ....:         transducers.Wait([0, 1])(
                ....:             transducers.CountSubblockOccurrences(
                ....:                 [i]*h,
                ....:                 alphabet))
                ....:         for i in alphabet]
                ....:     all_counter = counters[0].cartesian_product(
                ....:         counters[1:])
                ....:     adder = transducers.add(input_alphabet=[0, 1],
                ....:         number_of_operands=r)
                ....:     threshold = transducers.map(
                ....:         f=lambda x: x >= j,
                ....:         input_alphabet=srange(r+1))
                ....:     probabilities = Transducer(
                ....:         [(0, 0, p[i], i) for i in alphabet],
                ....:         initial_states=[0],
                ....:         final_states=[0],
                ....:         on_duplicate_transition=\
                ....:             duplicate_transition_add_input)
                ....:     chain = threshold(adder(all_counter(
                ....:         probabilities)))
                ....:     result = chain.moments_waiting_time(
                ....:         is_zero=is_zero,
                ....:         expectation_only=True)
                ....:     R_v = PolynomialRing(
                ....:             QQ,
                ....:             names=['p_%d' % i for i in range(r)])
                ....:     v = R_v.gens()
                ....:     S = 1/(1 - sum(v[i]/(1+v[i])
                ....:                    for i in range(r)))
                ....:     alpha = [(p[i] - p[i]^h)/(1 - p[i])
                ....:              for i in range(r)]
                ....:     gamma = [p[i]/(1 - p[i]) for i in range(r)]
                ....:     alphabet_set = set(alphabet)
                ....:     expectation = 0
                ....:     for q in range(j):
                ....:         for M in Subsets(alphabet_set, q):
                ....:             summand = S
                ....:             for i in M:
                ....:                 summand = summand.subs(
                ....:                     {v[i]: gamma[i]}) -\
                ....:                     summand.subs({v[i]: alpha[i]})
                ....:             for i in alphabet_set - set(M):
                ....:                 summand = summand.subs(
                ....:                     {v[i]: alpha[i]})
                ....:             expectation += summand
                ....:     return is_zero((result - expectation).\
                ....:             numerator())
                sage: test(2, 3, 2)
                True

        REFERENCES:

        .. [FGT1992] Philippe Flajolet, Danièle Gardy, Loÿs Thimonier,
           *Birthday paradox, coupon collectors, caching algorithms and
           self-organizing search*, Discrete Appl. Math. 39 (1992),
           207--229, :doi:`10.1016/0166-218X(92)90177-C`.

        .. [FHP2015] Uta Freiberg, Clemens Heuberger, Helmut Prodinger,
           *Application of Smirnov Words to Waiting Time Distributions
           of Runs*, :arxiv:`1503.08096`.

        .. [S1986] Gábor J. Székely, *Paradoxes in Probability Theory
           and Mathematical Statistics*, D. Reidel Publishing Company.

        TESTS:

        Only Markov chains are acceptable::

            sage: T = transducers.Identity([0, 1, 2])
            sage: T.moments_waiting_time()
            Traceback (most recent call last):
            ...
            ValueError: Only Markov chains can compute
            moments_waiting_time.

        There must be a unique initial state::

            sage: T = Transducer([(0, 1, 1, 1), (1, 0, 1, 0)],
            ....:                on_duplicate_transition=\
            ....:                    duplicate_transition_add_input)
            sage: T.moments_waiting_time()
            Traceback (most recent call last):
            ...
            ValueError: Unique initial state is required.

        Using `0` as initial state in this example, a `1` is written in
        the first step with probability `1`, so the waiting time is
        always `1`::

            sage: T.state(0).is_initial = True
            sage: T.moments_waiting_time()
            {'expectation': 1, 'variance': 0}

        Using both `0` and `1` as initial states again yields an error
        message::

            sage: T.state(1).is_initial = True
            sage: T.moments_waiting_time()
            Traceback (most recent call last):
            ...
            ValueError: Unique initial state is required.

        Detection of infinite waiting time for symbolic probabilities::

            sage: R.<p, q> = PolynomialRing(QQ)
            sage: T = Transducer([(0, 0, p, 0), (0, 0, q, 0)],
            ....:                initial_states=[0],
            ....:                on_duplicate_transition=\
            ....:                    duplicate_transition_add_input)
            sage: T.moments_waiting_time(
            ....:     is_zero=lambda e: e in (p + q - 1)*R)
            {'expectation': +Infinity, 'variance': +Infinity}
        """
        from sage.modules.free_module_element import vector
        from sage.matrix.constructor import identity_matrix
        from sage.rings.polynomial.polynomial_ring_constructor import\
            PolynomialRing

        def default_is_zero(expression):
            return expression.is_zero()

        is_zero_function = default_is_zero
        if is_zero is not None:
            is_zero_function = is_zero

        if not self.is_Markov_chain(is_zero):
            raise ValueError("Only Markov chains can compute "
                             "moments_waiting_time.")

        if len(self.initial_states()) != 1:
            raise ValueError("Unique initial state is required.")

        def entry(transition):
            word_out = transition.word_out
            if not word_out or (len(word_out) == 1 and not test(word_out[0])):
                return transition.word_in[0]
            else:
                return 0

        relabeled = self.relabeled()
        n = len(relabeled.states())
        assert [s.label() for s in relabeled.states()] == list(range(n))
        entry_vector = vector(ZZ(s.is_initial)
                              for s in relabeled.states())
        exit_vector = vector([1] * n)
        transition_matrix = relabeled.adjacency_matrix(entry=entry)
        # transition_matrix is the probability transition matrix
        # of the part of the transducer before the occurrence of true
        # output.
        # We cannot use the input parameter of adjacency_matrix
        # because we want to check for "true" input in the sense
        # of python's boolean conversion. So we cannot give
        # input=[False] as this might lead to strange phenomena.
        if all(map(is_zero_function,
                   transition_matrix * exit_vector - exit_vector)):
            import sage.rings.infinity
            expectation = sage.rings.infinity.PlusInfinity()
            variance = sage.rings.infinity.PlusInfinity()
        else:
            if expectation_only:
                system_matrix = identity_matrix(n) - transition_matrix
                expectation = entry_vector * \
                    system_matrix.solve_right(exit_vector)
            else:
                base_ring = transition_matrix.parent().base_ring()
                from sage.rings.polynomial.multi_polynomial_ring \
                    import MPolynomialRing_base
                if isinstance(base_ring, MPolynomialRing_base):
                    # if base_ring is already a multivariate polynomial
                    # ring, extend it instead of creating a univariate
                    # polynomial ring over a polynomial ring.  This
                    # should improve performance.
                    R = PolynomialRing(base_ring.base_ring(),
                                       base_ring.variable_names()
                                       + ('Z_waiting_time',))
                else:
                    R = PolynomialRing(base_ring, 'Z_waiting_time')
                Z = R.gens()[-1]
                system_matrix = identity_matrix(n) - Z * transition_matrix
                G = entry_vector * system_matrix.solve_right(exit_vector)
                expectation = G.subs({Z: 1})
                variance = 2 * G.derivative(Z).subs({Z: 1}) \
                    + expectation \
                    - expectation**2

        if expectation_only:
            return expectation
        else:
            return {'expectation': expectation,
                    'variance': variance}

    def is_monochromatic(self):
        """
        Check whether the colors of all states are equal.

        OUTPUT: boolean

        EXAMPLES::

            sage: G = transducers.GrayCode()
            sage: [s.color for s in G.iter_states()]
            [None, None, None]
            sage: G.is_monochromatic()
            True
            sage: G.state(1).color = 'blue'
            sage: G.is_monochromatic()
            False
        """
        return equal(s.color for s in self.iter_states())

    def language(self, max_length=None, **kwargs):
        r"""
        Return all words that can be written by this transducer.

        INPUT:

        - ``max_length`` -- integer or ``None`` (default). Only
          output words which come from inputs of length at most
          ``max_length`` will be considered. If ``None``, then this
          iterates over all possible words without length restrictions.

        - ``kwargs`` -- will be passed on to the :class:`process
          iterator <FSMProcessIterator>`. See :meth:`process` for a
          description.

        OUTPUT: an iterator

        EXAMPLES::

            sage: NAF = Transducer([('I', 0, 0, None), ('I', 1, 1, None),
            ....:                   (0, 0, 0, 0), (0, 1, 1, 0),
            ....:                   (1, 0, 0, 1), (1, 2, 1, -1),
            ....:                   (2, 1, 0, 0), (2, 2, 1, 0)],
            ....:                  initial_states=['I'], final_states=[0],
            ....:                  input_alphabet=[0, 1])
            sage: sorted(NAF.language(4),
            ....:        key=lambda o: (ZZ(o, base=2), len(o)))
            [[], [0], [0, 0], [0, 0, 0],
             [1], [1, 0], [1, 0, 0],
             [0, 1], [0, 1, 0],
             [-1, 0, 1],
             [0, 0, 1],
             [1, 0, 1]]

        ::

            sage: iterator = NAF.language()
            sage: next(iterator)
            []
            sage: next(iterator)
            [0]
            sage: next(iterator)
            [1]
            sage: next(iterator)
            [0, 0]
            sage: next(iterator)
            [0, 1]

        .. SEEALSO::

            :meth:`Automaton.language`,
            :meth:`process`.

        TESTS::

            sage: T = Transducer([(0, 1, 0, 'a'), (1, 2, 1, 'b')],
            ....:                initial_states=[0], final_states=[0, 1, 2])
            sage: T.determine_alphabets()
            sage: list(T.language(2))
            [[], ['a'], ['a', 'b']]
            sage: list(T.language(3))
            [[], ['a'], ['a', 'b']]
            sage: from sage.combinat.finite_state_machine import  _FSMProcessIteratorAll_
            sage: it = T.iter_process(
            ....:     process_iterator_class=_FSMProcessIteratorAll_,
            ....:     max_length=3,
            ....:     process_all_prefixes_of_input=True)
            sage: for current in it:
            ....:     print(current)
            ....:     print("finished: {}".format([branch.output for branch in it._finished_]))
            process (1 branch)
            + at state 1
            +-- tape at 1, [['a']]
            finished: [[]]
            process (1 branch)
            + at state 2
            +-- tape at 2, [['a', 'b']]
            finished: [[], ['a']]
            process (0 branches)
            finished: [[], ['a'], ['a', 'b']]
        """
        kwargs['process_iterator_class'] = _FSMProcessIteratorAll_
        kwargs['max_length'] = max_length
        kwargs['process_all_prefixes_of_input'] = True
        it = self.iter_process(**kwargs)
        for _ in it:
            for branch in it._finished_:
                if branch.accept:
                    yield branch.output
            it._finished_ = []


# ****************************************************************************


def is_Automaton(FSM):
    """
    Test whether or not ``FSM`` inherits from :class:`Automaton`.

    TESTS::

        sage: from sage.combinat.finite_state_machine import is_FiniteStateMachine, is_Automaton
        sage: is_Automaton(FiniteStateMachine())
        doctest:warning...
        DeprecationWarning: The function is_Automaton is deprecated; use 'isinstance(..., Automaton)' instead.
        See https://github.com/sagemath/sage/issues/38032 for details.
        False
        sage: is_Automaton(Automaton())
        True
        sage: is_FiniteStateMachine(Automaton())
        doctest:warning...
        DeprecationWarning: The function is_FiniteStateMachine is deprecated; use 'isinstance(..., FiniteStateMachine)' instead.
        See https://github.com/sagemath/sage/issues/38032 for details.
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(38032, "The function is_Automaton is deprecated; use 'isinstance(..., Automaton)' instead.")
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
        sage: A([0])
        True
        sage: A([1, 1, 0])
        True
        sage: A([1, 0, 1])
        False

    Note that the full output of the commands can be obtained by
    calling :meth:`.process` and looks like this::

        sage: A.process([1, 0, 1])
        (False, 'P')

    TESTS::

        sage: Automaton()
        Empty automaton
    """

    def __init__(self, *args, **kwargs):
        """
        Initialize an automaton. See :class:`Automaton` and its parent
        :class:`FiniteStateMachine` for more information.

        TESTS::

            sage: Transducer()._allow_composition_
            True
            sage: Automaton()._allow_composition_
            False
        """
        super().__init__(*args, **kwargs)
        self._allow_composition_ = False

    def _repr_(self):
        """
        Represent the finite state machine as "Automaton with n
        states" where n is the number of states.

        OUTPUT: string

        EXAMPLES::

            sage: Automaton()._repr_()
            'Empty automaton'


        TESTS::

            sage: A = Automaton()
            sage: A
            Empty automaton
            sage: A.add_state(42)
            42
            sage: A
            Automaton with 1 state
            sage: A.add_state(43)
            43
            sage: A
            Automaton with 2 states
        """
        if not self._states_:
            return "Empty automaton"
        if len(self._states_) == 1:
            return "Automaton with 1 state"
        else:
            return "Automaton with %s states" % len(self._states_)

    def _latex_transition_label_(self, transition,
                                 format_function=None):
        r"""
        Return the proper transition label.

        INPUT:

        - ``transition`` -- a transition

        - ``format_function`` -- a function formatting the labels

        OUTPUT: string

        EXAMPLES::

            sage: F = Automaton([('A', 'B', 1)])
            sage: print(latex(F))  # indirect doctest
            \begin{tikzpicture}[auto, initial text=, >=latex]
            \node[state] (v0) at (3.000000, 0.000000) {$\text{\texttt{A}}$};
            \node[state] (v1) at (-3.000000, 0.000000) {$\text{\texttt{B}}$};
            \path[->] (v0) edge node[rotate=360.00, anchor=south] {$1$} (v1);
            \end{tikzpicture}

        TESTS::

            sage: F = Automaton([('A', 'B', 0, 1)])
            sage: t = F.transitions()[0]
            sage: F._latex_transition_label_(t)
            \left[0\right]
        """
        if format_function is None:
            format_function = latex
        return format_function(transition.word_in)

    def intersection(self, other, only_accessible_components=True):
        """
        Return a new automaton which accepts an input if it is
        accepted by both given automata.

        INPUT:

        - ``other`` -- an automaton

        - ``only_accessible_components`` -- if ``True`` (default), then
          the result is piped through :meth:`.accessible_components`. If no
          ``new_input_alphabet`` is given, it is determined by
          :meth:`.determine_alphabets`.

        OUTPUT:

        A new automaton which computes the intersection
        (see below) of the languages of ``self`` and ``other``.

        The set of states of the new automaton is the Cartesian product of the
        set of states of both given automata. There is a transition `((A, B),
        (C, D), a)` in the new automaton if there are transitions `(A, C, a)`
        and `(B, D, a)` in the old automata.

        The methods :meth:`.intersection` and
        :meth:`.cartesian_product` are the same (for automata).

        EXAMPLES::

            sage: aut1 = Automaton([('1', '2', 1),
            ....:                   ('2', '2', 1),
            ....:                   ('2', '2', 0)],
            ....:                  initial_states=['1'],
            ....:                  final_states=['2'],
            ....:                  determine_alphabets=True)
            sage: aut2 = Automaton([('A', 'A', 1),
            ....:                   ('A', 'B', 0),
            ....:                   ('B', 'B', 0),
            ....:                   ('B', 'A', 1)],
            ....:                  initial_states=['A'],
            ....:                  final_states=['B'],
            ....:                  determine_alphabets=True)
            sage: res = aut1.intersection(aut2)
            sage: (aut1([1, 1]), aut2([1, 1]), res([1, 1]))
            (True, False, False)
            sage: (aut1([1, 0]), aut2([1, 0]), res([1, 0]))
            (True, True, True)
            sage: res.transitions()
            [Transition from ('1', 'A') to ('2', 'A'): 1|-,
             Transition from ('2', 'A') to ('2', 'B'): 0|-,
             Transition from ('2', 'A') to ('2', 'A'): 1|-,
             Transition from ('2', 'B') to ('2', 'B'): 0|-,
             Transition from ('2', 'B') to ('2', 'A'): 1|-]

        For automata with epsilon-transitions, intersection is not well
        defined. But for any finite state machine, epsilon-transitions can be
        removed by :meth:`.remove_epsilon_transitions`.

        ::

            sage: a1 = Automaton([(0, 0, 0),
            ....:                 (0, 1, None),
            ....:                 (1, 1, 1),
            ....:                 (1, 2, 1)],
            ....:                 initial_states=[0],
            ....:                 final_states=[1],
            ....:                 determine_alphabets=True)
            sage: a2 = Automaton([(0, 0, 0), (0, 1, 1), (1, 1, 1)],
            ....:                 initial_states=[0],
            ....:                 final_states=[1],
            ....:                 determine_alphabets=True)
            sage: a1.intersection(a2)
            Traceback (most recent call last):
            ...
            ValueError: An epsilon-transition (with empty input)
            was found.
            sage: a1.remove_epsilon_transitions()  # not tested (since not implemented yet)
            sage: a1.intersection(a2)  # not tested
        """
        if not isinstance(other, Automaton):
            raise TypeError(
                 "Only an automaton can be intersected with an automaton.")

        def function(transition1, transition2):
            if not transition1.word_in or not transition2.word_in:
                raise ValueError(
                    "An epsilon-transition (with empty input) was found.")
            if transition1.word_in == transition2.word_in:
                return (transition1.word_in, None)
            else:
                raise LookupError

        return self.product_FiniteStateMachine(
            other,
            function,
            only_accessible_components=only_accessible_components)

    cartesian_product = intersection

    def determinisation(self):
        """
        Return a deterministic automaton which accepts the same input
        words as the original one.

        OUTPUT: a new automaton, which is deterministic

        The labels of the states of the new automaton are frozensets
        of states of ``self``. The color of a new state is the
        frozenset of colors of the constituent states of ``self``.
        Therefore, the colors of the constituent states have to be
        hashable. However, if all constituent states have color
        ``None``, then the resulting color is ``None``, too.

        The input alphabet must be specified.

        EXAMPLES::

            sage: aut = Automaton([('A', 'A', 0), ('A', 'B', 1), ('B', 'B', 1)],
            ....:                 initial_states=['A'], final_states=['B'])
            sage: aut.determinisation().transitions()
            [Transition from frozenset({'A'}) to frozenset({'A'}): 0|-,
             Transition from frozenset({'A'}) to frozenset({'B'}): 1|-,
             Transition from frozenset({'B'}) to frozenset(): 0|-,
             Transition from frozenset({'B'}) to frozenset({'B'}): 1|-,
             Transition from frozenset() to frozenset(): 0|-,
             Transition from frozenset() to frozenset(): 1|-]

        ::

            sage: A = Automaton([('A', 'A', 1), ('A', 'A', 0), ('A', 'B', 1),
            ....:                ('B', 'C', 0), ('C', 'C', 1), ('C', 'C', 0)],
            ....:               initial_states=['A'], final_states=['C'])
            sage: A.determinisation().states()
            [frozenset({'A'}),
             frozenset({'A', 'B'}),
             frozenset({'A', 'C'}),
             frozenset({'A', 'B', 'C'})]

        ::

            sage: A = Automaton([(0, 1, 1), (0, 2, [1, 1]), (0, 3, [1, 1, 1]),
            ....:                (1, 0, -1), (2, 0, -2), (3, 0, -3)],
            ....:               initial_states=[0], final_states=[0, 1, 2, 3])
            sage: B = A.determinisation().relabeled().coaccessible_components()
            sage: sorted(B.transitions())
            [Transition from 0 to 1: 1|-,
             Transition from 1 to 0: -1|-,
             Transition from 1 to 3: 1|-,
             Transition from 3 to 0: -2|-,
             Transition from 3 to 4: 1|-,
             Transition from 4 to 0: -3|-]

        Note that colors of states have to be hashable::

            sage: A = Automaton([[0, 0, 0]], initial_states=[0])
            sage: A.state(0).color = []
            sage: A.determinisation()
            Traceback (most recent call last):
            ...
            TypeError: unhashable type: 'list'
            sage: A.state(0).color = ()
            sage: A.determinisation()
            Automaton with 1 state

        If the colors of all constituent states are ``None``,
        the resulting color is ``None``, too (:issue:`19199`)::

            sage: A = Automaton([(0, 0, 0)],
            ....:               initial_states=[0],
            ....:               final_states=[0])
            sage: [s.color for s in A.determinisation().iter_states()]
            [None]

        TESTS:

        This is from :issue:`15078`, comment 13.

        ::

            sage: D = {'A': [('A', 'a'), ('B', 'a'), ('A', 'b')],
            ....:      'C': [], 'B': [('C', 'b')]}
            sage: auto = Automaton(D, initial_states=['A'], final_states=['C'])
            sage: auto.is_deterministic()
            False
            sage: auto.process(list('aaab'))
            [(False, 'A'), (True, 'C')]
            sage: auto.states()
            ['A', 'C', 'B']
            sage: Ddet = auto.determinisation()
            sage: Ddet
            Automaton with 3 states
            sage: Ddet.is_deterministic()
            True
            sage: sorted(Ddet.transitions())
            [Transition from frozenset({'A'}) to frozenset({'A', 'B'}): 'a'|-,
             Transition from frozenset({'A'}) to frozenset({'A'}): 'b'|-,
             Transition from frozenset({'A', 'B'}) to frozenset({'A', 'B'}): 'a'|-,
             Transition from frozenset({'A', 'B'}) to frozenset({'A', 'C'}): 'b'|-,
             Transition from frozenset({'A', 'C'}) to frozenset({'A', 'B'}): 'a'|-,
             Transition from frozenset({'A', 'C'}) to frozenset({'A'}): 'b'|-]
            sage: Ddet.initial_states()
            [frozenset({'A'})]
            sage: Ddet.final_states()
            [frozenset({'A', 'C'})]
            sage: Ddet.process(list('aaab'))
            (True, frozenset({'A', 'C'}))

        Test that :issue:`18992` is fixed::

            sage: A = Automaton([(0, 1, []), (1, 1, 0)],
            ....:               initial_states=[0], final_states=[1])
            sage: B = A.determinisation()
            sage: B.initial_states()
            [frozenset({0, 1})]
            sage: B.final_states()
            [frozenset({0, 1}), frozenset({1})]
            sage: B.transitions()
            [Transition from frozenset({0, 1}) to frozenset({1}): 0|-,
             Transition from frozenset({1}) to frozenset({1}): 0|-]
            sage: C = B.minimization().relabeled()
            sage: C.initial_states()
            [0]
            sage: C.final_states()
            [0]
            sage: C.transitions()
            [Transition from 0 to 0: 0|-]
        """
        if any(len(t.word_in) > 1 for t in self.iter_transitions()):
            return self.split_transitions().determinisation()

        epsilon_successors = {}
        direct_epsilon_successors = {}
        for state in self.iter_states():
            direct_epsilon_successors[state] = set(
                t.to_state
                for t in self.iter_transitions(state)
                if not t.word_in)
            epsilon_successors[state] = set([state])

        old_count_epsilon_successors = 0
        count_epsilon_successors = len(epsilon_successors)

        while old_count_epsilon_successors < count_epsilon_successors:
            old_count_epsilon_successors = count_epsilon_successors
            count_epsilon_successors = 0
            for state in self.iter_states():
                for direct_successor in direct_epsilon_successors[state]:
                    epsilon_successors[state] = epsilon_successors[state].union(epsilon_successors[direct_successor])
                count_epsilon_successors += len(epsilon_successors[state])

        def set_transition(states, letter):
            result = set()
            for state in states:
                for transition in self.iter_transitions(state):
                    if transition.word_in == [letter]:
                        result.add(transition.to_state)
            result = result.union(*(epsilon_successors[s] for s in result))
            return (frozenset(result), [])

        result = self.empty_copy()
        new_initial_states = [frozenset(set().union(
                    *(epsilon_successors[s]
                      for s in self.iter_initial_states()
                      )))]
        result.add_from_transition_function(set_transition,
                                            initial_states=new_initial_states)

        for state in result.iter_states():
            state.is_final = any(s.is_final for s in state.label())
            if all(s.color is None for s in state.label()):
                state.color = None
            else:
                state.color = frozenset(s.color for s in state.label())

        return result

    def minimization(self, algorithm=None):
        """
        Return the minimization of the input automaton as a new automaton.

        INPUT:

        - ``algorithm`` -- either Moore's algorithm (by
          ``algorithm='Moore'`` or as default for deterministic
          automata) or Brzozowski's algorithm (when
          ``algorithm='Brzozowski'`` or when the automaton is not
          deterministic) is used.

        OUTPUT: a new automaton

        The resulting automaton is deterministic and has a minimal
        number of states.

        EXAMPLES::

            sage: A = Automaton([('A', 'A', 1), ('A', 'A', 0), ('A', 'B', 1),
            ....:                ('B', 'C', 0), ('C', 'C', 1), ('C', 'C', 0)],
            ....:               initial_states=['A'], final_states=['C'])
            sage: B = A.minimization(algorithm='Brzozowski')
            sage: B_trans = B.transitions(B.states()[1])
            sage: B_trans # random
            [Transition from frozenset({frozenset({'B', 'C'}),
                                        frozenset({'A', 'C'}),
                                        frozenset({'A', 'B', 'C'})})
                        to frozenset({frozenset({'C'}),
                                      frozenset({'B', 'C'}),
                                      frozenset({'A', 'C'}),
                                      frozenset({'A', 'B', 'C'})}):
                        0|-,
             Transition from frozenset({frozenset({'B', 'C'}),
                                        frozenset({'A', 'C'}),
                                        frozenset({'A', 'B', 'C'})})
                        to frozenset({frozenset({'B', 'C'}),
                                      frozenset({'A', 'C'}),
                                      frozenset({'A', 'B', 'C'})}):
                        1|-]
            sage: len(B.states())
            3
            sage: C = A.minimization(algorithm='Brzozowski')
            sage: C_trans = C.transitions(C.states()[1])
            sage: B_trans == C_trans
            True
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
        deterministic = self.is_deterministic()

        if algorithm == "Moore" or (algorithm is None and deterministic):
            return self._minimization_Moore_()
        elif algorithm == "Brzozowski" or (algorithm is None and not deterministic):
            return self._minimization_Brzozowski_()
        else:
            raise NotImplementedError("Algorithm '%s' is not implemented. Choose 'Moore' or 'Brzozowski'" % algorithm)

    def _minimization_Brzozowski_(self):
        """
        Return a minimized automaton by using Brzozowski's algorithm.

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
        Return a minimized automaton by using Moore's algorithm.

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
        if self.is_deterministic():
            return self.quotient(self.equivalence_classes())
        else:
            raise NotImplementedError("Minimization via Moore's Algorithm is only "
                                      "implemented for deterministic finite state machines")

    def complement(self):
        r"""
        Return the complement of this automaton.

        OUTPUT: an :class:`Automaton`

        If this automaton recognizes language `\mathcal{L}` over an
        input alphabet `\mathcal{A}`, then the complement recognizes
        `\mathcal{A}\setminus\mathcal{L}`.

        EXAMPLES::

            sage: A = automata.Word([0, 1])
            sage: [w for w in ([], [0], [1], [0, 0], [0, 1], [1, 0], [1, 1])
            ....:  if A(w)]
            [[0, 1]]
            sage: Ac = A.complement()
            sage: Ac.transitions()
            [Transition from 0 to 1: 0|-,
             Transition from 0 to 3: 1|-,
             Transition from 2 to 3: 0|-,
             Transition from 2 to 3: 1|-,
             Transition from 1 to 2: 1|-,
             Transition from 1 to 3: 0|-,
             Transition from 3 to 3: 0|-,
             Transition from 3 to 3: 1|-]
            sage: [w for w in ([], [0], [1], [0, 0], [0, 1], [1, 0], [1, 1])
            ....:  if Ac(w)]
            [[], [0], [1], [0, 0], [1, 0], [1, 1]]

        The automaton must be deterministic::

            sage: A = automata.Word([0]) * automata.Word([1])
            sage: A.complement()
            Traceback (most recent call last):
            ...
            ValueError: The finite state machine must be deterministic.
            sage: Ac = A.determinisation().complement()
            sage: [w for w in ([], [0], [1], [0, 0], [0, 1], [1, 0], [1, 1])
            ....:  if Ac(w)]
            [[], [0], [1], [0, 0], [1, 0], [1, 1]]
        """
        result = self.completion()
        for state in result.iter_states():
            state.is_final = not state.is_final

        return result

    def is_equivalent(self, other):
        """
        Test whether two automata are equivalent, i.e., accept the same
        language.

        INPUT:

        - ``other`` -- an :class:`Automaton`

        EXAMPLES::

            sage: A = Automaton([(0, 0, 0), (0, 1, 1), (1, 0, 1)],
            ....:               initial_states=[0],
            ....:               final_states=[0])
            sage: B = Automaton([('a', 'a', 0), ('a', 'b', 1), ('b', 'a', 1)],
            ....:               initial_states=['a'],
            ....:               final_states=['a'])
            sage: A.is_equivalent(B)
            True
            sage: B.add_transition('b', 'a', 0)
            Transition from 'b' to 'a': 0|-
            sage: A.is_equivalent(B)
            False
        """
        A = self.minimization().relabeled()
        [initial] = A.initial_states()
        address = {initial: ()}
        for v in A.digraph().breadth_first_search(initial.label()):
            state = A.state(v)
            state_address = address[state]
            for t in A.iter_transitions(state):
                if t.to_state not in address:
                    address[t.to_state] = state_address + tuple(t.word_in)

        B = other.minimization().relabeled()
        labels = {B.process(path)[1].label(): state.label()
                  for (state, path) in address.items()}
        try:
            return A == B.relabeled(labels=labels)
        except KeyError:
            return False

    def process(self, *args, **kwargs):
        """
        Return whether the automaton accepts the input and the state
        where the computation stops.

        INPUT:

        - ``input_tape`` -- the input tape can be a list or an
          iterable with entries from the input alphabet. If we are
          working with a multi-tape machine (see parameter
          ``use_multitape_input`` and notes below), then the tape is a
          list or tuple of tracks, each of which can be a list or an
          iterable with entries from the input alphabet.

        - ``initial_state`` or ``initial_states`` -- the initial
          state(s) in which the machine starts. Either specify a
          single one with ``initial_state`` or a list of them with
          ``initial_states``. If both are given, ``initial_state``
          will be appended to ``initial_states``. If neither is
          specified, the initial states of the finite state machine
          are taken.

        - ``list_of_outputs`` -- (default: ``None``) a boolean or
          ``None``. If ``True``, then the outputs are given in list form
          (even if we have no or only one single output). If
          ``False``, then the result is never a list (an exception is
          raised if the result cannot be returned). If
          ``list_of_outputs=None`` the method determines automatically
          what to do (e.g. if a non-deterministic machine returns more
          than one path, then the output is returned in list form).

        - ``only_accepted`` -- boolean (default: ``False``); if set,
          then the first argument in the output is guaranteed to be
          ``True`` (if the output is a list, then the first argument
          of each element will be ``True``).

        - ``full_output`` -- boolean (default: ``True``); if set,
          then the full output is given, otherwise only whether the
          sequence is accepted or not (the first entry below only).

        - ``always_include_output`` -- if set (not by default), always
          return a triple containing the (non-existing) output. This
          is in order to obtain output compatible with that of
          :meth:`FiniteStateMachine.process`. If this parameter is set,
          ``full_output`` has no effect.

        - ``format_output`` -- a function that translates the written
          output (which is in form of a list) to something more
          readable. By default (``None``) identity is used here.

        - ``check_epsilon_transitions`` -- boolean (default: ``True``); a
          boolean. If ``False``, then epsilon transitions are not taken into
          consideration during process.

        - ``write_final_word_out`` -- boolean (default: ``True``); whether the
          final output words should be written or not

        - ``use_multitape_input`` -- boolean (default: ``False``); if ``True``,
          then the multi-tape mode of the process iterator is activated. See
          also the notes below for multi-tape machines.

        - ``process_all_prefixes_of_input`` -- boolean (default: ``False``); if
          ``True``, then each prefix of the input word is processed (instead of
          processing the whole input word at once). Consequently, there is an
          output generated for each of these prefixes.

        - ``process_iterator_class`` -- (default: ``None``) a class
          inherited from :class:`FSMProcessIterator`. If ``None``,
          then :class:`FSMProcessIterator` is taken. An instance of this
          class is created and is used during the processing.

        OUTPUT:

        The full output is a pair (or a list of pairs,
        cf. parameter ``list_of_outputs``), where

        - the first entry is ``True`` if the input string is accepted and

        - the second gives the state reached after processing the
          input tape (This is a state with label ``None`` if the input
          could not be processed, i.e., if at one point no
          transition to go on could be found.).

        If ``full_output`` is ``False``, then only the first entry
        is returned.

        If ``always_include_output`` is set, an additional third entry
        ``[]`` is included.

        Note that in the case the automaton is not
        deterministic, all possible paths are taken into account.
        You can use :meth:`.determinisation` to get a deterministic
        automaton machine.

        This function uses an iterator which, in its simplest form, goes
        from one state to another in each step. To decide which way to
        go, it uses the input words of the outgoing transitions and
        compares them to the input tape. More precisely, in each step,
        the iterator takes an outgoing transition of the current state,
        whose input label equals the input letter of the tape.

        If the choice of the outgoing transition is not unique (i.e.,
        we have a non-deterministic finite state machine), all
        possibilities are followed. This is done by splitting the
        process into several branches, one for each of the possible
        outgoing transitions.

        The process (iteration) stops if all branches are finished,
        i.e., for no branch, there is any transition whose input word
        coincides with the processed input tape. This can simply
        happen when the entire tape was read.

        Also see :meth:`~FiniteStateMachine.__call__` for a
        version of :meth:`.process` with shortened output.

        Internally this function creates and works with an instance of
        :class:`FSMProcessIterator`. This iterator can also be obtained
        with :meth:`~FiniteStateMachine.iter_process`.

        If working with multi-tape finite state machines, all input
        words of transitions are words of `k`-tuples of letters.
        Moreover, the input tape has to consist of `k` tracks, i.e.,
        be a list or tuple of `k` iterators, one for each track.

        .. WARNING::

            Working with multi-tape finite state machines is still
            experimental and can lead to wrong outputs.

        EXAMPLES:

        In the following examples, we construct an automaton which
        accepts non-adjacent forms (see also the example on
        :ref:`non-adjacent forms <finite_state_machine_recognizing_NAFs_example>`
        in the documentation of the module
        :doc:`finite_state_machine`)
        and then test it by feeding it with several binary digit
        expansions.

        ::

            sage: NAF = Automaton(
            ....:     {'_': [('_', 0), ('1', 1)], '1': [('_', 0)]},
            ....:     initial_states=['_'], final_states=['_', '1'])
            sage: [NAF.process(w) for w in [[0], [0, 1], [1, 1], [0, 1, 0, 1],
            ....:                           [0, 1, 1, 1, 0], [1, 0, 0, 1, 1]]]
            [(True, '_'), (True, '1'), (False, None),
             (True, '1'), (False, None), (False, None)]

        If we just want a condensed output, we use::

            sage: [NAF.process(w, full_output=False)
            ....:     for w in [[0], [0, 1], [1, 1], [0, 1, 0, 1],
            ....:               [0, 1, 1, 1, 0], [1, 0, 0, 1, 1]]]
            [True, True, False, True, False, False]

        It is equivalent to::

            sage: [NAF(w) for w in [[0], [0, 1], [1, 1], [0, 1, 0, 1],
            ....:                   [0, 1, 1, 1, 0], [1, 0, 0, 1, 1]]]
            [True, True, False, True, False, False]

        The following example illustrates the difference between
        non-existing paths and reaching a non-final state::

            sage: NAF.process([2])
            (False, None)
            sage: NAF.add_transition(('_', 's', 2))
            Transition from '_' to 's': 2|-
            sage: NAF.process([2])
            (False, 's')

        A simple example of a (non-deterministic) multi-tape automaton is the
        following: It checks whether the two input tapes have the same number
        of ones::

            sage: M = Automaton([('=', '=', (1, 1)),
            ....:                ('=', '=', (None, 0)),
            ....:                ('=', '=', (0, None)),
            ....:                ('=', '<', (None, 1)),
            ....:                ('<', '<', (None, 1)),
            ....:                ('<', '<', (None, 0)),
            ....:                ('=', '>', (1, None)),
            ....:                ('>', '>', (1, None)),
            ....:                ('>', '>', (0, None))],
            ....:               initial_states=['='],
            ....:               final_states=['='])
            sage: M.process(([1, 0, 1], [1, 0]), use_multitape_input=True)
            (False, '>')
            sage: M.process(([0, 1, 0], [0, 1, 1]), use_multitape_input=True)
            (False, '<')
            sage: M.process(([1, 1, 0, 1], [0, 0, 1, 0, 1, 1]),
            ....:           use_multitape_input=True)
            (True, '=')

        Alternatively, we can use the following (non-deterministic)
        multi-tape automaton for the same check::

            sage: N = Automaton([('=', '=', (0, 0)),
            ....:                ('=', '<', (None, 1)),
            ....:                ('<', '<', (0, None)),
            ....:                ('<', '=', (1, None)),
            ....:                ('=', '>', (1, None)),
            ....:                ('>', '>', (None, 0)),
            ....:                ('>', '=', (None, 1))],
            ....:               initial_states=['='],
            ....:               final_states=['='])
            sage: N.process(([1, 0, 1], [1, 0]), use_multitape_input=True)
            (False, '>')
            sage: N.process(([0, 1, 0], [0, 1, 1]), use_multitape_input=True)
            (False, '<')
            sage: N.process(([1, 1, 0, 1], [0, 0, 1, 0, 1, 1]),
            ....:           use_multitape_input=True)
            (True, '=')

        .. SEEALSO::

            :meth:`FiniteStateMachine.process`,
            :meth:`Transducer.process`,
            :meth:`~FiniteStateMachine.iter_process`,
            :meth:`~FiniteStateMachine.__call__`,
            :class:`FSMProcessIterator`.
        """
        # set default values
        options = copy(self._process_default_options_)
        options.update(kwargs)

        condensed_output = (options['list_of_outputs'] is False and
                            not options['full_output'])

        if condensed_output:
            options['list_of_outputs'] = True
            options['only_accepted'] = True

        result = super().process(*args, **options)

        if condensed_output:
            return any(result)
        return result

    def _process_convert_output_(self, output_data, **kwargs):
        """
        Helper function which converts the output of
        :meth:`FiniteStateMachine.process` to one suitable for
        automata.

        INPUT:

        - ``output_data`` -- a triple

        - ``full_output`` -- boolean

        - ``always_include_output`` -- if set (not by default), always
          return a triple containing the (non-existing) output. This
          is for compatibility with transducers.

        OUTPUT: the converted output

        TESTS::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: A = Automaton()
            sage: A._process_convert_output_((True, FSMState('a'), [1, 0, 1]),
            ....:                            full_output=False,
            ....:                            always_include_output=False)
            True
            sage: A._process_convert_output_((True, FSMState('a'), [1, 0, 1]),
            ....:                            full_output=True,
            ....:                            always_include_output=False)
            (True, 'a')
            sage: A._process_convert_output_((True, FSMState('a'), [1, 0, 1]),
            ....:                            full_output=False,
            ....:                            always_include_output=True)
            (True, 'a', [1, 0, 1])
        """
        if kwargs['always_include_output']:
            return super()._process_convert_output_(
                output_data, **kwargs)
        accept_input, current_state, _ = output_data
        if kwargs['full_output']:
            return (accept_input, current_state)
        else:
            return accept_input

    def shannon_parry_markov_chain(self):
        """
        Compute a time homogeneous Markov chain such that all words of a
        given length recognized by the original automaton occur as the
        output with the same weight; the transition probabilities
        correspond to the Parry measure.

        OUTPUT:

        A Markov chain. Its input labels are the transition probabilities, the
        output labels the labels of the original automaton. In order to obtain
        equal weight for all words of the same length, an "exit weight" is
        needed. It is stored in the attribute ``color`` of the states of the
        Markov chain. The weights of the words of the same length sum up to one
        up to an exponentially small error.

        The stationary distribution of this Markov chain is
        saved as the initial probabilities of the states.

        The transition probabilities correspond to the Parry measure
        (see [S1948]_ and [P1964]_).

        The automaton is assumed to be deterministic, irreducible and
        aperiodic. All states must be final.

        EXAMPLES::

            sage: NAF = Automaton([(0, 0, 0), (0, 1, 1), (0, 1, -1),
            ....:                  (1, 0, 0)], initial_states=[0],
            ....:                 final_states=[0, 1])
            sage: P_NAF = NAF.shannon_parry_markov_chain()                              # needs sage.symbolic
            sage: P_NAF.transitions()                                                   # needs sage.symbolic
            [Transition from 0 to 0: 1/2|0,
             Transition from 0 to 1: 1/4|1,
             Transition from 0 to 1: 1/4|-1,
             Transition from 1 to 0: 1|0]
            sage: for s in P_NAF.iter_states():                                         # needs sage.symbolic
            ....:     print(s.color)
            3/4
            3/2

        The stationary distribution is also computed and saved as the
        initial probabilities of the returned Markov chain::

            sage: for s in P_NAF.states():                                              # needs sage.symbolic
            ....:     print("{} {}".format(s, s.initial_probability))
            0 2/3
            1 1/3

        The automaton is assumed to be deterministic, irreducible and aperiodic::

            sage: A = Automaton([(0, 0, 0), (0, 1, 1), (1, 1, 1), (1, 1, 0)],
            ....:               initial_states=[0])
            sage: A.shannon_parry_markov_chain()
            Traceback (most recent call last):
            ...
            NotImplementedError: Automaton must be strongly connected.
            sage: A = Automaton([(0, 0, 0), (0, 1, 0)],
            ....:               initial_states=[0])
            sage: A.shannon_parry_markov_chain()
            Traceback (most recent call last):
            ...
            NotImplementedError: Automaton must be deterministic.
            sage: A = Automaton([(0, 1, 0), (1, 0, 0)],
            ....:               initial_states=[0])
            sage: A.shannon_parry_markov_chain()
            Traceback (most recent call last):
            ...
            NotImplementedError: Automaton must be aperiodic.

        All states must be final::

            sage: A = Automaton([(0, 1, 0), (0, 0, 1), (1, 0, 0)],
            ....:               initial_states=[0])
            sage: A.shannon_parry_markov_chain()
            Traceback (most recent call last):
            ...
            NotImplementedError: All states must be final.

        ALGORITHM:

        See [HKP2015a]_, Lemma 4.1.

        REFERENCES:

        .. [HKP2015a] Clemens Heuberger, Sara Kropf, and Helmut
           Prodinger, *Analysis of Carries in Signed Digit Expansions*,
           :arxiv:`1503.08816`.
        .. [P1964] William Parry, *Intrinsic Markov chains*, Transactions
           of the American Mathematical Society 112, 1964, pp. 55-66.
           :doi:`10.1090/S0002-9947-1964-0161372-1`.
        .. [S1948] Claude E. Shannon, *A mathematical theory of communication*,
           The Bell System Technical Journal 27, 1948, 379-423,
           :doi:`10.1002/j.1538-7305.1948.tb01338.x`.
        """
        from sage.modules.free_module_element import vector
        if not self.is_deterministic():
            raise NotImplementedError("Automaton must be deterministic.")
        if not self.digraph().is_aperiodic():
            raise NotImplementedError("Automaton must be aperiodic.")
        if not self.digraph().is_strongly_connected():
            raise NotImplementedError("Automaton must be strongly connected.")
        if not all(s.is_final for s in self.iter_states()):
            raise NotImplementedError("All states must be final.")
        M = self.adjacency_matrix().change_ring(ZZ)
        states = {state: i for i, state in enumerate(self.iter_states())}
        w_all = sorted(M.eigenvectors_right(),
                       key=lambda x: abs(x[0]),
                       reverse=True)
        w = w_all[0][1][0]
        mu = w_all[0][0]
        u_all = sorted(M.eigenvectors_left(),
                       key=lambda x: abs(x[0]),
                       reverse=True)
        u = u_all[0][1][0]
        u = 1/(u*w) * u
        final = vector(int(s.is_final) for s in self.iter_states())
        ff = u*final

        assert u*w == 1
        P = Transducer(initial_states=[s.label() for s in self.iter_initial_states()],
                       final_states=[s.label() for s in self.iter_final_states()],
                       on_duplicate_transition=duplicate_transition_add_input)
        for t in self.iter_transitions():
            P.add_transition(t.from_state.label(),
                             t.to_state.label(),
                             w[states[t.to_state]]/w[states[t.from_state]]/mu,
                             t.word_in)
        for s in self.iter_states():
            P.state(s.label()).color = 1/(w[states[s]] * ff)
            P.state(s.label()).initial_probability = w[states[s]] * u[states[s]]
        return P

    def with_output(self, word_out_function=None):
        r"""
        Construct a transducer out of this automaton.

        INPUT:

        - ``word_out_function`` -- (default: ``None``) a function. It
          transforms a :class:`transition <FSMTransition>` to the
          output word for this transition.

          If this is ``None``, then the output word will be equal to
          the input word of each transition.

        OUTPUT: a transducer

        EXAMPLES::

            sage: A = Automaton([(0, 0, 'A'), (0, 1, 'B'), (1, 2, 'C')])
            sage: T = A.with_output(); T
            Transducer with 3 states
            sage: T.transitions()
            [Transition from 0 to 0: 'A'|'A',
             Transition from 0 to 1: 'B'|'B',
             Transition from 1 to 2: 'C'|'C']

        This result is in contrast to::

            sage: Transducer(A).transitions()
            [Transition from 0 to 0: 'A'|-,
             Transition from 0 to 1: 'B'|-,
             Transition from 1 to 2: 'C'|-]

        where no output labels are created.

        Here is another example::

            sage: T2 = A.with_output(lambda t: [c.lower() for c in t.word_in])
            sage: T2.transitions()
            [Transition from 0 to 0: 'A'|'a',
             Transition from 0 to 1: 'B'|'b',
             Transition from 1 to 2: 'C'|'c']

        We can obtain the same result by composing two transducers. As inner
        transducer of the composition, we use :meth:`.with_output`
        without the optional argument
        ``word_out_function`` (which makes the output of each
        transition equal to its input); as outer transducer we use a
        :meth:`map-transducer
        <sage.combinat.finite_state_machine_generators.TransducerGenerators.map>`
        (for converting to lower case).
        This gives

        ::

            sage: L = transducers.map(lambda x: x.lower(), ['A', 'B', 'C'])
            sage: L.composition(A.with_output()).relabeled().transitions()
            [Transition from 0 to 0: 'A'|'a',
             Transition from 0 to 1: 'B'|'b',
             Transition from 1 to 2: 'C'|'c']

        .. SEEALSO::

           :meth:`.input_projection`,
           :meth:`.output_projection`,
           :class:`Transducer`,
           :meth:`transducers.map()
           <sage.combinat.finite_state_machine_generators.TransducerGenerators.map>`.

        TESTS::

            sage: A.with_output().input_projection() == A
            True
            sage: NAF = Automaton(
            ....:     {'A': [('A', 0), ('B', 1), ('B', -1)], 'B': [('A', 0)]})
            sage: NAF.with_output().input_projection() == NAF
            True
            sage: B = Automaton(
            ....:     {0: [(0, 'a'), (1, ['b', 'c']), (2, ['d', 'e'])],
            ....:      1: [(0, ['f', 'g']), (1, 'h'), (2, None)],
            ....:      2: [(0, None), (1, None), (2, ['i', 'j'])]},
            ....:     initial_states=[1, 2], final_states=[0])
            sage: B.with_output(lambda t: [c.upper() for c in t.word_in]).input_projection() == B
            True
        """
        if word_out_function is None:
            word_out_function = lambda transition: copy(transition.word_in)
        new = Transducer()
        memo = dict()
        new._copy_from_other_(self, memo=memo)
        for t in new.iter_transitions():
            t.word_out = word_out_function(t)
        return new

    def language(self, max_length=None, **kwargs):
        r"""
        Return all words accepted by this automaton.

        INPUT:

        - ``max_length`` -- integer or ``None`` (default). Only
          inputs of length at most ``max_length`` will be
          considered. If ``None``, then this iterates over all
          possible words without length restrictions.

        - ``kwargs`` -- will be passed on to the :class:`process
          iterator <FSMProcessIterator>`. See :meth:`process` for a
          description.

        OUTPUT: an iterator

        EXAMPLES::

            sage: NAF = Automaton(
            ....:     {'A': [('A', 0), ('B', 1), ('B', -1)],
            ....:      'B': [('A', 0)]},
            ....:     initial_states=['A'], final_states=['A', 'B'])
            sage: list(NAF.language(3))
            [[],
             [0], [-1], [1],
             [-1, 0], [0, 0], [1, 0], [0, -1], [0, 1],
             [-1, 0, 0], [0, -1, 0], [0, 0, 0], [0, 1, 0], [1, 0, 0],
             [-1, 0, -1], [-1, 0, 1], [0, 0, -1],
             [0, 0, 1], [1, 0, -1], [1, 0, 1]]

        .. SEEALSO::

            :meth:`FiniteStateMachine.language`,
            :meth:`process`.

        TESTS::

            sage: def R(ell):
            ....:     return (2^(ell+2)-(-1)^ell)/3
            sage: import itertools
            sage: all(len(list(NAFs)) == R(ell) for ell, NAFs in
            ....:     itertools.groupby(NAF.language(5), key=len))
            True
        """
        T = self.with_output()
        return T.language(max_length)


# ****************************************************************************


def is_Transducer(FSM):
    """
    Test whether or not ``FSM`` inherits from :class:`Transducer`.

    TESTS::

        sage: from sage.combinat.finite_state_machine import is_FiniteStateMachine, is_Transducer
        sage: is_Transducer(FiniteStateMachine())
        doctest:warning...
        DeprecationWarning: The function is_Transducer is deprecated; use 'isinstance(..., Transducer)' instead.
        See https://github.com/sagemath/sage/issues/38032 for details.
        False
        sage: is_Transducer(Transducer())
        True
        sage: is_FiniteStateMachine(Transducer())
        doctest:warning...
        DeprecationWarning: The function is_FiniteStateMachine is deprecated; use 'isinstance(..., FiniteStateMachine)' instead.
        See https://github.com/sagemath/sage/issues/38032 for details.
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(38032, "The function is_Transducer is deprecated; use 'isinstance(..., Transducer)' instead.")
    return isinstance(FSM, Transducer)


class Transducer(FiniteStateMachine):
    """
    This creates a transducer, which is a finite state machine, whose
    transitions have input and output labels.

    A transducer has additional features like creating a simplified
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
        [1]
        sage: T([1,1,0])
        [0, 0, 1]
        sage: ZZ(T(15.digits(base=2)+[0]), base=2)
        16

    Note that we have padded the binary input sequence by a `0` so
    that the transducer can reach its final state.

    TESTS::

        sage: Transducer()
        Empty transducer
    """

    def _repr_(self):
        """
        Represent the transducer as "Transducer with n states" where
        n is the number of states.

        OUTPUT: string

        EXAMPLES::

            sage: Transducer()._repr_()
            'Empty transducer'

        TESTS::

            sage: T = Transducer()
            sage: T
            Empty transducer
            sage: T.add_state(42)
            42
            sage: T
            Transducer with 1 state
            sage: T.add_state(43)
            43
            sage: T
            Transducer with 2 states
        """
        if not self._states_:
            return "Empty transducer"
        if len(self._states_) == 1:
            return "Transducer with 1 state"
        else:
            return "Transducer with %s states" % len(self._states_)

    def _latex_transition_label_(self, transition,
                                 format_function=None):
        r"""
        Return the proper transition label.

        INPUT:

        - ``transition`` -- a transition

        - ``format_function`` -- a function formatting the labels

        OUTPUT: string

        EXAMPLES::

            sage: F = Transducer([('A', 'B', 1, 2)])
            sage: print(latex(F))  # indirect doctest
            \begin{tikzpicture}[auto, initial text=, >=latex]
            \node[state] (v0) at (3.000000, 0.000000) {$\text{\texttt{A}}$};
            \node[state] (v1) at (-3.000000, 0.000000) {$\text{\texttt{B}}$};
            \path[->] (v0) edge node[rotate=360.00, anchor=south] {$1\mid 2$} (v1);
            \end{tikzpicture}

        TESTS::

            sage: F = Transducer([('A', 'B', 0, 1)])
            sage: t = F.transitions()[0]
            sage: F._latex_transition_label_(t)
            \left[0\right] \mid \left[1\right]
        """
        if format_function is None:
            format_function = latex
        return (format_function(transition.word_in) + "\\mid "
                + format_function(transition.word_out))

    def intersection(self, other, only_accessible_components=True):
        """
        Return a new transducer which accepts an input if it is accepted by
        both given finite state machines producing the same output.

        INPUT:

        - ``other`` -- a transducer

        - ``only_accessible_components`` -- if ``True`` (default), then
          the result is piped through :meth:`.accessible_components`. If no
          ``new_input_alphabet`` is given, it is determined by
          :meth:`.determine_alphabets`.

        OUTPUT:

        A new transducer which computes the intersection
        (see below) of the languages of ``self`` and ``other``.

        The set of states of the transducer is the Cartesian product of the
        set of states of both given transducer. There is a transition `((A,
        B), (C, D), a, b)` in the new transducer if there are
        transitions `(A, C, a, b)` and `(B, D, a, b)` in the old transducers.

        EXAMPLES::

            sage: transducer1 = Transducer([('1', '2', 1, 0),
            ....:                           ('2', '2', 1, 0),
            ....:                           ('2', '2', 0, 1)],
            ....:                          initial_states=['1'],
            ....:                          final_states=['2'])
            sage: transducer2 = Transducer([('A', 'A', 1, 0),
            ....:                           ('A', 'B', 0, 0),
            ....:                           ('B', 'B', 0, 1),
            ....:                           ('B', 'A', 1, 1)],
            ....:                          initial_states=['A'],
            ....:                          final_states=['B'])
            sage: res = transducer1.intersection(transducer2)
            sage: res.transitions()
            [Transition from ('1', 'A') to ('2', 'A'): 1|0,
             Transition from ('2', 'A') to ('2', 'A'): 1|0]

        In general, transducers are not closed under intersection. But
        for transducer which do not have epsilon-transitions, the
        intersection is well defined (cf. [BaWo2012]_). However, in
        the next example the intersection of the two transducers is
        not well defined. The intersection of the languages consists
        of `(a^n, b^n c^n)`. This set is not recognizable by a
        *finite* transducer.

        ::

            sage: t1 = Transducer([(0, 0, 'a', 'b'),
            ....:                  (0, 1, None, 'c'),
            ....:                  (1, 1, None, 'c')],
            ....:                 initial_states=[0],
            ....:                 final_states=[0, 1])
            sage: t2 = Transducer([('A', 'A', None, 'b'),
            ....:                  ('A', 'B', 'a', 'c'),
            ....:                  ('B', 'B', 'a', 'c')],
            ....:                 initial_states=['A'],
            ....:                 final_states=['A', 'B'])
            sage: t2.intersection(t1)
            Traceback (most recent call last):
            ...
            ValueError: An epsilon-transition (with empty input or output)
            was found.

        TESTS::

            sage: transducer1 = Transducer([('1', '2', 1, 0)],
            ....:                          initial_states=['1'],
            ....:                          final_states=['2'])
            sage: transducer2 = Transducer([('A', 'B', 1, 0)],
            ....:                          initial_states=['A'],
            ....:                          final_states=['B'])
            sage: res = transducer1.intersection(transducer2)
            sage: res.final_states()
            [('2', 'B')]
            sage: transducer1.state('2').final_word_out = 1
            sage: transducer2.state('B').final_word_out = 2
            sage: res = transducer1.intersection(transducer2)
            sage: res.final_states()
            []

        REFERENCES:

        .. [BaWo2012] Javier Baliosian and Dina Wonsever, *Finite State
           Transducers*, chapter in *Handbook of Finite State Based Models and
           Applications*, edited by Jiacun Wang, Chapman and Hall/CRC, 2012.
        """
        if not isinstance(other, Transducer):
            raise TypeError(
                "Only a transducer can be intersected with a transducer.")

        def function(transition1, transition2):
            if not transition1.word_in or not transition2.word_in \
                    or not transition1.word_out or not transition2.word_out:
                raise ValueError("An epsilon-transition "
                                 "(with empty input or output) was found.")
            if transition1.word_in == transition2.word_in \
                    and transition1.word_out == transition2.word_out:
                return (transition1.word_in, transition1.word_out)
            else:
                raise LookupError

        new = self.product_FiniteStateMachine(
               other,
               function,
               only_accessible_components=only_accessible_components,
               final_function=lambda s1, s2: s1.final_word_out)

        for state in new.iter_final_states():
            state0 = self.state(state.label()[0])
            state1 = other.state(state.label()[1])
            if state0.final_word_out != state1.final_word_out:
                state.final_word_out = None
                state.is_final = False

        return new

    def cartesian_product(self, other, only_accessible_components=True):
        r"""
        Return a new transducer which can simultaneously process an
        input with the machines ``self`` and ``other`` where the
        output labels are `d`-tuples of the original output labels.

        INPUT:

        - ``other`` -- a finite state machine (if `d=2`) or a list (or
          other iterable) of `d-1` finite state machines

        - ``only_accessible_components`` -- if ``True`` (default), then
          the result is piped through :meth:`.accessible_components`. If no
          ``new_input_alphabet`` is given, it is determined by
          :meth:`.determine_alphabets`.

        OUTPUT:

        A transducer which can simultaneously process an input with ``self``
        and the machine(s) in ``other``.

        The set of states of the new transducer is the Cartesian product of
        the set of states of ``self`` and ``other``.

        Let `(A_j, B_j, a_j, b_j)` for `j\in\{1, \ldots, d\}` be
        transitions in the machines ``self`` and in ``other``. Then
        there is a transition `((A_1, \ldots, A_d), (B_1, \ldots,
        B_d), a, (b_1, \ldots, b_d))` in the new transducer if `a_1 =
        \cdots = a_d =: a`.

        EXAMPLES::

            sage: transducer1 = Transducer([('A', 'A', 0, 0),
            ....:                           ('A', 'A', 1, 1)],
            ....:                          initial_states=['A'],
            ....:                          final_states=['A'],
            ....:                          determine_alphabets=True)
            sage: transducer2 = Transducer([(0, 1, 0, ['b', 'c']),
            ....:                           (0, 0, 1, 'b'),
            ....:                           (1, 1, 0, 'a')],
            ....:                          initial_states=[0],
            ....:                          final_states=[1],
            ....:                          determine_alphabets=True)
            sage: result = transducer1.cartesian_product(transducer2)
            sage: result
            Transducer with 2 states
            sage: result.transitions()
            [Transition from ('A', 0) to ('A', 1): 0|(0, 'b'),(None, 'c'),
             Transition from ('A', 0) to ('A', 0): 1|(1, 'b'),
             Transition from ('A', 1) to ('A', 1): 0|(0, 'a')]
            sage: result([1, 0, 0])
            [(1, 'b'), (0, 'b'), (None, 'c'),  (0, 'a')]
            sage: (transducer1([1, 0, 0]), transducer2([1, 0, 0]))
            ([1, 0, 0], ['b', 'b', 'c', 'a'])

        Also final output words are correctly processed::

            sage: transducer1.state('A').final_word_out = 2
            sage: result = transducer1.cartesian_product(transducer2)
            sage: result.final_states()[0].final_word_out
            [(2, None)]

        The following transducer counts the number of 11 blocks minus
        the number of 10 blocks over the alphabet ``[0, 1]``.

        ::

            sage: count_11 = transducers.CountSubblockOccurrences(
            ....:     [1, 1],
            ....:     input_alphabet=[0, 1])
            sage: count_10 = transducers.CountSubblockOccurrences(
            ....:     [1, 0],
            ....:     input_alphabet=[0, 1])
            sage: count_11x10 = count_11.cartesian_product(count_10)
            sage: difference = transducers.sub([0, 1])(count_11x10)
            sage: T = difference.simplification().relabeled()
            sage: T.initial_states()
            [1]
            sage: sorted(T.transitions())
            [Transition from 0 to 1: 0|-1,
             Transition from 0 to 0: 1|1,
             Transition from 1 to 1: 0|0,
             Transition from 1 to 0: 1|0]
            sage: input =  [0, 1, 1,  0, 1,  0, 0, 0, 1, 1, 1,  0]
            sage: output = [0, 0, 1, -1, 0, -1, 0, 0, 0, 1, 1, -1]
            sage: T(input) == output
            True

        If ``other`` is an automaton, then :meth:`.cartesian_product` returns
        ``self`` where the input is restricted to the input accepted by
        ``other``.

        For example, if the transducer transforms the standard
        binary expansion into the non-adjacent form and the automaton
        recognizes the binary expansion without adjacent ones, then the
        Cartesian product of these two is a transducer which does not change
        the input (except for changing ``a`` to ``(a, None)`` and ignoring a
        leading `0`).

        ::

            sage: NAF = Transducer([(0, 1, 0, None),
            ....:                   (0, 2, 1, None),
            ....:                   (1, 1, 0, 0),
            ....:                   (1, 2, 1, 0),
            ....:                   (2, 1, 0, 1),
            ....:                   (2, 3, 1, -1),
            ....:                   (3, 2, 0, 0),
            ....:                   (3, 3, 1, 0)],
            ....:                  initial_states=[0],
            ....:                  final_states=[1],
            ....:                  determine_alphabets=True)
            sage: aut11 = Automaton([(0, 0, 0), (0, 1, 1), (1, 0, 0)],
            ....:                   initial_states=[0],
            ....:                   final_states=[0, 1],
            ....:                   determine_alphabets=True)
            sage: res = NAF.cartesian_product(aut11)
            sage: res([1, 0, 0, 1, 0, 1, 0])
            [(1, None), (0, None), (0, None), (1, None), (0, None), (1, None)]

        This is obvious because if the standard binary expansion does not have
        adjacent ones, then it is the same as the non-adjacent form.

        Be aware that :meth:`.cartesian_product` is not commutative.

        ::

            sage: aut11.cartesian_product(NAF)
            Traceback (most recent call last):
            ...
            TypeError: Only an automaton can be intersected with an automaton.

        The Cartesian product of more than two finite state machines can also
        be computed::

            sage: T0 = transducers.CountSubblockOccurrences([0, 0], [0, 1, 2])
            sage: T1 = transducers.CountSubblockOccurrences([1, 1], [0, 1, 2])
            sage: T2 = transducers.CountSubblockOccurrences([2, 2], [0, 1, 2])
            sage: T = T0.cartesian_product([T1, T2])
            sage: T.transitions()
            [Transition from ((), (), ()) to ((0,), (), ()): 0|(0, 0, 0),
             Transition from ((), (), ()) to ((), (1,), ()): 1|(0, 0, 0),
             Transition from ((), (), ()) to ((), (), (2,)): 2|(0, 0, 0),
             Transition from ((0,), (), ()) to ((0,), (), ()): 0|(1, 0, 0),
             Transition from ((0,), (), ()) to ((), (1,), ()): 1|(0, 0, 0),
             Transition from ((0,), (), ()) to ((), (), (2,)): 2|(0, 0, 0),
             Transition from ((), (1,), ()) to ((0,), (), ()): 0|(0, 0, 0),
             Transition from ((), (1,), ()) to ((), (1,), ()): 1|(0, 1, 0),
             Transition from ((), (1,), ()) to ((), (), (2,)): 2|(0, 0, 0),
             Transition from ((), (), (2,)) to ((0,), (), ()): 0|(0, 0, 0),
             Transition from ((), (), (2,)) to ((), (1,), ()): 1|(0, 0, 0),
             Transition from ((), (), (2,)) to ((), (), (2,)): 2|(0, 0, 1)]
            sage: T([0, 0, 1, 1, 2, 2, 0, 1, 2, 2])
            [(0, 0, 0),
             (1, 0, 0),
             (0, 0, 0),
             (0, 1, 0),
             (0, 0, 0),
             (0, 0, 1),
             (0, 0, 0),
             (0, 0, 0),
             (0, 0, 0),
             (0, 0, 1)]
        """
        def function(*transitions):
            if equal(t.word_in for t in transitions):
                return (transitions[0].word_in,
                        list(itertools.zip_longest(
                            *(t.word_out for t in transitions)
                             )))
            else:
                raise LookupError

        def final_function(*states):
            return list(itertools.zip_longest(*(s.final_word_out
                                                 for s in states)))

        return self.product_FiniteStateMachine(
            other,
            function,
            final_function=final_function,
            only_accessible_components=only_accessible_components)

    def simplification(self):
        """
        Return a simplified transducer.

        OUTPUT: a new transducer

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
            [Transition from ('B', 'D') to ('A', 'C'): 0|0,
             Transition from ('B', 'D') to ('A', 'C'): 1|1,
             Transition from ('A', 'C') to ('B', 'D'): 0|1,
             Transition from ('A', 'C') to ('B', 'D'): 1|0]
            sage: fsms.relabeled().transitions()
            [Transition from 0 to 1: 0|0,
             Transition from 0 to 1: 1|1,
             Transition from 1 to 0: 0|1,
             Transition from 1 to 0: 1|0]

        ::

            sage: fsm = Transducer([("A", "A", 0, 0),
            ....:                   ("A", "B", 1, 1),
            ....:                   ("A", "C", 1, -1),
            ....:                   ("B", "A", 2, 0),
            ....:                   ("C", "A", 2, 0)])
            sage: fsm_simplified = fsm.simplification()
            sage: fsm_simplified
            Transducer with 2 states
            sage: fsm_simplified.transitions()
            [Transition from ('A',) to ('A',): 0|0,
             Transition from ('A',) to ('B', 'C'): 1|1,0,
             Transition from ('A',) to ('B', 'C'): 1|-1,0,
             Transition from ('B', 'C') to ('A',): 2|-]

        ::

            sage: from sage.combinat.finite_state_machine import duplicate_transition_add_input
            sage: T = Transducer([('A', 'A', 1/2, 0),
            ....:                 ('A', 'B', 1/4, 1),
            ....:                 ('A', 'C', 1/4, 1),
            ....:                 ('B', 'A', 1, 0),
            ....:                 ('C', 'A', 1, 0)],
            ....:                initial_states=[0],
            ....:                final_states=['A', 'B', 'C'],
            ....:                on_duplicate_transition=duplicate_transition_add_input)
            sage: sorted(T.simplification().transitions())
            [Transition from ('A',) to ('A',): 1/2|0,
             Transition from ('A',) to ('B', 'C'): 1/2|1,
             Transition from ('B', 'C') to ('A',): 1|0]

        Illustrating the use of colors in order to avoid identification of states::

            sage: T = Transducer( [[0,0,0,0], [0,1,1,1],
            ....:                  [1,0,0,0], [1,1,1,1]],
            ....:                 initial_states=[0],
            ....:                 final_states=[0,1])
            sage: sorted(T.simplification().transitions())
            [Transition from (0, 1) to (0, 1): 0|0,
             Transition from (0, 1) to (0, 1): 1|1]
            sage: T.state(0).color = 0
            sage: T.state(0).color = 1
            sage: sorted(T.simplification().transitions())
            [Transition from (0,) to (0,): 0|0,
             Transition from (0,) to (1,): 1|1,
             Transition from (1,) to (0,): 0|0,
             Transition from (1,) to (1,): 1|1]
        """
        fsm = deepcopy(self)
        fsm.prepone_output()
        return fsm.quotient(fsm.equivalence_classes())

    def process(self, *args, **kwargs):
        """
        Return whether the transducer accepts the input, the state
        where the computation stops and which output is generated.

        INPUT:

        - ``input_tape`` -- the input tape can be a list or an
          iterable with entries from the input alphabet. If we are
          working with a multi-tape machine (see parameter
          ``use_multitape_input`` and notes below), then the tape is a
          list or tuple of tracks, each of which can be a list or an
          iterable with entries from the input alphabet.

        - ``initial_state`` or ``initial_states`` -- the initial
          state(s) in which the machine starts. Either specify a
          single one with ``initial_state`` or a list of them with
          ``initial_states``. If both are given, ``initial_state``
          will be appended to ``initial_states``. If neither is
          specified, the initial states of the finite state machine
          are taken.

        - ``list_of_outputs`` -- (default: ``None``) a boolean or
          ``None``. If ``True``, then the outputs are given in list form
          (even if we have no or only one single output). If
          ``False``, then the result is never a list (an exception is
          raised if the result cannot be returned). If
          ``list_of_outputs=None`` the method determines automatically
          what to do (e.g. if a non-deterministic machine returns more
          than one path, then the output is returned in list form).

        - ``only_accepted`` -- boolean (default: ``False``); if set,
          then the first argument in the output is guaranteed to be
          ``True`` (if the output is a list, then the first argument
          of each element will be ``True``).

        - ``full_output`` -- boolean (default: ``True``); if set,
          then the full output is given, otherwise only the generated
          output (the third entry below only). If the input is not
          accepted, a :exc:`ValueError` is raised.

        - ``always_include_output`` -- if set (not by default), always
          include the output. This is inconsequential for a
          :class:`Transducer`, but can be used in other classes
          derived from :class:`FiniteStateMachine` where the output is
          suppressed by default, cf. :meth:`Automaton.process`.

        - ``format_output`` -- a function that translates the written
          output (which is in form of a list) to something more
          readable. By default (``None``) identity is used here.

        - ``check_epsilon_transitions`` -- boolean (default: ``True``); if
          ``False``, then epsilon transitions are not taken into consideration
          during process

        - ``write_final_word_out`` -- boolean (default: ``True``); whether the
          final output words should be written or not

        - ``use_multitape_input`` -- boolean (default: ``False``); if ``True``,
          then the multi-tape mode of the process iterator is activated. See
          also the notes below for multi-tape machines.

        - ``process_all_prefixes_of_input`` -- boolean (default: ``False``); if
          ``True``, then each prefix of the input word is processed (instead of
          processing the whole input word at once). Consequently, there is an
          output generated for each of these prefixes.

        - ``process_iterator_class`` -- (default: ``None``) a class
          inherited from :class:`FSMProcessIterator`. If ``None``,
          then :class:`FSMProcessIterator` is taken. An instance of this
          class is created and is used during the processing.

        - ``automatic_output_type`` -- boolean (default: ``False``); if set and
          the input has a parent, then the output will have the same parent. If
          the input does not have a parent, then the output will be of the same
          type as the input.

        OUTPUT:

        The full output is a triple (or a list of triples,
        cf. parameter ``list_of_outputs``), where

        - the first entry is ``True`` if the input string is accepted,

        - the second gives the reached state after processing the
          input tape (This is a state with label ``None`` if the input
          could not be processed, i.e., if at one point no
          transition to go on could be found.), and

        - the third gives a list of the output labels written during
          processing.

        If ``full_output`` is ``False``, then only the third entry
        is returned.

        Note that in the case the transducer is not
        deterministic, all possible paths are taken into account.

        This function uses an iterator which, in its simplest form, goes
        from one state to another in each step. To decide which way to
        go, it uses the input words of the outgoing transitions and
        compares them to the input tape. More precisely, in each step,
        the iterator takes an outgoing transition of the current state,
        whose input label equals the input letter of the tape. The
        output label of the transition, if present, is written on the
        output tape.

        If the choice of the outgoing transition is not unique (i.e.,
        we have a non-deterministic finite state machine), all
        possibilities are followed. This is done by splitting the
        process into several branches, one for each of the possible
        outgoing transitions.

        The process (iteration) stops if all branches are finished,
        i.e., for no branch, there is any transition whose input word
        coincides with the processed input tape. This can simply
        happen when the entire tape was read.

        Also see :meth:`~FiniteStateMachine.__call__` for a version of
        :meth:`.process` with shortened output.

        Internally this function creates and works with an instance of
        :class:`FSMProcessIterator`. This iterator can also be obtained
        with :meth:`~FiniteStateMachine.iter_process`.

        If working with multi-tape finite state machines, all input
        words of transitions are words of `k`-tuples of letters.
        Moreover, the input tape has to consist of `k` tracks, i.e.,
        be a list or tuple of `k` iterators, one for each track.

        .. WARNING::

            Working with multi-tape finite state machines is still
            experimental and can lead to wrong outputs.

        EXAMPLES::

            sage: binary_inverter = Transducer({'A': [('A', 0, 1), ('A', 1, 0)]},
            ....:                              initial_states=['A'], final_states=['A'])
            sage: binary_inverter.process([0, 1, 0, 0, 1, 1])
            (True, 'A', [1, 0, 1, 1, 0, 0])

        If we are only interested in the output, we can also use::

            sage: binary_inverter([0, 1, 0, 0, 1, 1])
            [1, 0, 1, 1, 0, 0]

        This can also be used with words as input::

            sage: # needs sage.combinat
            sage: W = Words([0, 1]); W
            Finite and infinite words over {0, 1}
            sage: w = W([0, 1, 0, 0, 1, 1]); w
            word: 010011
            sage: binary_inverter(w)
            word: 101100

        In this case it is automatically determined that the output is
        a word. The call above is equivalent to::

            sage: binary_inverter.process(w,                                            # needs sage.combinat
            ....:                         full_output=False,
            ....:                         list_of_outputs=False,
            ....:                         automatic_output_type=True)
            word: 101100

        The following transducer transforms `0^n 1` to `1^n 2`::

            sage: T = Transducer([(0, 0, 0, 1), (0, 1, 1, 2)])
            sage: T.state(0).is_initial = True
            sage: T.state(1).is_final = True

        We can see the different possibilities of the output by::

            sage: [T.process(w) for w in [[1], [0, 1], [0, 0, 1], [0, 1, 1],
            ....:                         [0], [0, 0], [2, 0], [0, 1, 2]]]
            [(True, 1, [2]), (True, 1, [1, 2]),
             (True, 1, [1, 1, 2]), (False, None, None),
             (False, 0, [1]), (False, 0, [1, 1]),
             (False, None, None), (False, None, None)]

        If we just want a condensed output, we use::

            sage: [T.process(w, full_output=False)
            ....:      for w in [[1], [0, 1], [0, 0, 1]]]
            [[2], [1, 2], [1, 1, 2]]
            sage: T.process([0], full_output=False)
            Traceback (most recent call last):
            ...
            ValueError: Invalid input sequence.
            sage: T.process([0, 1, 2], full_output=False)
            Traceback (most recent call last):
            ...
            ValueError: Invalid input sequence.

        It is equivalent to::

            sage: [T(w) for w in [[1], [0, 1], [0, 0, 1]]]
            [[2], [1, 2], [1, 1, 2]]
            sage: T([0])
            Traceback (most recent call last):
            ...
            ValueError: Invalid input sequence.
            sage: T([0, 1, 2])
            Traceback (most recent call last):
            ...
            ValueError: Invalid input sequence.

        A cycle with empty input and empty output is correctly processed::

            sage: T = Transducer([(0, 1, None, None), (1, 0, None, None)],
            ....:                initial_states=[0], final_states=[1])
            sage: T.process([])
            [(False, 0, []), (True, 1, [])]
            sage: _ = T.add_transition(-1, 0, 0, 'r')
            sage: T.state(-1).is_initial = True
            sage: T.state(0).is_initial = False
            sage: T.process([0])
            [(False, 0, ['r']), (True, 1, ['r'])]

        If there is a cycle with empty input but non-empty output, the
        possible outputs would be an infinite set::

            sage: T = Transducer([(0, 1, None, 'z'), (1, 0, None, None)],
            ....:                initial_states=[0], final_states=[1])
            sage: T.process([])
            Traceback (most recent call last):
            ...
            RuntimeError: State 0 is in an epsilon cycle (no input),
            but output is written.

        But if this cycle with empty input and non-empty output is not
        reached, the correct output is produced::

            sage: _ = T.add_transition(-1, 0, 0, 'r')
            sage: T.state(-1).is_initial = True
            sage: T.state(0).is_initial = False
            sage: T.process([])
            (False, -1, [])
            sage: T.process([0])
            Traceback (most recent call last):
            ...
            RuntimeError: State 0 is in an epsilon cycle (no input),
            but output is written.

        If we set ``check_epsilon_transitions=False``, then no
        transitions with empty input are considered
        anymore. Thus cycles with empty input are no problem anymore::

            sage: T.process([0], check_epsilon_transitions=False)
            (False, 0, ['r'])

        A simple example of a multi-tape transducer is the
        following: It writes the length of the first tape many letters ``a``
        and then the length of the second tape many letters ``b``::

            sage: M = Transducer([(0, 0, (1, None), 'a'),
            ....:                 (0, 1, [], []),
            ....:                 (1, 1, (None, 1), 'b')],
            ....:                initial_states=[0],
            ....:                final_states=[1])
            sage: M.process(([1, 1], [1]), use_multitape_input=True)
            (True, 1, ['a', 'a', 'b'])

        .. SEEALSO::

            :meth:`FiniteStateMachine.process`,
            :meth:`Automaton.process`,
            :meth:`~FiniteStateMachine.iter_process`,
            :meth:`~FiniteStateMachine.__call__`,
            :class:`FSMProcessIterator`.

        TESTS::

            sage: T = Transducer([(0, 1, 1, 'a'), (1, 0, 1, 'b')],
            ....:                initial_states=[0, 1], final_states=[1])
            sage: T.process([1, 1])
            [(False, 0, ['a', 'b']), (True, 1, ['b', 'a'])]
            sage: T.process([1, 1], T.state(0))
            (False, 0, ['a', 'b'])
            sage: T.state(1).final_word_out = 'c'
            sage: T.process([1, 1], T.state(1))
            (True, 1, ['b', 'a', 'c'])
            sage: T.process([1, 1], T.state(1), write_final_word_out=False)
            (True, 1, ['b', 'a'])

        The parameter ``input_tape`` is required::

            sage: T.process()
            Traceback (most recent call last):
            ...
            TypeError: No input tape given.
        """
        # set default values
        options = copy(self._process_default_options_)
        options.update(kwargs)

        condensed_output = (options['list_of_outputs'] is False and
                            not options['full_output'])

        if condensed_output:
            options['list_of_outputs'] = True
            options['only_accepted'] = True

        result = super().process(*args, **options)

        if (condensed_output and not result or
              not options['full_output'] and result is None):
            raise ValueError("Invalid input sequence.")
        if condensed_output and len(result) >= 2:
            raise ValueError("Found more than one accepting path.")

        if condensed_output:
            return result[0]
        return result

    def _process_convert_output_(self, output_data, **kwargs):
        """
        Helper function which converts the output of
        :meth:`FiniteStateMachine.process` to one suitable for
        transducers.

        INPUT:

        - ``output_data`` -- a triple

        - ``full_output`` -- boolean

        OUTPUT: the converted output

        TESTS::

            sage: from sage.combinat.finite_state_machine import FSMState
            sage: T = Transducer()
            sage: T._process_convert_output_((True, FSMState('a'), [1, 0, 1]),
            ....:                            full_output=False)
            [1, 0, 1]
            sage: T._process_convert_output_((True, FSMState('a'), [1, 0, 1]),
            ....:                            full_output=True)
            (True, 'a', [1, 0, 1])
        """
        accept_input, current_state, output = output_data
        if kwargs['full_output']:
            if current_state.label() is None:
                return (accept_input, current_state, None)
            else:
                return (accept_input, current_state, output)
        else:
            if not accept_input:
                return None
            return output


# ****************************************************************************


class _FSMTapeCache_(SageObject):
    """
    This is a class for caching an input tape. It is used in
    :class:`FSMProcessIterator`.

    INPUT:

    - ``tape_cache_manager`` -- list of the existing instances of
      :class:`_FSMTapeCache_`. ``self`` will be appended to this list

    - ``tape`` -- tuple or list of the input tracks (iterables)

    - ``tape_ended`` -- list of booleans (one for each track of the
      tape), which indicate whether the track iterator has already raised
      a ``StopIteration`` exception

    - ``position`` -- tuple of pairs `(p, t)` marking the current
      positions of each of the input tracks. There `p` is the number
      of letter read from track `t`. The pairs of ``position`` are
      sorted first by `p` (smallest first) and then by `t`, i.e.,
      lexicographically.

    - ``is_multitape`` -- if ``True`` each entry of the
      input-word-tuple of a transition is interpreted as word for the
      corresponding input track. If ``False`` input-words are
      interpreted as an iterable of letters.

    OUTPUT: a tape-cache

    TESTS::

        sage: from sage.combinat.finite_state_machine import _FSMTapeCache_
        sage: TC1 = _FSMTapeCache_([], (xsrange(37, 42),),
        ....:                      [False], ((0, 0),), False)
        sage: TC2 = _FSMTapeCache_([], (xsrange(37, 42), xsrange(11,15)),
        ....:                      [False, False], ((0, 0), (0, 1)), True)
        sage: TC1
        tape at 0
        sage: TC1.tape_cache_manager
        [tape at 0]
        sage: TC2
        multi-tape at (0, 0)
        sage: TC2.tape_cache_manager
        [multi-tape at (0, 0)]
    """

    def __init__(self, tape_cache_manager, tape, tape_ended,
                 position, is_multitape):
        """
        See :class:`_FSMTapeCache_` for more details.

        TESTS::

            sage: from sage.combinat.finite_state_machine import _FSMTapeCache_
            sage: TC1 = _FSMTapeCache_([], (xsrange(37, 42),),
            ....:                      [False], ((0, 0),), False)
            sage: TC2 = _FSMTapeCache_([], (xsrange(37, 42), xsrange(11,15)),
            ....:                      [False, False], ((0, 0), (0, 1)), True)
            sage: TC1m = _FSMTapeCache_([], (xsrange(37, 42),),
            ....:                       [False], ((0, 0),), True)
            sage: TC3 = _FSMTapeCache_([], (xsrange(37, 42), xsrange(11,15)),
            ....:                      [False], ((0, 0),), False)
            Traceback (most recent call last):
            ...
            TypeError: The lengths of the inputs do not match
            sage: TC4 = _FSMTapeCache_([], (xsrange(37, 42), xsrange(11,15)),
            ....:                      [False, False], ((0, 0),), False)
            Traceback (most recent call last):
            ...
            TypeError: The lengths of the inputs do not match
            sage: TC5 = _FSMTapeCache_([], (xsrange(37, 42),),
            ....:                      [False, False], ((0, 0), (0, 1)), True)
            Traceback (most recent call last):
            ...
            TypeError: The lengths of the inputs do not match
            sage: _FSMTapeCache_([], (xsrange(37, 42), xsrange(11,15)),
            ....:                [False, False], ((0, 2), (0, 1)), True)
            Traceback (most recent call last):
            ...
            TypeError: Tape position ((0, 2), (0, 1)) wrong.
        """
        if not len(tape) == len(position) == len(tape_ended):
            raise TypeError('The lengths of the inputs do not match')
        if sorted(p[1] for p in position) != list(range(len(tape))):
            raise TypeError('Tape position %s wrong.' % (position,))
        self.position = position
        self.tape = tape
        self.tape_ended = tape_ended
        self.is_multitape = is_multitape

        self.tape_cache_manager = tape_cache_manager
        self.tape_cache_manager.append(self)
        self.cache = tuple(deque() for _ in self.tape)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        OUTPUT: string

        Note that this representation depends on the parameter
        ``is_multitape`` of ``self``.

        TESTS::

            sage: from sage.combinat.finite_state_machine import _FSMTapeCache_
            sage: TC1 = _FSMTapeCache_([], (xsrange(37, 42),),
            ....:                      [False], ((0, 0),), False)
            sage: TC2 = _FSMTapeCache_([], (xsrange(37, 42), xsrange(11,15)),
            ....:                      [False, False], ((0, 0), (0, 1)), True)
            sage: TC1m = _FSMTapeCache_([], (xsrange(37, 42),),
            ....:                       [False], ((0, 0),), True)
            sage: repr(TC1)  # indirect doctest
            'tape at 0'
            sage: repr(TC1m)  # indirect doctest
            'multi-tape at (0,)'
            sage: repr(TC2)  # indirect doctest
            'multi-tape at (0, 0)'
        """
        if self.is_multitape:
            pos = tuple(p for p, t in sorted(self.position, key=lambda x: x[1]))
            return 'multi-tape at %s' % (pos,)
        else:
            return 'tape at %s' % (self.position[0][0],)

    def __deepcopy__(self, memo):
        """
        See :meth:`.deepcopy` for details.

        TESTS::

            sage: from sage.combinat.finite_state_machine import _FSMTapeCache_
            sage: TC2 = _FSMTapeCache_([], (xsrange(37, 42), xsrange(11,15)),
            ....:                      [False, False], ((0, 0), (0, 1)), True)
            sage: TC3 = deepcopy(TC2)  # indirect doctest
            sage: TC3
            multi-tape at (0, 0)
            sage: TC2.tape_cache_manager
            [multi-tape at (0, 0), multi-tape at (0, 0)]
            sage: TC2.tape_cache_manager is TC3.tape_cache_manager
            True
        """
        new = type(self)(self.tape_cache_manager,
                         self.tape, self.tape_ended,
                         self.position, self.is_multitape)
        new.cache = deepcopy(self.cache, memo)
        return new

    def deepcopy(self, memo=None):
        """
        Return a deepcopy of ``self``.

        INPUT:

        - ``memo`` -- dictionary

        OUTPUT: an instance of ``_FSMCacheTape_``

        TESTS::

            sage: from sage.combinat.finite_state_machine import _FSMTapeCache_
            sage: TC2 = _FSMTapeCache_([], (xsrange(37, 42), xsrange(11,15)),
            ....:                      [False, False], ((0, 0), (0, 1)), True)
            sage: TC3 = deepcopy(TC2)  # indirect doctest
            sage: TC2
            multi-tape at (0, 0)
            sage: TC3
            multi-tape at (0, 0)
            sage: TC2.read(0), TC2.read(1), TC2.read(1)
            ((True, 37), (True, 11), (True, 12))
            sage: TC2.preview_word()
            (37, 11)
            sage: TC2.cache is TC3.cache
            False
        """
        return deepcopy(self, memo)

    def read(self, track_number):
        """
        Reads one letter from the given track of the input tape into
        the cache.

        INPUT:

        - ``track_number`` -- integer

        OUTPUT:

        ``(True, letter)`` if reading was successful (``letter`` was
        read), otherwise ``(False, None)``.

        Note that this updates the cache of all tapes in
        ``self.tape_cache_manager``.

        TESTS::

            sage: from sage.combinat.finite_state_machine import _FSMTapeCache_
            sage: TC2 = _FSMTapeCache_([], (xsrange(37, 42), xsrange(11,15)),
            ....:                      [False, False], ((0, 0), (0, 1)), True)
            sage: TC2.read(0), TC2.read(1), TC2.read(1)
            ((True, 37), (True, 11), (True, 12))
            sage: TC2.preview_word()
            (37, 11)
            sage: TC3 = deepcopy(TC2)
            sage: TC2.cache, TC3.cache
            ((deque([37]), deque([11, 12])), (deque([37]), deque([11, 12])))
            sage: TC3.read(1)
            (True, 13)
            sage: TC2.cache, TC3.cache
            ((deque([37]), deque([11, 12, 13])),
             (deque([37]), deque([11, 12, 13])))
            sage: TC2.read(1), TC2.read(1)
            ((True, 14), (False, None))
            sage: TC2.cache
            (deque([37]), deque([11, 12, 13, 14]))
            sage: TC2.tape_ended
            [False, True]
            sage: TC2.read(1)
            (False, None)
        """
        try:
            newval = next(self.tape[track_number])
        except StopIteration:
            self.tape_ended[track_number] = True
            return (False, None)

        # update all tapes
        for tape in self.tape_cache_manager:
            tape.cache[track_number].append(newval)

        return (True, newval)

    def finished(self, track_number=None):
        r"""
        Return whether the tape (or a particular track) has reached an
        end, i.e., there are no more letters in the cache and nothing
        more to read on the original tape.

        INPUT:

        - ``track_number`` -- integer or ``None`` (default); if ``None``,
          then ``True`` is returned if all tracks are finished

        OUTPUT: boolean

        TESTS::

            sage: from sage.combinat.finite_state_machine import (
            ....:     _FSMTapeCache_, FSMTransition)
            sage: TC2 = _FSMTapeCache_([], (xsrange(37, 42), xsrange(11,15)),
            ....:                      [False, False], ((0, 0), (0, 1)), True)
            sage: while True:
            ....:     try:
            ....:         word = TC2.preview_word(return_word=True)
            ....:     except RuntimeError:
            ....:         print('stop')
            ....:         break
            ....:     print('cache: {} {}'.format(TC2.cache, TC2))
            ....:     print('finished: {} {} {}'.format(TC2.finished(),
            ....:           TC2.finished(0), TC2.finished(1)))
            ....:     TC2.forward(
            ....:         FSMTransition(0, 0, word))
            cache: (deque([37]), deque([11])) multi-tape at (0, 0)
            finished: False False False
            cache: (deque([38]), deque([12])) multi-tape at (1, 1)
            finished: False False False
            cache: (deque([39]), deque([13])) multi-tape at (2, 2)
            finished: False False False
            cache: (deque([40]), deque([14])) multi-tape at (3, 3)
            finished: False False False
            stop
            sage: print('cache: {} {}'.format(TC2.cache, TC2))
            cache: (deque([41]), deque([])) multi-tape at (4, 4)
            sage: print('finished: {} {} {}'.format(TC2.finished(),
            ....:       TC2.finished(0), TC2.finished(1)))
            finished: False False True
            sage: TC2.preview_word()
            Traceback (most recent call last):
            ...
            RuntimeError: tape reached the end
            sage: print('cache: {} {}'.format(TC2.cache, TC2))
            cache: (deque([41]), deque([])) multi-tape at (4, 4)
            sage: TC2.read(0)
            (False, None)
            sage: TC2.forward(FSMTransition(0, 0, [(0, None)]))
            sage: print('finished: {} {} {}'.format(TC2.finished(),
            ....:       TC2.finished(0), TC2.finished(1)))
            finished: True True True
        """
        if track_number is None:
            return all(self.finished(n) for n, _ in enumerate(self.cache))
        if not self.cache[track_number]:
            self.read(track_number)  # to make sure tape_ended is correct
        return self.tape_ended[track_number] and not self.cache[track_number]

    def preview_word(self, track_number=None, length=1, return_word=False):
        """
        Reads a word from the input tape.

        INPUT:

        - ``track_number`` -- integer (default: ``None``); if ``None``,
          then a tuple of words (one from each track) is returned

        - ``length`` -- (default: ``1``) the length of the word(s)

        - ``return_word`` -- boolean (default: ``False``); if set, then a word
          is returned, otherwise a single letter (in which case ``length`` has
          to be ``1``)

        OUTPUT: a single letter or a word

        A :python:`RuntimeError<library/exceptions.html#exceptions.RuntimeError>`
        is thrown if the tape (at least one track) has reached its end.

        Typically, this method is called from a hook-function of a
        state.

        The attribute ``position`` is not changed.

        TESTS::

            sage: from sage.combinat.finite_state_machine import (
            ....:     _FSMTapeCache_, FSMTransition)
            sage: TC2 = _FSMTapeCache_([], (xsrange(37, 42), xsrange(11,15)),
            ....:                      [False, False], ((0, 0), (0, 1)), True)
            sage: TC2.preview_word(), TC2.preview_word()
            ((37, 11), (37, 11))
            sage: while True:
            ....:     try:
            ....:         word = TC2.preview_word(return_word=True)
            ....:     except RuntimeError:
            ....:         print('stop')
            ....:         break
            ....:     print('read: {}'.format(word))
            ....:     print('cache: {} {}'.format(TC2.cache, TC2))
            ....:     TC2.forward(
            ....:         FSMTransition(0, 0, word))
            ....:     print('cache: {} {}'.format(TC2.cache, TC2))
            read: [(37, 11)]
            cache: (deque([37]), deque([11])) multi-tape at (0, 0)
            cache: (deque([]), deque([])) multi-tape at (1, 1)
            read: [(38, 12)]
            cache: (deque([38]), deque([12])) multi-tape at (1, 1)
            cache: (deque([]), deque([])) multi-tape at (2, 2)
            read: [(39, 13)]
            cache: (deque([39]), deque([13])) multi-tape at (2, 2)
            cache: (deque([]), deque([])) multi-tape at (3, 3)
            read: [(40, 14)]
            cache: (deque([40]), deque([14])) multi-tape at (3, 3)
            cache: (deque([]), deque([])) multi-tape at (4, 4)
            stop
            sage: print('cache: {} {}'.format(TC2.cache, TC2))
            cache: (deque([41]), deque([])) multi-tape at (4, 4)
            sage: TC2.preview_word()
            Traceback (most recent call last):
            ...
            RuntimeError: tape reached the end
            sage: print('cache: {} {}'.format(TC2.cache, TC2))
            cache: (deque([41]), deque([])) multi-tape at (4, 4)
            sage: TC2.preview_word(0)
            41
            sage: print('cache: {} {}'.format(TC2.cache, TC2))
            cache: (deque([41]), deque([])) multi-tape at (4, 4)
            sage: TC2.forward(FSMTransition(0, 0, [(41, None)]))
            sage: print('cache: {} {}'.format(TC2.cache, TC2))
            cache: (deque([]), deque([])) multi-tape at (5, 4)
        """
        if not return_word and length != 1:
            raise ValueError("Should return a letter, but parameter "
                             "length is not 1.")
        if track_number is None:
            if self.is_multitape:
                result = tuple(self.preview_word(n, length, return_word)
                               for n, _ in enumerate(self.cache))
                if len(result) != len(self.cache):
                    raise RuntimeError('tape reached the end')
                if return_word:
                    return tupleofwords_to_wordoftuples(result)
                else:
                    return result
            else:
                return self.preview_word(0, length, return_word)

        track_cache = self.cache[track_number]
        while len(track_cache) < length:
            if not self.read(track_number)[0]:
                raise RuntimeError('tape reached the end')
        if return_word:
            return list(itertools.islice(track_cache, 0, length))
        else:
            return track_cache[0]

    def compare_to_tape(self, track_number, word):
        """
        Return whether it is possible to read ``word`` from the given
        track successfully.

        INPUT:

        - ``track_number`` -- integer

        - ``word`` -- tuple or list of letters

        OUTPUT: boolean

        TESTS::

            sage: from sage.combinat.finite_state_machine import _FSMTapeCache_
            sage: TC2 = _FSMTapeCache_([], (xsrange(37, 42), xsrange(11,15)),
            ....:                      [False, False], ((0, 0), (0, 1)), True)
            sage: TC2.compare_to_tape(0, [37])
            True
            sage: TC2.compare_to_tape(1, [37])
            False
            sage: TC2.compare_to_tape(0, [37, 38])
            True
            sage: TC2.compare_to_tape(1, srange(11,15))
            True
            sage: TC2.compare_to_tape(1, srange(11,16))
            False
            sage: TC2.compare_to_tape(1, [])
            True
        """
        track_cache = self.cache[track_number]
        it_word = iter(word)

        # check letters in cache
        if any(letter_on_track != letter_in_word
               for letter_on_track, letter_in_word in zip(track_cache, it_word)):
            return False

        # check letters not already cached
        for letter_in_word in it_word:
            successful, letter_on_track = self.read(track_number)
            if not successful:
                return False
            if letter_in_word != letter_on_track:
                return False
        return True

    def forward(self, transition):
        """
        Forwards the tape according to the given transition.

        INPUT:

        - ``transition`` -- a transition of a finite state machine

        OUTPUT: nothing

        If ``self.is_multitape`` is ``False``, then this function
        forwards ``self`` (track `0`) by the number of entries of
        ``transition.word_in`` different from ``None``.
        Otherwise (if ``self.is_multitape`` is
        ``True``), this function forwards each track of ``self`` by
        the length of each entry of ``transition.word_in``. Note that
        the actual values in the input word do not play a role
        (just the length).

        This function changes the attribute ``position``.

        TESTS::

            sage: from sage.combinat.finite_state_machine import (
            ....:     _FSMTapeCache_, FSMTransition,
            ....:     tupleofwords_to_wordoftuples)
            sage: TC2 = _FSMTapeCache_([], (xsrange(37, 42), xsrange(11,15)),
            ....:                      [False, False], ((0, 0), (0, 1)), True)
            sage: TC2, TC2.cache
            (multi-tape at (0, 0), (deque([]), deque([])))
            sage: letter = TC2.preview_word(); letter
            (37, 11)
            sage: TC2, TC2.cache
            (multi-tape at (0, 0), (deque([37]), deque([11])))
            sage: TC2.forward(FSMTransition(0, 0, [letter]))
            sage: TC2, TC2.cache
            (multi-tape at (1, 1), (deque([]), deque([])))
            sage: TC2.forward(FSMTransition(0, 0, [(0, 0), (None, 0)]))
            sage: TC2, TC2.cache
            (multi-tape at (2, 3), (deque([]), deque([])))
            sage: letter = TC2.preview_word(); letter
            (39, 14)
            sage: TC2, TC2.cache
            (multi-tape at (2, 3), (deque([39]), deque([14])))
            sage: word_in = tupleofwords_to_wordoftuples([[None], [None, None]])
            sage: TC2.forward(FSMTransition(0, 0, word_in))
            sage: TC2, TC2.cache
            (multi-tape at (2, 3), (deque([39]), deque([14])))
            sage: TC2.forward(FSMTransition(0, 0, [[0, None], [None, 0]]))
            sage: TC2, TC2.cache
            (multi-tape at (3, 4), (deque([]), deque([])))
            sage: TC2.forward(FSMTransition(0, 0, [(0, 0)]))
            Traceback (most recent call last):
            ...
            ValueError: forwarding tape is not possible
        """
        def length(word):
            return len(tuple(letter for letter in word if letter is not None))

        if self.is_multitape:
            increments = tuple(length(word) for word in
                               zip(*transition.word_in))
        else:
            increments = (length(transition.word_in),)

        for track_number, (track_cache, inc) in \
                enumerate(zip(self.cache, increments)):
            for _ in range(inc):
                if not track_cache:
                    if not self.read(track_number)[0]:
                        raise ValueError('forwarding tape is not possible')
                track_cache.popleft()
        position = [(p + increments[t], t)
                    for p, t in self.position]
        self.position = tuple(sorted(position))

    def transition_possible(self, transition):
        """
        Test whether the input word of ``transition`` can be read
        from the tape.

        INPUT:

        - ``transition`` -- a transition of a finite state machine

        OUTPUT: boolean

        TESTS::

            sage: from sage.combinat.finite_state_machine import (
            ....:     _FSMTapeCache_, FSMTransition)
            sage: TC2 = _FSMTapeCache_([], (xsrange(37, 42), xsrange(11,15)),
            ....:                      [False, False], ((0, 0), (0, 1)), True)
            sage: TC2, TC2.cache
            (multi-tape at (0, 0), (deque([]), deque([])))
            sage: TC2.transition_possible(
            ....:     FSMTransition(0, 0, [(37, 11), (38, 12), (None, 13)]))
            True
            sage: TC2.transition_possible(
            ....:     FSMTransition(0, 0, [(37, 11), (38, 13)]))
            False
            sage: TC2.transition_possible(
            ....:     FSMTransition(0, 0, [(37,), (38,)]))
            Traceback (most recent call last):
            ...
            TypeError: Transition from 0 to 0: (37,),(38,)|- has bad
            input word (entries should be tuples of size 2).
        """
        if self.is_multitape:
            word_in = transition.word_in
        else:
            word_in = tupleofwords_to_wordoftuples((transition.word_in,))
        if any(len(t) != len(self.cache) for t in word_in):
            raise TypeError('%s has bad input word (entries should be '
                            'tuples of size %s).' % (transition,
                                                     len(self.cache)))
        return self._transition_possible_test_(word_in)

    def _transition_possible_epsilon_(self, word_in):
        """
        This helper function tests whether ``word_in`` equals ``epsilon``,
        i.e., whether it is the empty word or consists only of letters ``None``.

        INPUT:

        - ``word_in`` -- an input word of a transition

        OUTPUT: boolean

        TESTS::

            sage: from sage.combinat.finite_state_machine import (
            ....:     _FSMTapeCache_)
            sage: TC2 = _FSMTapeCache_([], (xsrange(37, 42), xsrange(11,15)),
            ....:                      [False, False], ((0, 0), (0, 1)), True)
            sage: TC2._transition_possible_epsilon_([])
            True
            sage: TC2._transition_possible_epsilon_([(None, None)])
            True
            sage: TC2._transition_possible_epsilon_(
            ....:     [(None, None), (None, None)])
            True
        """
        # Note that this function does not need self, but it is given
        # to be consistent with the other _transition_possible_*_
        # functions.
        return all(letter is None for t in word_in for letter in t)

    def _transition_possible_test_(self, word_in):
        """
        This helper function tests whether ``word_in`` can be read
        from the tape.

        INPUT:

        - ``word_in`` -- an input word of a transition

        OUTPUT: boolean

        This method is usually overridden in inherited classes,
        cf. :class:`_FSMTapeCacheDetectEpsilon_` and
        :class:`_FSMTapeCacheDetectAll_`.

        TESTS::

            sage: from sage.combinat.finite_state_machine import (
            ....:     _FSMTapeCache_, tupleofwords_to_wordoftuples)
            sage: TC2 = _FSMTapeCache_([], (xsrange(37, 42), xsrange(11,15)),
            ....:                      [False, False], ((0, 0), (0, 1)), True)
            sage: TC2, TC2.cache
            (multi-tape at (0, 0), (deque([]), deque([])))
            sage: word_in = tupleofwords_to_wordoftuples(
            ....:     [(37, 38), (11, 12, 13)])
            sage: TC2._transition_possible_test_(word_in)
            True
            sage: word_in = tupleofwords_to_wordoftuples(
            ....:     [(37, 38), (11, 13)])
            sage: TC2._transition_possible_test_(word_in)
            False

        Note that this function does not perform a check whether the
        input word is correct or not. This is done by the higher-level
        method :meth:`.transition_possible`::

            sage: TC2._transition_possible_test_([(37,), (38,)])
            True

        This function does not accept words of epsilon-transitions::

            sage: TC2._transition_possible_test_([])
            False
            sage: TC2._transition_possible_test_([(None, None)])
            False
        """
        if self._transition_possible_epsilon_(word_in):
            return False
        word_in_transposed = wordoftuples_to_tupleofwords(word_in)
        return all(self.compare_to_tape(track_number, word)
                   for track_number, word in enumerate(word_in_transposed))


# ****************************************************************************


class _FSMTapeCacheDetectEpsilon_(_FSMTapeCache_):
    """
    This is a class is similar to :class:`_FSMTapeCache_` but accepts
    only epsilon transitions.
    """

    def __init__(self, *args, **kwargs):
        """
        See :class:`_FSMTapeCache_` for more details.

        TESTS::

            sage: from sage.combinat.finite_state_machine import _FSMTapeCacheDetectEpsilon_
            sage: _FSMTapeCacheDetectEpsilon_([], (xsrange(37, 42),),
            ....:                             [False], ((0, 0),), False)
            tape at 0
        """
        super().__init__(*args, **kwargs)
        self._visited_states_ = set()

    def __deepcopy__(self, memo):
        """
        See :meth:`_FSMTapeCache_.deepcopy` for details.

        TESTS::

            sage: from sage.combinat.finite_state_machine import _FSMTapeCacheDetectEpsilon_
            sage: TC2 = _FSMTapeCacheDetectEpsilon_([], (xsrange(37, 42),),
            ....:                                   [False], ((0, 0),), True)
            sage: TC2._visited_states_.add(1)
            sage: TC3 = deepcopy(TC2)  # indirect doctest
            sage: TC3._visited_states_
            {1}
        """
        new = super().__deepcopy__(memo)
        new._visited_states_ = copy(self._visited_states_)
        return new

    def _transition_possible_test_(self, word_in):
        """
        This helper function tests whether ``word_in`` equals ``epsilon``,
        i.e., whether it is the empty word or consists only of letters ``None``.

        INPUT:

        - ``word_in`` -- an input word of a transition

        OUTPUT: boolean

        TESTS::

            sage: from sage.combinat.finite_state_machine import (
            ....:     _FSMTapeCacheDetectEpsilon_)
            sage: TCE = _FSMTapeCacheDetectEpsilon_([], (xsrange(37, 42), xsrange(11,15)),
            ....:                      [False, False], ((0, 0), (0, 1)), True)
            sage: TCE._transition_possible_test_([])
            True
            sage: TCE._transition_possible_test_([(None, None)])
            True
            sage: TCE._transition_possible_test_([(37, 11), (38, 12)])
            False
            sage: TCE._transition_possible_test_([(37, 11), (38, 13)])
            False
        """
        return self._transition_possible_epsilon_(word_in)


# ****************************************************************************


class _FSMTapeCacheDetectAll_(_FSMTapeCache_):
    """
    This is a class is similar to :class:`_FSMTapeCache_` but accepts
    each transition.
    """

    def compare_to_tape(self, track_number, word):
        """
        Return whether it is possible to read a word of the same length
        as ``word`` (but ignoring its actual content)
        from the given track successfully.

        INPUT:

        - ``track_number`` -- integer

        - ``word`` -- tuple or list of letters; only its length is used

        OUTPUT: boolean

        Note that this method usually returns ``True``. ``False`` can
        only be returned at the end of the input tape.

        TESTS::

            sage: from sage.combinat.finite_state_machine import _FSMTapeCacheDetectAll_
            sage: TC = _FSMTapeCacheDetectAll_(
            ....:     [], (iter((11, 12)),),
            ....:     [False], ((0, 0),), False)
            sage: TC.compare_to_tape(0, [])
            True
            sage: TC.compare_to_tape(0, [37])
            True
            sage: TC.compare_to_tape(0, [37, 38])
            True
            sage: TC.compare_to_tape(0, [37, 38, 39])
            False
            sage: TC.compare_to_tape(0, [1, 2])
            True
        """
        track_cache = self.cache[track_number]
        it_word = iter(word)

        # process letters in cache
        for _ in track_cache:
            next(it_word)

        # check letters not already cached
        for _ in it_word:
            successful, _ = self.read(track_number)
            if not successful:
                return False
        return True


# ****************************************************************************


def tupleofwords_to_wordoftuples(tupleofwords):
    """
    Transposes a tuple of words over the alphabet to a word of tuples.

    INPUT:

    - ``tupleofwords`` -- tuple of a list of letters

    OUTPUT: list of tuples

    Missing letters in the words are padded with the letter ``None``
    (from the empty word).

    EXAMPLES::

        sage: from sage.combinat.finite_state_machine import (
        ....:     tupleofwords_to_wordoftuples)
        sage: tupleofwords_to_wordoftuples(
        ....:     ([1, 2], [3, 4, 5, 6], [7]))
        [(1, 3, 7), (2, 4, None), (None, 5, None), (None, 6, None)]
    """
    return list(itertools.zip_longest(*tupleofwords, fillvalue=None))


def wordoftuples_to_tupleofwords(wordoftuples):
    """
    Transposes a word of tuples to a tuple of words over the alphabet.

    INPUT:

    - ``wordoftuples`` -- list of tuples of letters

    OUTPUT: a tuple of lists

    Letters ``None`` (empty word) are removed from each word in the output.

    EXAMPLES::

        sage: from sage.combinat.finite_state_machine import (
        ....:     wordoftuples_to_tupleofwords)
        sage: wordoftuples_to_tupleofwords(
        ....:     [(1, 2), (1, None), (1, None), (1, 2), (None, 2)])
        ([1, 1, 1, 1], [2, 2, 2])
    """
    if not equal(len(t) for t in wordoftuples):
        raise ValueError("Not all entries of input have the same length.")

    def remove_empty_letters(word):
        return [letter for letter in word if letter is not None]
    return tuple(remove_empty_letters(word)
                 for word in zip(*wordoftuples))


# ****************************************************************************


def is_FSMProcessIterator(PI):
    """
    Test whether or not ``PI`` inherits from :class:`FSMProcessIterator`.

    TESTS::

        sage: from sage.combinat.finite_state_machine import is_FSMProcessIterator, FSMProcessIterator
        sage: is_FSMProcessIterator(FSMProcessIterator(FiniteStateMachine([[0, 0, 0, 0]], initial_states=[0]), []))
        doctest:warning...
        DeprecationWarning: The function is_FSMProcessIterator is deprecated; use 'isinstance(..., FSMProcessIterator)' instead.
        See https://github.com/sagemath/sage/issues/38032 for details.
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(38032, "The function is_FSMProcessIterator is deprecated; use 'isinstance(..., FSMProcessIterator)' instead.")
    return isinstance(PI, FSMProcessIterator)


# ****************************************************************************


class FSMProcessIterator(SageObject, Iterator):
    """
    This class takes an input, feeds it into a finite state machine
    (automaton or transducer, in particular), tests whether this was
    successful and calculates the written output.

    INPUT:

    - ``fsm`` -- the finite state machine on which the input should be
      processed

    - ``input_tape`` -- the input tape can be a list or an
      iterable with entries from the input alphabet. If we are
      working with a multi-tape machine (see parameter
      ``use_multitape_input`` and notes below), then the tape is a
      list or tuple of tracks, each of which can be a list or an
      iterable with entries from the input alphabet.

    - ``initial_state`` or ``initial_states`` -- the initial
      state(s) in which the machine starts. Either specify a
      single one with ``initial_state`` or a list of them with
      ``initial_states``. If both are given, ``initial_state``
      will be appended to ``initial_states``. If neither is
      specified, the initial states of the finite state machine
      are taken.

    - ``format_output`` -- a function that translates the written
      output (which is in form of a list) to something more
      readable. By default (``None``) identity is used here.

    - ``check_epsilon_transitions`` -- boolean (default: ``True``); if
      ``False``, then epsilon transitions are not taken into consideration
      during process

    - ``write_final_word_out`` -- boolean (default: ``True``); whether the
      final output words should be written or not

    - ``use_multitape_input`` -- boolean (default: ``False``); if ``True``,
      then the multi-tape mode of the process iterator is activated. See also
      the notes below for multi-tape machines.

    - ``process_all_prefixes_of_input`` -- boolean (default: ``False``); if
      ``True``, then each prefix of the input word is processed (instead of
      processing the whole input word at once). Consequently, there is an
      output generated for each of these prefixes.

    OUTPUT: an iterator

    In its simplest form, it behaves like an iterator which, in
    each step, goes from one state to another. To decide which way
    to go, it uses the input words of the outgoing transitions and
    compares them to the input tape. More precisely, in each step,
    the process iterator takes an outgoing transition of the
    current state, whose input label equals the input letter of
    the tape. The output label of the transition, if present, is
    written on the output tape.

    If the choice of the outgoing transition is not unique (i.e.,
    we have a non-deterministic finite state machine), all
    possibilities are followed. This is done by splitting the
    process into several branches, one for each of the possible
    outgoing transitions.

    The process (iteration) stops if all branches are finished,
    i.e., for no branch, there is any transition whose input word
    coincides with the processed input tape. This can simply
    happen when the entire tape was read.
    When the process stops, a ``StopIteration`` exception is thrown.

    .. WARNING::

        Processing an input tape of length `n` usually takes at least `n+1`
        iterations, since there will be `n+1` states visited (in the
        case the taken transitions have input words consisting of single
        letters).

    An instance of this class is generated when
    :meth:`FiniteStateMachine.process` or
    :meth:`FiniteStateMachine.iter_process` of a finite state machine,
    an automaton, or a transducer is invoked.

    When working with multi-tape finite state machines, all input
    words of transitions are words of `k`-tuples of letters.
    Moreover, the input tape has to consist of `k` tracks, i.e.,
    be a list or tuple of `k` iterators, one for each track.

    .. WARNING::

        Working with multi-tape finite state machines is still
        experimental and can lead to wrong outputs.

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

    The function :meth:`FiniteStateMachine.process` (internally) uses a
    :class:`FSMProcessIterator`. We can do that manually, too, and get full
    access to the iteration process::

        sage: from sage.combinat.finite_state_machine import FSMProcessIterator
        sage: it = FSMProcessIterator(T, input_tape=input)
        sage: for current in it:
        ....:     print(current)
        process (1 branch)
        + at state 'B'
        +-- tape at 1, [[]]
        process (1 branch)
        + at state 'B'
        +-- tape at 2, [[]]
        process (1 branch)
        + at state 'A'
        +-- tape at 3, [[1, 0]]
        process (1 branch)
        + at state 'A'
        +-- tape at 4, [[1, 0, 0]]
        process (1 branch)
        + at state 'B'
        +-- tape at 5, [[1, 0, 0]]
        process (1 branch)
        + at state 'A'
        +-- tape at 6, [[1, 0, 0, 1, 0]]
        process (1 branch)
        + at state 'B'
        +-- tape at 7, [[1, 0, 0, 1, 0]]
        process (1 branch)
        + at state 'B'
        +-- tape at 8, [[1, 0, 0, 1, 0]]
        process (1 branch)
        + at state 'B'
        +-- tape at 9, [[1, 0, 0, 1, 0]]
        process (1 branch)
        + at state 'A'
        +-- tape at 10, [[1, 0, 0, 1, 0, 1, 0]]
        process (0 branches)
        sage: it.result()
        [Branch(accept=True, state='A', output=[1, 0, 0, 1, 0, 1, 0])]

    ::

        sage: T = Transducer([(0, 0, 0, 'a'), (0, 1, 0, 'b'),
        ....:                 (1, 2, 1, 'c'), (2, 0, 0, 'd'),
        ....:                 (2, 1, None, 'd')],
        ....:                initial_states=[0], final_states=[2])
        sage: sorted(T.process([0, 0, 1], format_output=lambda o: ''.join(o)))
        [(False, 1, 'abcd'), (True, 2, 'abc')]
        sage: it = FSMProcessIterator(T, input_tape=[0, 0, 1],
        ....:                         format_output=lambda o: ''.join(o))
        sage: for current in it:
        ....:     print(current)
        process (2 branches)
        + at state 0
        +-- tape at 1, [['a']]
        + at state 1
        +-- tape at 1, [['b']]
        process (2 branches)
        + at state 0
        +-- tape at 2, [['a', 'a']]
        + at state 1
        +-- tape at 2, [['a', 'b']]
        process (2 branches)
        + at state 1
        +-- tape at 3, [['a', 'b', 'c', 'd']]
        + at state 2
        +-- tape at 3, [['a', 'b', 'c']]
        process (0 branches)
        sage: sorted(it.result())
        [Branch(accept=False, state=1, output='abcd'),
         Branch(accept=True, state=2, output='abc')]

    .. SEEALSO::

        :meth:`FiniteStateMachine.process`,
        :meth:`Automaton.process`,
        :meth:`Transducer.process`,
        :meth:`FiniteStateMachine.iter_process`,
        :meth:`FiniteStateMachine.__call__`,
        :meth:`next`.

    TESTS::

        sage: T = Transducer([[0, 0, 0, 0]])
        sage: T.process([])
        Traceback (most recent call last):
        ...
        ValueError: No state is initial.

    ::

        sage: T = Transducer([[0, 1, 0, 0]], initial_states=[0, 1])
        sage: T.process([])
        [(False, 0, []), (False, 1, [])]

    ::

        sage: T = Transducer([[0, 0, 0, 0]],
        ....:                initial_states=[0], final_states=[0])
        sage: T.state(0).final_word_out = [42]
        sage: T.process([0])
        (True, 0, [0, 42])
        sage: T.process([0], write_final_word_out=False)
        (True, 0, [0])
    """

    class Current(dict):
        """
        This class stores the branches which have to be processed
        during iteration and provides a nicer formatting of them.

        This class is derived from ``dict``. It is returned by the
        ``next``-function during iteration.

        EXAMPLES:

        In the following example you can see the dict directly and
        then the nicer output provided by this class::

            sage: from sage.combinat.finite_state_machine import FSMProcessIterator
            sage: inverter = Transducer({'A': [('A', 0, 1), ('A', 1, 0)]},
            ....:     initial_states=['A'], final_states=['A'])
            sage: it = FSMProcessIterator(inverter, input_tape=[0, 1])
            sage: for current in it:
            ....:     print(dict(current))
            ....:     print(current)
            {((1, 0),): {'A': Branch(tape_cache=tape at 1, outputs=[[1]])}}
            process (1 branch)
            + at state 'A'
            +-- tape at 1, [[1]]
            {((2, 0),): {'A': Branch(tape_cache=tape at 2, outputs=[[1, 0]])}}
            process (1 branch)
            + at state 'A'
            +-- tape at 2, [[1, 0]]
            {}
            process (0 branches)
        """

        def __repr__(self):
            """
            Return a nice representation of ``self``.

            OUTPUT: string

            TESTS::

                sage: from sage.combinat.finite_state_machine import FSMProcessIterator
                sage: T = Transducer([(0, 0, 0, 0)],
                ....:     initial_states=[0], final_states=[0])
                sage: it = FSMProcessIterator(T, input_tape=[0, 0])
                sage: for current in it:                             # indirect doctest
                ....:     print(current)
                process (1 branch)
                + at state 0
                +-- tape at 1, [[0]]
                process (1 branch)
                + at state 0
                +-- tape at 2, [[0, 0]]
                process (0 branches)
            """
            data = sorted(
                (state, pos, tape_cache, outputs)
                for pos, states in self.items()
                for state, (tape_cache, outputs) in states.items())
            branch = "branch" if len(data) == 1 else "branches"
            result = "process (%s %s)" % (len(data), branch)
            for s, sdata in itertools.groupby(data, lambda x: x[0]):
                result += "\n+ at state %s" % (s,)
                for state, pos, tape_cache, outputs in sdata:
                    result += "\n+-- %s, %s" % (tape_cache, outputs)
            return result

    FinishedBranch = namedtuple('Branch', 'accept, state, output')
    r"""
    A :func:`named tuple <collections.namedtuple>` representing the
    attributes of a branch, once
    it is fully processed.
    """

    def __init__(self, fsm,
                 input_tape=None,
                 initial_state=None, initial_states=[],
                 use_multitape_input=False,
                 check_epsilon_transitions=True,
                 write_final_word_out=True,
                 format_output=None,
                 process_all_prefixes_of_input=False,
                 **kwargs):
        """
        See :class:`FSMProcessIterator` for more information.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMProcessIterator
            sage: inverter = Transducer({'A': [('A', 0, 1), ('A', 1, 0)]},
            ....:     initial_states=['A'], final_states=['A'])
            sage: it = FSMProcessIterator(inverter, input_tape=[0, 1])
            sage: for current in it:
            ....:     print(current)
            process (1 branch)
            + at state 'A'
            +-- tape at 1, [[1]]
            process (1 branch)
            + at state 'A'
            +-- tape at 2, [[1, 0]]
            process (0 branches)
            sage: it.result()
            [Branch(accept=True, state='A', output=[1, 0])]
        """
        # FSM
        self.fsm = fsm

        # multi-tape flag
        self.is_multitape = use_multitape_input

        # initial states
        initial_states = list(initial_states)
        if initial_state is not None:
            initial_states.append(initial_state)
        if not initial_states:
            initial_states = self.fsm.initial_states()
            if not initial_states:
                raise ValueError("No state is initial.")

        # input tapes
        tape = []
        if input_tape is not None:
            if self.is_multitape:
                tape.extend(input_tape)
            else:
                tape.append(input_tape)
        if not tape:
            raise TypeError('No input tape given.')
        if not all(isinstance(track, Iterable) for track in tape):
            raise TypeError('Given input tape is not iterable.')
        self._input_tape_ = tuple(iter(track) for track in tape)
        self._input_tape_ended_ = [False for _ in tape]

        # other options
        if format_output is None:
            self.format_output = list
        else:
            self.format_output = format_output

        self.check_epsilon_transitions = check_epsilon_transitions
        self.write_final_word_out = write_final_word_out
        self.process_all_prefixes_of_input = process_all_prefixes_of_input

        # init branches
        self._current_ = self.Current()
        self._current_positions_ = []  # a heap queue of the keys of _current_
        self._tape_cache_manager_ = []
        position_zero = tuple((0, t) for t, _ in enumerate(self._input_tape_))
        if not hasattr(self, 'TapeCache'):
            self.TapeCache = _FSMTapeCache_

        for state in initial_states:
            tape_cache = self.TapeCache(self._tape_cache_manager_,
                                        self._input_tape_,
                                        self._input_tape_ended_,
                                        position_zero,
                                        self.is_multitape)
            self._push_branches_(state, tape_cache, [[]])

        self._finished_ = []  # contains (accept, state, output)

    _branch_ = namedtuple('Branch', 'tape_cache, outputs')
    r"""
    A :func:`named tuple <collections.namedtuple>` representing the
    attributes of a branch at a particular state during processing.
    """

    def _push_branch_(self, state, tape_cache, outputs):
        """
        This helper function pushes a ``state`` together with
        ``tape_cache`` and ``outputs`` (i.e. a branch) to the queue
        ``self._current_``. See also :meth:`._push_branches_`.

        INPUT:

        - ``state`` -- state which has to be processed

        - ``tape_cache`` -- an instance of :class:`_FSMTapeCache_` (storing
          information what to read next)

        - ``outputs`` -- list of output tapes on each of which words
          were written until reaching ``state``

        OUTPUT: nothing

        .. NOTE::

            ``tape_cache`` is discarded if ``self.__current__`` already
            contains a branch with the same position and state.

        TESTS::

            sage: from sage.combinat.finite_state_machine import FSMProcessIterator
            sage: A = Automaton({'a': [('a', 0), ('b', 1), ('c', 1)],
            ....:                'c': [('b', None)], 'b': [('c', None)]},
            ....:     initial_states=['a'], final_states=['b', 'c'])
            sage: it = FSMProcessIterator(A, input_tape=[0, 1, 2])  # indirect doctest
            sage: it._current_
            process (1 branch)
            + at state 'a'
            +-- tape at 0, [[]]
            sage: it._push_branch_(
            ....:     A.state('b'),
            ....:     deepcopy(it._current_[((0, 0),)][A.state('a')][0]),
            ....:     [[]])
            sage: it._current_
            process (2 branches)
            + at state 'a'
            +-- tape at 0, [[]]
            + at state 'b'
            +-- tape at 0, [[]]
            sage: it._push_branches_(                               # indirect doctest
            ....:     A.state('c'),
            ....:     deepcopy(it._current_[((0, 0),)][A.state('a')][0]),
            ....:     [[]])
            sage: it._current_
            process (3 branches)
            + at state 'a'
            +-- tape at 0, [[]]
            + at state 'b'
            +-- tape at 0, [[]]
            + at state 'c'
            +-- tape at 0, [[]]

        ::

            sage: T = Transducer([(0, 1, 0, 'd'), (0, 1, 0, 'a'),
            ....:                 (0, 1, 0, 'a'), (0, 1, 0, 'a'),
            ....:                 (0, 1, 0, 'n'), (0, 1, 0, 'i'),
            ....:                 (0, 1, 0, 'e'), (0, 1, 0, 'l'),
            ....:                 (0, 1, 0, 'l'), (0, 1, 0, 'l'),
            ....:                 (1, 2, 0, ':'), (2, 3, 0, ')')],
            ....:                initial_states=[0], final_states=[3])
            sage: T.process([0, 0, 0], format_output=lambda o: ''.join(o))
            [(True, 3, 'a:)'), (True, 3, 'd:)'), (True, 3, 'e:)'),
             (True, 3, 'i:)'), (True, 3, 'l:)'), (True, 3, 'n:)')]
        """
        import heapq

        if tape_cache.position in self._current_:
            states = self._current_[tape_cache.position]
        else:
            states = self._current_[tape_cache.position] = {}
            heapq.heappush(self._current_positions_, tape_cache.position)

        if state in states:
            existing = states[state]
            new_outputs = existing.outputs
            new_outputs.extend(outputs)
            new_outputs = [t for t, _ in
                           itertools.groupby(sorted(new_outputs))]
            states[state] = FSMProcessIterator._branch_(
                existing.tape_cache, new_outputs)
        else:
            states[state] = FSMProcessIterator._branch_(tape_cache, outputs)

    def _push_branches_(self, state, tape_cache, outputs):
        """
        This function pushes a branch (consisting of a ``state``, an
        input ``tape_cache`` and ``outputs``) and one other branch for
        each epsilon successor of ``state`` to the queue (containing
        branches to process).

        INPUT:

        - ``state`` -- state which has to be processed (i.e., the
          current state, this branch is in)

        - ``tape_cache`` -- an instance of :class:`_FSMTapeCache_` (storing
          information what to read next)

        - ``outputs`` -- list of output tapes on each of which words
          were written until reaching ``state``

        OUTPUT: nothing

        When this function is called, a branch is updated, which
        means, stored for further processing. If the state has epsilon
        successors, then a new branch for each epsilon successor is
        created. All these branches start on the same position on the
        tape and get the same (more precisely, a deepcopy of the) list
        of output tapes.

        Note that ``self._current_`` contains all states which have to
        be visited in the next steps during processing. The actual
        adding of the data is done in the helper function
        :meth:`._push_branch_`.

        TESTS::

            sage: from sage.combinat.finite_state_machine import FSMProcessIterator
            sage: A = Automaton({'a': [('a', 0), ('b', 1), ('c', None)],
            ....:                'c': [('b', 2)]},
            ....:     initial_states=['a'], final_states=['b', 'c'])
            sage: it = FSMProcessIterator(A, input_tape=[0, 1, 2])  # indirect doctest
            sage: it._current_
            process (2 branches)
            + at state 'a'
            +-- tape at 0, [[]]
            + at state 'c'
            +-- tape at 0, [[]]
            sage: it._push_branches_(
            ....:     A.state('b'),
            ....:     deepcopy(it._current_[((0, 0),)][A.state('a')][0]),
            ....:     [[]])
            sage: it._current_
            process (3 branches)
            + at state 'a'
            +-- tape at 0, [[]]
            + at state 'b'
            +-- tape at 0, [[]]
            + at state 'c'
            +-- tape at 0, [[]]
        """
        self._push_branch_(state, tape_cache, outputs)
        if not self.check_epsilon_transitions:
            return
        if state._in_epsilon_cycle_(self.fsm):
            if not state._epsilon_cycle_output_empty_(self.fsm):
                raise RuntimeError(
                    'State %s is in an epsilon cycle (no input), '
                    'but output is written.' % (state,))

        for eps_state, eps_outputs in \
                state._epsilon_successors_(self.fsm).items():
            if eps_state == state:
                continue
                # "eps_state == state" means epsilon cycle
                # Since we excluded epsilon cycles where
                # output is written, this has to be one
                # which does not write output; therefore
                # skipped.
            for eps_out in eps_outputs:
                new_out = [o + list(eps_out) for o in outputs]
                self._push_branch_(eps_state, deepcopy(tape_cache), new_out)

    def __next__(self):
        """
        Makes one step in processing the input tape.

        OUTPUT:

        It returns the current status of the iterator (see below). A
        ``StopIteration`` exception is thrown when there is/was
        nothing to do (i.e. all branches ended with previous call
        of :meth:`.next`).

        The current status is a dictionary (encapsulated into an instance of
        :class:`~FSMProcessIterator.Current`).
        The keys are positions on
        the tape. The value corresponding to such a position is again
        a dictionary, where each entry represents a branch of the
        process. This dictionary maps the current state of a branch to
        a pair consisting of a tape cache and a list of output words,
        which were written during reaching this current state.

        EXAMPLES::

            sage: from sage.combinat.finite_state_machine import FSMProcessIterator
            sage: inverter = Transducer({'A': [('A', 0, 1), ('A', 1, 0)]},
            ....:     initial_states=['A'], final_states=['A'])
            sage: it = FSMProcessIterator(inverter, input_tape=[0, 1])
            sage: next(it)
            process (1 branch)
            + at state 'A'
            +-- tape at 1, [[1]]
            sage: next(it)
            process (1 branch)
            + at state 'A'
            +-- tape at 2, [[1, 0]]
            sage: next(it)
            process (0 branches)
            sage: next(it)
            Traceback (most recent call last):
            ...
            StopIteration

        .. SEEALSO::

            :meth:`FiniteStateMachine.process`,
            :meth:`Automaton.process`,
            :meth:`Transducer.process`,
            :meth:`FiniteStateMachine.iter_process`,
            :meth:`FiniteStateMachine.__call__`,
            :class:`FSMProcessIterator`.

        TESTS::

            sage: Z = Transducer()
            sage: s = Z.add_state(0)
            sage: s.is_initial = True
            sage: s.is_final = True
            sage: s.final_word_out = [1, 2]
            sage: Z.process([])
            (True, 0, [1, 2])
            sage: it = FSMProcessIterator(Z, input_tape=[])
            sage: next(it)
            process (0 branches)
            sage: next(it)
            Traceback (most recent call last):
            ...
            StopIteration

        ::

            sage: N = Transducer([(0, 0, 0, 1)], initial_states=[0])
            sage: def h_old(state, process):
            ....:     print("{} {}".format(state, process))
            sage: N.state(0).hook = h_old
            sage: N.process([0, 0])
            Traceback (most recent call last):
            ...
            ValueError: invalid input

            sage: def h_new(process, state, outputs):
            ....:     print("{} {}".format(state, outputs))
            sage: N.state(0).hook = h_new
            sage: N.process([0, 0], check_epsilon_transitions=False)
            0 [[]]
            0 [[1]]
            0 [[1, 1]]
            (False, 0, [1, 1])
        """
        import heapq

        if not self._current_:
            raise StopIteration

        def write_word(outputs, word):
            for o in outputs:
                o.extend(word)

        def step(current_state, input_tape, outputs):
            # process current state
            next_transitions = None
            state_said_finished = False
            if hasattr(current_state, 'hook'):
                if len(sage_getargspec(current_state.hook).args) == 2:
                    raise ValueError('invalid input')
                else:
                    try:
                        self._current_branch_input_tape_ = input_tape  # for preview_word
                        next_transitions = current_state.hook(
                            self, current_state, outputs)
                    except StopIteration:
                        next_transitions = []
                        state_said_finished = True
            if isinstance(next_transitions, FSMTransition):
                next_transitions = [next_transitions]
            if next_transitions is not None and \
                    not isinstance(next_transitions, Iterable):
                raise ValueError('hook of state should return a '
                                 'transition or '
                                 'a list/tuple of transitions.')

            # write output word of state
            write_word(outputs, current_state.word_out)

            # get next
            if next_transitions is None:
                next_transitions = \
                    [transition for transition in current_state.transitions
                     if input_tape.transition_possible(transition)]

            if not next_transitions:
                # this branch has to end here...
                if not (input_tape.finished() or
                        state_said_finished or
                        self.process_all_prefixes_of_input):
                    return

            if not next_transitions or self.process_all_prefixes_of_input:
                # this branch has to end here... (continued)
                successful = current_state.is_final
                if self.process_all_prefixes_of_input:
                    write_outputs = deepcopy(outputs)
                else:
                    write_outputs = outputs
                if successful and self.write_final_word_out:
                    write_word(write_outputs, current_state.final_word_out)
                for o in write_outputs:
                    self._finished_.append(
                        FSMProcessIterator.FinishedBranch(
                            accept=successful,
                            state=current_state,
                            output=self.format_output(o)))

            if not next_transitions:
                # this branch has to end here... (continued)
                return

            # at this point we know that there is at least one
            # outgoing transition to take

            new_currents = [(input_tape, outputs)]
            if len(next_transitions) > 1:
                new_currents.extend(
                    [deepcopy(new_currents[0])
                     for _ in range(len(next_transitions) - 1)])

            # process transitions
            for transition, (tape, out) in zip(next_transitions, new_currents):
                if hasattr(transition, 'hook'):
                    transition.hook(transition, self)
                write_word(out, transition.word_out)

                # go to next state
                state = transition.to_state
                tape.forward(transition)
                self._push_branches_(state, tape, out)
            return

        states_dict = self._current_.pop(heapq.heappop(self._current_positions_))
        for state, branch in states_dict.items():
            step(state, branch.tape_cache, branch.outputs)

        return self._current_

    next = __next__

    def result(self, format_output=None):
        """
        Return the already finished branches during process.

        INPUT:

        - ``format_output`` -- a function converting the output from
          list form to something more readable (default: output the
          list directly).

        OUTPUT:

        A list of triples ``(accepted, state, output)``.

        See also the parameter ``format_output`` of
        :class:`FSMProcessIterator`.

        EXAMPLES::

            sage: inverter = Transducer({'A': [('A', 0, 'one'), ('A', 1, 'zero')]},
            ....:     initial_states=['A'], final_states=['A'])
            sage: it = inverter.iter_process(input_tape=[0, 1, 1])
            sage: for _ in it:
            ....:     pass
            sage: it.result()
            [Branch(accept=True, state='A', output=['one', 'zero', 'zero'])]
            sage: it.result(lambda L: ', '.join(L))
            [(True, 'A', 'one, zero, zero')]

        Using both the parameter ``format_output`` of
        :class:`FSMProcessIterator` and the parameter ``format_output``
        of :meth:`.result` leads to concatenation of the two
        functions::

            sage: it = inverter.iter_process(input_tape=[0, 1, 1],
            ....:                            format_output=lambda L: ', '.join(L))
            sage: for _ in it:
            ....:     pass
            sage: it.result()
            [Branch(accept=True, state='A', output='one, zero, zero')]
            sage: it.result(lambda L: ', '.join(L))
            [(True, 'A', 'o, n, e, ,,  , z, e, r, o, ,,  , z, e, r, o')]
        """
        if format_output is None:
            return self._finished_
        return [r[:2] + (format_output(r[2]),) for r in self._finished_]

    def preview_word(self, track_number=None, length=1, return_word=False):
        """
        Read a word from the input tape.

        INPUT:

        - ``track_number`` -- integer (default: ``None``); if ``None``, then
          a tuple of words (one from each track) is returned

        - ``length`` -- (default: ``1``) the length of the word(s)

        - ``return_word`` -- boolean (default: ``False``); if set, then a word
          is returned, otherwise a single letter (in which case ``length`` has
          to be ``1``)

        OUTPUT: a single letter or a word

        An exception ``StopIteration`` is thrown if the tape (at least
        one track) has reached its end.

        Typically, this method is called from a hook-function of a
        state.

        EXAMPLES::

            sage: inverter = Transducer({'A': [('A', 0, 'one'),
            ....:                              ('A', 1, 'zero')]},
            ....:     initial_states=['A'], final_states=['A'])
            sage: def state_hook(process, state, output):
            ....:     print("We are now in state %s." % (state.label(),))
            ....:     try:
            ....:         w = process.preview_word()
            ....:     except RuntimeError:
            ....:         raise StopIteration
            ....:     print("Next on the tape is a %s." % (w,))
            sage: inverter.state('A').hook = state_hook
            sage: it = inverter.iter_process(
            ....:     input_tape=[0, 1, 1],
            ....:     check_epsilon_transitions=False)
            sage: for _ in it:
            ....:     pass
            We are now in state A.
            Next on the tape is a 0.
            We are now in state A.
            Next on the tape is a 1.
            We are now in state A.
            Next on the tape is a 1.
            We are now in state A.
            sage: it.result()
            [Branch(accept=True, state='A', output=['one', 'zero', 'zero'])]
        """
        return self._current_branch_input_tape_.preview_word(
            track_number, length, return_word)


# ****************************************************************************


class _FSMProcessIteratorEpsilon_(FSMProcessIterator):
    """
    This class is similar to :class:`FSMProcessIterator`, but only
    accepts epsilon transitions during process. See
    :class:`FSMProcessIterator` for more information.

    EXAMPLES::

        sage: T = Transducer([(0, 1, 0, 'a'), (0, 2, None, 'b'),
        ....:                 (2, 1, None, 'c')])
        sage: from sage.combinat.finite_state_machine import _FSMProcessIteratorEpsilon_
        sage: it = _FSMProcessIteratorEpsilon_(T, initial_state=T.state(0),
        ....:                                format_output=lambda o: ''.join(o))

    To see what is going on, we let the transducer run::

        sage: for current in it:
        ....:     print(current)
        process (1 branch)
        + at state 2
        +-- tape at 0, [['b']]
        process (1 branch)
        + at state 1
        +-- tape at 0, [['b', 'c']]
        process (0 branches)

    This class has the additional attribute ``visited_states``::

        sage: it.visited_states
        {0: [''], 1: ['bc'], 2: ['b']}

    This means the following (let us skip the state `0` for a moment):
    State `1` can be reached by a epsilon path which write ``'bc'`` as
    output. Similarly, state `2` can be reached by writing ``'b'``. We
    started in state `0`, so this is included in visited states as
    well (``''`` means that nothing was written, which is clear, since
    no path had to be taken).

    We continue with the other states as initial states::

        sage: it = _FSMProcessIteratorEpsilon_(T, initial_state=T.state(1),
        ....:                                  format_output=lambda o: ''.join(o))
        sage: for current in it:
        ....:     print(current)
        process (0 branches)
        sage: it.visited_states
        {1: ['']}
        sage: it = _FSMProcessIteratorEpsilon_(T, initial_state=T.state(2),
        ....:                                  format_output=lambda o: ''.join(o))
        sage: for current in it:
        ....:     print(current)
        process (1 branch)
        + at state 1
        +-- tape at 0, [['c']]
        process (0 branches)
        sage: it.visited_states
        {1: ['c'], 2: ['']}

    TESTS::

        sage: A = Automaton([(0, 1, 0), (1, 2, None), (2, 3, None),
        ....:                (3, 1, None), (3, 4, None), (1, 4, None)])
        sage: it = _FSMProcessIteratorEpsilon_(A, initial_state=A.state(0))
        sage: for current in it:
        ....:     print(current)
        process (0 branches)
        sage: it.visited_states
        {0: [[]]}
        sage: it = _FSMProcessIteratorEpsilon_(A, initial_state=A.state(1))
        sage: for current in it:
        ....:     print(current)
        process (2 branches)
        + at state 2
        +-- tape at 0, [[]]
        + at state 4
        +-- tape at 0, [[]]
        process (1 branch)
        + at state 3
        +-- tape at 0, [[]]
        process (1 branch)
        + at state 4
        +-- tape at 0, [[]]
        process (0 branches)
        sage: it.visited_states
        {1: [[], []], 2: [[]], 3: [[]], 4: [[], []]}

    At this point note that in the previous output, state `1` (from
    which we started) was also reached by a non-trivial
    path. Moreover, there are two different paths from `1` to `4`.

    Let us continue with the other initial states::

        sage: it = _FSMProcessIteratorEpsilon_(A, initial_state=A.state(2))
        sage: for current in it:
        ....:     print(current)
        process (1 branch)
        + at state 3
        +-- tape at 0, [[]]
        process (2 branches)
        + at state 1
        +-- tape at 0, [[]]
        + at state 4
        +-- tape at 0, [[]]
        process (1 branch)
        + at state 4
        +-- tape at 0, [[]]
        process (0 branches)
        sage: it.visited_states
        {1: [[]], 2: [[], []], 3: [[]], 4: [[], []]}
        sage: it = _FSMProcessIteratorEpsilon_(A, initial_state=A.state(3))
        sage: for current in it:
        ....:     print(current)
        process (2 branches)
        + at state 1
        +-- tape at 0, [[]]
        + at state 4
        +-- tape at 0, [[]]
        process (2 branches)
        + at state 2
        +-- tape at 0, [[]]
        + at state 4
        +-- tape at 0, [[]]
        process (0 branches)
        sage: it.visited_states
        {1: [[]], 2: [[]], 3: [[], []], 4: [[], []]}
        sage: it = _FSMProcessIteratorEpsilon_(A, initial_state=A.state(4))
        sage: for current in it:
        ....:     print(current)
        process (0 branches)
        sage: it.visited_states
        {4: [[]]}

    ::

        sage: T = Transducer([(0, 1, 0, 'a'), (1, 2, None, 'b'),
        ....:                 (2, 3, None, 'c'), (3, 1, None, 'd'),
        ....:                 (3, 4, None, 'e'), (1, 4, None, 'f')])
        sage: it = _FSMProcessIteratorEpsilon_(T, initial_state=T.state(0),
        ....:                                  format_output=lambda o: ''.join(o))
        sage: for current in it:
        ....:     print(current)
        process (0 branches)
        sage: it.visited_states
        {0: ['']}
        sage: it = _FSMProcessIteratorEpsilon_(T, initial_state=T.state(1),
        ....:                                  format_output=lambda o: ''.join(o))
        sage: for current in it:
        ....:     print(current)
        process (2 branches)
        + at state 2
        +-- tape at 0, [['b']]
        + at state 4
        +-- tape at 0, [['f']]
        process (1 branch)
        + at state 3
        +-- tape at 0, [['b', 'c']]
        process (1 branch)
        + at state 4
        +-- tape at 0, [['b', 'c', 'e']]
        process (0 branches)
        sage: it.visited_states
        {1: ['', 'bcd'], 2: ['b'],
         3: ['bc'], 4: ['f', 'bce']}
        sage: it = _FSMProcessIteratorEpsilon_(T, initial_state=T.state(2),
        ....:                                  format_output=lambda o: ''.join(o))
        sage: for current in it:
        ....:     print(current)
        process (1 branch)
        + at state 3
        +-- tape at 0, [['c']]
        process (2 branches)
        + at state 1
        +-- tape at 0, [['c', 'd']]
        + at state 4
        +-- tape at 0, [['c', 'e']]
        process (1 branch)
        + at state 4
        +-- tape at 0, [['c', 'd', 'f']]
        process (0 branches)
        sage: it.visited_states
        {1: ['cd'], 2: ['', 'cdb'],
         3: ['c'], 4: ['ce', 'cdf']}
        sage: it = _FSMProcessIteratorEpsilon_(T, initial_state=T.state(3),
        ....:                                  format_output=lambda o: ''.join(o))
        sage: for current in it:
        ....:     print(current)
        process (2 branches)
        + at state 1
        +-- tape at 0, [['d']]
        + at state 4
        +-- tape at 0, [['e']]
        process (2 branches)
        + at state 2
        +-- tape at 0, [['d', 'b']]
        + at state 4
        +-- tape at 0, [['d', 'f']]
        process (0 branches)
        sage: it.visited_states
        {1: ['d'], 2: ['db'],
         3: ['', 'dbc'], 4: ['e', 'df']}
        sage: it = _FSMProcessIteratorEpsilon_(T, initial_state=T.state(4),
        ....:                                  format_output=lambda o: ''.join(o))
        sage: for current in it:
        ....:     print(current)
        process (0 branches)
        sage: it.visited_states
        {4: ['']}

    ::

        sage: T = Transducer([(0, 1, None, 'a'), (0, 2, None, 'b'),
        ....:                 (1, 3, None, 'c'), (2, 3, None, 'd'),
        ....:                 (3, 0, None, 'e')])
        sage: it = _FSMProcessIteratorEpsilon_(T, initial_state=T.state(0),
        ....:                                  format_output=lambda o: ''.join(o))
        sage: for current in it:
        ....:     print(current)
        process (2 branches)
        + at state 1
        +-- tape at 0, [['a']]
        + at state 2
        +-- tape at 0, [['b']]
        process (1 branch)
        + at state 3
        +-- tape at 0, [['a', 'c'], ['b', 'd']]
        process (0 branches)
        sage: it.visited_states
        {0: ['', 'ace', 'bde'], 1: ['a'], 2: ['b'], 3: ['ac', 'bd']}

    ::

        sage: T = Transducer([(0, 1, None, None), (0, 2, None, 'b'),
        ....:                 (1, 3, None, None), (2, 3, None, 'd'),
        ....:                 (3, 0, None, None)])
        sage: it = _FSMProcessIteratorEpsilon_(T, initial_state=T.state(0),
        ....:                                  format_output=lambda o: ''.join(o))
        sage: for current in it:
        ....:     print(current)
        process (2 branches)
        + at state 1
        +-- tape at 0, [[]]
        + at state 2
        +-- tape at 0, [['b']]
        process (1 branch)
        + at state 3
        +-- tape at 0, [[], ['b', 'd']]
        process (0 branches)
        sage: it.visited_states
        {0: ['', '', 'bd'], 1: [''], 2: ['b'], 3: ['', 'bd']}
        sage: T.state(0)._epsilon_cycle_output_empty_(T)
        False

    ::

        sage: T = Transducer([(0, 1, None, 'a'), (1, 2, None, 'b'),
        ....:                 (0, 2, None, 'c'), (2, 3, None, 'd'),
        ....:                 (3, 0, None, 'e')])
        sage: it = _FSMProcessIteratorEpsilon_(T, initial_state=T.state(0),
        ....:                                  format_output=lambda o: ''.join(o))
        sage: for current in it:
        ....:     print(current)
        process (2 branches)
        + at state 1
        +-- tape at 0, [['a']]
        + at state 2
        +-- tape at 0, [['c']]
        process (2 branches)
        + at state 2
        +-- tape at 0, [['a', 'b']]
        + at state 3
        +-- tape at 0, [['c', 'd']]
        process (1 branch)
        + at state 3
        +-- tape at 0, [['a', 'b', 'd']]
        process (0 branches)
        sage: it.visited_states
        {0: ['', 'cde', 'abde'], 1: ['a'], 2: ['c', 'ab'], 3: ['cd', 'abd']}

    ::

        sage: T = Transducer([(0, 1, None, 'a'), (0, 2, None, 'b'),
        ....:                 (0, 2, None, 'c'), (2, 3, None, 'd'),
        ....:                 (3, 0, None, 'e')])
        sage: it = _FSMProcessIteratorEpsilon_(T, initial_state=T.state(0),
        ....:                                  format_output=lambda o: ''.join(o))
        sage: for current in it:
        ....:     print(current)
        process (2 branches)
        + at state 1
        +-- tape at 0, [['a']]
        + at state 2
        +-- tape at 0, [['b'], ['c']]
        process (1 branch)
        + at state 3
        +-- tape at 0, [['b', 'd'], ['c', 'd']]
        process (0 branches)
        sage: it.visited_states
        {0: ['', 'bde', 'cde'], 1: ['a'], 2: ['b', 'c'], 3: ['bd', 'cd']}
    """

    def __init__(self, *args, **kwargs):
        """
        See :class:`_FSMProcessIteratorEpsilon_` and
        :class:`FSMProcessIterator` for more information.

        TESTS::

            sage: T = Transducer([(0, 1, None, 'a'), (1, 2, None, 'b')])
            sage: T.state(0)._epsilon_successors_(T)  # indirect doctest
            {1: [['a']], 2: [['a', 'b']]}
        """
        kwargs['input_tape'] = iter([])
        self.TapeCache = _FSMTapeCacheDetectEpsilon_
        self.visited_states = {}
        kwargs['check_epsilon_transitions'] = False
        return super().__init__(*args, **kwargs)

    def _push_branch_(self, state, tape_cache, outputs):
        """
        This helper function does the actual adding of a ``state`` to
        ``self._current_`` (during the update of a branch), but,
        in contrast to :meth:`FSMProcessIterator._push_branch_`, it
        skips adding when the state was already visited in this branch
        (i.e. detects whether ``state`` is in an epsilon cycle).

        INPUT:

        - ``state`` -- state which has to be processed

        - ``tape_cache`` -- an instance of :class:`_FSMTapeCache_` (storing
          information what to read next)

        - ``outputs`` -- list of output tapes on each of which words
          were written until reaching ``state``

        OUTPUT: nothing

        TESTS::

            sage: T = Transducer([(0, 1, None, 'a'), (1, 2, None, 'b'),
            ....:                 (2, 0, None, 'c')])
            sage: T.state(0)._epsilon_successors_(T)  # indirect doctest
            {0: [['a', 'b', 'c']], 1: [['a']], 2: [['a', 'b']]}
            sage: T.state(1)._epsilon_successors_(T)  # indirect doctest
            {0: [['b', 'c']], 1: [['b', 'c', 'a']], 2: [['b']]}
            sage: T.state(2)._epsilon_successors_(T)  # indirect doctest
            {0: [['c']], 1: [['c', 'a']], 2: [['c', 'a', 'b']]}
        """
        if state not in self.visited_states:
            self.visited_states[state] = []
        self.visited_states[state].extend(
            self.format_output(o) for o in outputs)

        found = state in tape_cache._visited_states_
        tape_cache._visited_states_.add(state)
        if found:
            return

        super()._push_branch_(
            state, tape_cache, outputs)

        # As tape_cache may have been discarded because current already
        # contains a branch at the same state, _visited_states_ is
        # updated manually.
        tape_at_state = self._current_[tape_cache.position][state].tape_cache
        tape_at_state._visited_states_.update(tape_cache._visited_states_)


class _FSMProcessIteratorAll_(FSMProcessIterator):
    r"""
    This class is similar to :class:`FSMProcessIterator`, but
    accepts all transitions during process. See
    :class:`FSMProcessIterator` for more information.

    This is used in :meth:`FiniteStateMachine.language`.

    EXAMPLES::

        sage: F = FiniteStateMachine(
        ....:     {'A': [('A', 0, 'z'), ('B', 1, 'o'), ('B', -1, 'm')],
        ....:      'B': [('A', 0, 'z')]},
        ....:     initial_states=['A'], final_states=['A', 'B'])
        sage: from sage.combinat.finite_state_machine import _FSMProcessIteratorAll_
        sage: it = _FSMProcessIteratorAll_(F, max_length=3,
        ....:                              format_output=lambda o: ''.join(o))
        sage: for current in it:
        ....:     print(current)
        process (2 branches)
        + at state 'A'
        +-- tape at 1, [['z']]
        + at state 'B'
        +-- tape at 1, [['m'], ['o']]
        process (2 branches)
        + at state 'A'
        +-- tape at 2, [['m', 'z'], ['o', 'z'], ['z', 'z']]
        + at state 'B'
        +-- tape at 2, [['z', 'm'], ['z', 'o']]
        process (2 branches)
        + at state 'A'
        +-- tape at 3, [['m', 'z', 'z'], ['o', 'z', 'z'], ['z', 'm', 'z'],
                        ['z', 'o', 'z'], ['z', 'z', 'z']]
        + at state 'B'
        +-- tape at 3, [['m', 'z', 'm'], ['m', 'z', 'o'], ['o', 'z', 'm'],
                        ['o', 'z', 'o'], ['z', 'z', 'm'], ['z', 'z', 'o']]
        process (0 branches)
        sage: it.result()
        [Branch(accept=True, state='A', output='mzz'),
         Branch(accept=True, state='A', output='ozz'),
         Branch(accept=True, state='A', output='zmz'),
         Branch(accept=True, state='A', output='zoz'),
         Branch(accept=True, state='A', output='zzz'),
         Branch(accept=True, state='B', output='mzm'),
         Branch(accept=True, state='B', output='mzo'),
         Branch(accept=True, state='B', output='ozm'),
         Branch(accept=True, state='B', output='ozo'),
         Branch(accept=True, state='B', output='zzm'),
         Branch(accept=True, state='B', output='zzo')]
    """

    def __init__(self, *args, **kwargs):
        """
        See :class:`_FSMProcessIteratorAll_` and
        :class:`FSMProcessIterator` for more information.

        TESTS::

            sage: T = Transducer([(0, 1, 0, 'a'), (1, 2, 1, 'b')],
            ....:                initial_states=[0], final_states=[0, 1, 2])
            sage: T.determine_alphabets()
            sage: list(T.language(2))  # indirect doctest
            [[], ['a'], ['a', 'b']]
        """
        max_length = kwargs.get('max_length')
        if max_length is None:
            kwargs['input_tape'] = itertools.count()
        else:
            kwargs['input_tape'] = iter(0 for _ in range(max_length))
        self.TapeCache = _FSMTapeCacheDetectAll_
        self.visited_states = {}
        kwargs['check_epsilon_transitions'] = False
        return super().__init__(*args, **kwargs)


# ****************************************************************************


@cached_function
def setup_latex_preamble():
    r"""
    This function adds the package ``tikz`` with support for automata
    to the preamble of Latex so that the finite state machines can be
    drawn nicely.

    See the section on :ref:`finite_state_machine_LaTeX_output`
    in the introductory examples of this module.

    TESTS::

        sage: from sage.combinat.finite_state_machine import setup_latex_preamble
        sage: setup_latex_preamble()
        sage: (r"\usepackage{tikz}" in latex.extra_preamble()) == latex.has_file("tikz.sty")
        True
    """
    latex.add_package_to_preamble_if_available('tikz')
    if latex.has_file("tikz.sty"):
        latex.add_to_preamble(r'\usetikzlibrary{automata}')
