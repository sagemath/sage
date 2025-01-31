# sage.doctest: needs sage.graphs sage.modules
r"""
Common automata and transducers (finite state machines generators)

Automata and transducers in Sage can be built through the
:class:`automata <AutomatonGenerators>`
and :class:`transducers <TransducerGenerators>` objects, respectively.
It contains generators for
common finite state machines. For example,

::

    sage: I = transducers.Identity([0, 1, 2])

generates an identity transducer on the alphabet `\{0, 1, 2\}`.

To construct automata and transducers manually, you can use the
classes :class:`Automaton` and :class:`Transducer`, respectively. See
:doc:`finite_state_machine` for more details and a lot
of :ref:`examples <finite_state_machine_examples>`.

**Automata**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~AutomatonGenerators.AnyLetter` | Return an automaton recognizing any letter.
    :meth:`~AutomatonGenerators.AnyWord` | Return an automaton recognizing any word.
    :meth:`~AutomatonGenerators.EmptyWord` | Return an automaton recognizing the empty word.
    :meth:`~AutomatonGenerators.Word` | Return an automaton recognizing the given word.
    :meth:`~AutomatonGenerators.ContainsWord` | Return an automaton recognizing words containing the given word.

**Transducers**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~TransducerGenerators.Identity` | Return a transducer realizing the identity map.
    :meth:`~TransducerGenerators.abs` | Return a transducer realizing absolute value.
    :meth:`~TransducerGenerators.map` | Return a transducer realizing a function.
    :meth:`~TransducerGenerators.operator` | Return a transducer realizing a binary operation.
    :meth:`~TransducerGenerators.all` | Return a transducer realizing logical ``and``.
    :meth:`~TransducerGenerators.any` | Return a transducer realizing logical ``or``.
    :meth:`~TransducerGenerators.add` | Return a transducer realizing addition.
    :meth:`~TransducerGenerators.sub` | Return a transducer realizing subtraction.
    :meth:`~TransducerGenerators.CountSubblockOccurrences` | Return a transducer counting the occurrences of a subblock.
    :meth:`~TransducerGenerators.Wait` | Return a transducer writing ``False`` until first (or `k`-th) true input is read.
    :meth:`~TransducerGenerators.weight` | Return a transducer realizing the Hamming weight.
    :meth:`~TransducerGenerators.GrayCode` | Return a transducer realizing binary Gray code.
    :meth:`~TransducerGenerators.Recursion` | Return a transducer defined by recursions.

AUTHORS:

- Clemens Heuberger (2014-04-07): initial version
- Sara Kropf (2014-04-10): some changes in TransducerGenerator
- Daniel Krenn (2014-04-15): improved common docstring during review
- Clemens Heuberger, Daniel Krenn, Sara Kropf (2014-04-16--2014-05-02):
  A couple of improvements. Details see
  :issue:`16141`, :issue:`16142`, :issue:`16143`, :issue:`16186`.
- Sara Kropf (2014-04-29): weight transducer
- Clemens Heuberger, Daniel Krenn (2014-07-18): transducers Wait, all,
  any
- Clemens Heuberger (2014-08-10): transducer Recursion
- Clemens Heuberger (2015-07-31): automaton word
- Daniel Krenn (2015-09-14): cleanup :issue:`18227`

ACKNOWLEDGEMENT:

- Clemens Heuberger, Daniel Krenn and Sara Kropf are supported by the
  Austrian Science Fund (FWF): P 24644-N26.

Functions and methods
---------------------
"""
# ****************************************************************************
#       Copyright (C) 2014--2015 Clemens Heuberger <clemens.heuberger@aau.at>
#                     2014--2015 Daniel Krenn <dev@danielkrenn.at>
#                     2014 Sara Kropf <sara.kropf@aau.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from collections import namedtuple
import operator

from sage.combinat.finite_state_machine import Automaton, Transducer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ


class AutomatonGenerators:
    r"""
    A collection of constructors for several common automata.

    A list of all automata in this database is available via tab
    completion. Type "``automata.``" and then hit tab to see which
    automata are available.

    The automata currently in this class include:

    - :meth:`~AnyLetter`
    - :meth:`~AnyWord`
    - :meth:`~EmptyWord`
    - :meth:`~Word`
    - :meth:`~ContainsWord`
    """

    def AnyLetter(self, input_alphabet):
        r"""
        Return an automaton recognizing any letter of the given
        input alphabet.

        INPUT:

        - ``input_alphabet`` -- list; the input alphabet

        OUTPUT: an :class:`~Automaton`

        EXAMPLES::

            sage: A = automata.AnyLetter([0, 1])
            sage: A([])
            False
            sage: A([0])
            True
            sage: A([1])
            True
            sage: A([0, 0])
            False

        .. SEEALSO::

            :meth:`AnyWord`
        """
        z = ZZ(0)
        o = ZZ(1)
        return Automaton([(z, o, _) for _ in input_alphabet],
                         initial_states=[z],
                         final_states=[o])

    def AnyWord(self, input_alphabet):
        r"""
        Return an automaton recognizing any word of the given
        input alphabet.

        INPUT:

        - ``input_alphabet`` -- list; the input alphabet

        OUTPUT: an :class:`~Automaton`

        EXAMPLES::

            sage: A = automata.AnyWord([0, 1])
            sage: A([0])
            True
            sage: A([1])
            True
            sage: A([0, 1])
            True
            sage: A([0, 2])
            False

        This is equivalent to taking the :meth:`~FiniteStateMachine.kleene_star`
        of :meth:`AnyLetter` and minimizing the result. This method
        immediately gives a minimized version::

           sage: B = automata.AnyLetter([0, 1]).kleene_star().minimization().relabeled()
           sage: B == A
           True

        .. SEEALSO::

            :meth:`AnyLetter`,
            :meth:`Word`.
        """
        z = ZZ(0)
        return Automaton([(z, z, _) for _ in input_alphabet],
                         initial_states=[z],
                         final_states=[z])

    def EmptyWord(self, input_alphabet=None):
        r"""
        Return an automaton recognizing the empty word.

        INPUT:

        - ``input_alphabet`` -- iterable or ``None`` (default: ``None``)

        OUTPUT: an :class:`~Automaton`

        EXAMPLES::

            sage: A = automata.EmptyWord()
            sage: A([])
            True
            sage: A([0])
            False

        .. SEEALSO::

            :meth:`AnyLetter`,
            :meth:`AnyWord`.
        """
        z = ZZ(0)
        return Automaton(initial_states=[z],
                         final_states=[z],
                         input_alphabet=input_alphabet)

    def Word(self, word, input_alphabet=None):
        r"""
        Return an automaton recognizing the given word.

        INPUT:

        - ``word`` -- an iterable

        - ``input_alphabet`` -- list or ``None``; if ``None``,
          then the letters occurring in the word are used

        OUTPUT: an :class:`~Automaton`

        EXAMPLES::

            sage: A = automata.Word([0])
            sage: A.transitions()
            [Transition from 0 to 1: 0|-]
            sage: [A(w) for w in ([], [0], [1])]
            [False, True, False]
            sage: A = automata.Word([0, 1, 0])
            sage: A.transitions()
            [Transition from 0 to 1: 0|-,
            Transition from 1 to 2: 1|-,
            Transition from 2 to 3: 0|-]
            sage: [A(w) for w in ([], [0], [0, 1], [0, 1, 1], [0, 1, 0])]
            [False, False, False, False, True]

        If the input alphabet is not given, it is derived from the given
        word. ::

            sage: A.input_alphabet
            [0, 1]
            sage: A = automata.Word([0, 1, 0], input_alphabet=[0, 1, 2])
            sage: A.input_alphabet
            [0, 1, 2]

        .. SEEALSO::

            :meth:`AnyWord`,
            :meth:`ContainsWord`.

        TESTS::

            sage: from sage.rings.integer import Integer
            sage: all(isinstance(s.label(), Integer) for s in A.states())
            True
        """
        letters = list(word)
        length = len(letters)
        from sage.rings.integer_ring import ZZ
        return Automaton([(ZZ(i), ZZ(i + 1), letter)
                          for i, letter in enumerate(letters)],
                         initial_states=[ZZ(0)],
                         final_states=[ZZ(length)],
                         input_alphabet=input_alphabet)

    def ContainsWord(self, word, input_alphabet):
        r"""
        Return an automaton recognizing the words containing
        the given word as a factor.

        INPUT:

        - ``word`` -- list (or other iterable) of letters; the
          word we are looking for

        - ``input_alphabet`` -- list or other iterable; the input
          alphabet

        OUTPUT: an :class:`~Automaton`

        EXAMPLES::

            sage: A = automata.ContainsWord([0, 1, 0, 1, 1],
            ....:                           input_alphabet=[0, 1])
            sage: A([1, 0, 1, 0, 1, 0, 1, 1, 0, 0])
            True
            sage: A([1, 0, 1, 0, 1, 0, 1, 0])
            False

        This is equivalent to taking the concatenation of :meth:`AnyWord`,
        :meth:`Word` and :meth:`AnyWord` and minimizing the result. This
        method immediately gives a minimized version::

            sage: B = (automata.AnyWord([0, 1]) *
            ....:     automata.Word([0, 1, 0, 1, 1], [0, 1]) *
            ....:     automata.AnyWord([0, 1])).minimization()
            sage: B.is_equivalent(A)
            True

        .. SEEALSO::

            :meth:`~TransducerGenerators.CountSubblockOccurrences`,
            :meth:`AnyWord`,
            :meth:`Word`.
        """
        word = tuple(word)

        def starts_with(what, pattern):
            return len(what) >= len(pattern) \
                and what[:len(pattern)] == pattern

        def transition_function(read, input):
            if read == word:
                return (word, None)
            current = read + (input,)
            k = 0
            while not starts_with(word, current[k:]):
                k += 1
            return (current[k:], None)

        return Automaton(
            transition_function,
            input_alphabet=input_alphabet,
            initial_states=[()],
            final_states=[word])


class TransducerGenerators:
    r"""
    A collection of constructors for several common transducers.

    A list of all transducers in this database is available via tab
    completion. Type "``transducers.``" and then hit tab to see which
    transducers are available.

    The transducers currently in this class include:

    - :meth:`~Identity`
    - :meth:`~abs`
    - :meth:`~TransducerGenerators.operator`
    - :meth:`~all`
    - :meth:`~any`
    - :meth:`~add`
    - :meth:`~sub`
    - :meth:`~CountSubblockOccurrences`
    - :meth:`~Wait`
    - :meth:`~GrayCode`
    - :meth:`~Recursion`
    """

    def Identity(self, input_alphabet):
        """
        Return the identity transducer realizing the identity map.

        INPUT:

        - ``input_alphabet`` -- list or other iterable

        OUTPUT: a transducer mapping each word over ``input_alphabet`` to
        itself

        EXAMPLES::

            sage: T = transducers.Identity([0, 1])
            sage: sorted(T.transitions())
            [Transition from 0 to 0: 0|0,
             Transition from 0 to 0: 1|1]
            sage: T.initial_states()
            [0]
            sage: T.final_states()
            [0]
            sage: T.input_alphabet
            [0, 1]
            sage: T.output_alphabet
            [0, 1]
            sage: T([0, 1, 0, 1, 1])
            [0, 1, 0, 1, 1]
        """
        return Transducer(
            [(0, 0, d, d) for d in input_alphabet],
            input_alphabet=input_alphabet,
            output_alphabet=input_alphabet,
            initial_states=[0],
            final_states=[0])

    def CountSubblockOccurrences(self, block, input_alphabet):
        r"""
        Return a transducer counting the number of (possibly
        overlapping) occurrences of a block in the input.

        INPUT:

        - ``block`` -- list (or other iterable) of letters

        - ``input_alphabet`` -- list or other iterable

        OUTPUT:

        A transducer counting (in unary) the number of occurrences of the given
        block in the input.  Overlapping occurrences are counted several
        times.

        Denoting the block by `b_0\ldots b_{k-1}`, the input word by
        `i_0\ldots i_L` and the output word by `o_0\ldots o_L`, we
        have `o_j = 1` if and only if `i_{j-k+1}\ldots i_{j} = b_0\ldots
        b_{k-1}`. Otherwise, `o_j = 0`.

        EXAMPLES:

        #.  Counting the number of ``10`` blocks over the alphabet
            ``[0, 1]``::

                sage: T = transducers.CountSubblockOccurrences(
                ....:     [1, 0],
                ....:     [0, 1])
                sage: sorted(T.transitions())
                [Transition from () to (): 0|0,
                 Transition from () to (1,): 1|0,
                 Transition from (1,) to (): 0|1,
                 Transition from (1,) to (1,): 1|0]
                sage: T.input_alphabet
                [0, 1]
                sage: T.output_alphabet
                [0, 1]
                sage: T.initial_states()
                [()]
                sage: T.final_states()
                [(), (1,)]

            Check some sequence::

                sage: T([0, 1, 0, 1, 1, 0])
                [0, 0, 1, 0, 0, 1]

        #.  Counting the number of ``11`` blocks over the alphabet
            ``[0, 1]``::

                sage: T = transducers.CountSubblockOccurrences(
                ....:     [1, 1],
                ....:     [0, 1])
                sage: sorted(T.transitions())
                [Transition from () to (): 0|0,
                 Transition from () to (1,): 1|0,
                 Transition from (1,) to (): 0|0,
                 Transition from (1,) to (1,): 1|1]

            Check some sequence::

                sage: T([0, 1, 0, 1, 1, 0])
                [0, 0, 0, 0, 1, 0]

        #.  Counting the number of ``1010`` blocks over the
            alphabet ``[0, 1, 2]``::

                sage: T = transducers.CountSubblockOccurrences(
                ....:     [1, 0, 1, 0],
                ....:     [0, 1, 2])
                sage: sorted(T.transitions())
                [Transition from () to (): 0|0,
                 Transition from () to (1,): 1|0,
                 Transition from () to (): 2|0,
                 Transition from (1,) to (1, 0): 0|0,
                 Transition from (1,) to (1,): 1|0,
                 Transition from (1,) to (): 2|0,
                 Transition from (1, 0) to (): 0|0,
                 Transition from (1, 0) to (1, 0, 1): 1|0,
                 Transition from (1, 0) to (): 2|0,
                 Transition from (1, 0, 1) to (1, 0): 0|1,
                 Transition from (1, 0, 1) to (1,): 1|0,
                 Transition from (1, 0, 1) to (): 2|0]
                sage: input =  [0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 2]
                sage: output = [0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0]
                sage: T(input) == output
                True

        .. SEEALSO::

            :meth:`~AutomatonGenerators.ContainsWord`
        """
        block_as_tuple = tuple(block)

        def starts_with(what, pattern):
            return len(what) >= len(pattern) \
                and what[:len(pattern)] == pattern

        def transition_function(read, input):
            current = read + (input, )
            if starts_with(block_as_tuple, current) \
                    and len(block_as_tuple) > len(current):
                return (current, 0)
            else:
                k = 1
                while not starts_with(block_as_tuple, current[k:]):
                    k += 1
                return (current[k:], int(block_as_tuple == current))

        T = Transducer(
            transition_function,
            input_alphabet=input_alphabet,
            output_alphabet=[0, 1],
            initial_states=[()])
        for s in T.iter_states():
            s.is_final = True
        return T

    def Wait(self, input_alphabet, threshold=1):
        r"""
        Writes ``False`` until reading the ``threshold``-th occurrence
        of a true input letter; then writes ``True``.

        INPUT:

        - ``input_alphabet`` -- list or other iterable

        - ``threshold`` -- positive integer specifying how many
          occurrences of ``True`` inputs are waited for

        OUTPUT:

        A transducer writing ``False`` until the ``threshold``-th true
        (Python's standard conversion to boolean is used to convert the
        actual input to boolean) input is read. Subsequently, the
        transducer writes ``True``.

        EXAMPLES::

            sage: T = transducers.Wait([0, 1])
            sage: T([0, 0, 1, 0, 1, 0])
            [False, False, True, True, True, True]
            sage: T2 = transducers.Wait([0, 1], threshold=2)
            sage: T2([0, 0, 1, 0, 1, 0])
            [False, False, False, False, True, True]
        """
        def transition(state, input):
            if state == threshold:
                return (threshold, True)
            if not input:
                return (state, False)
            return (state + 1, state + 1 == threshold)

        T = Transducer(transition,
                       input_alphabet=input_alphabet,
                       initial_states=[0])
        for s in T.iter_states():
            s.is_final = True

        return T

    def map(self, f, input_alphabet):
        r"""
        Return a transducer which realizes a function
        on the alphabet.

        INPUT:

        - ``f`` -- function to realize

        - ``input_alphabet`` -- list or other iterable

        OUTPUT:

        A transducer mapping an input letter `x` to
        `f(x)`.

        EXAMPLES:

        The following binary transducer realizes component-wise
        absolute value (this transducer is also available as :meth:`.abs`)::

            sage: T = transducers.map(abs, [-1, 0, 1])
            sage: T.transitions()
            [Transition from 0 to 0: -1|1,
             Transition from 0 to 0: 0|0,
             Transition from 0 to 0: 1|1]
            sage: T.input_alphabet
            [-1, 0, 1]
            sage: T.initial_states()
            [0]
            sage: T.final_states()
            [0]
            sage: T([-1, 1, 0, 1])
            [1, 1, 0, 1]

        .. SEEALSO::

            :meth:`Automaton.with_output()
            <sage.combinat.finite_state_machine.Automaton.with_output>`.
        """
        return Transducer(lambda state, input: (0, f(input)),
                          input_alphabet=input_alphabet,
                          initial_states=[0],
                          final_states=[0])

    def operator(self, operator, input_alphabet, number_of_operands=2):
        r"""
        Return a transducer which realizes an operation
        on tuples over the given input alphabet.

        INPUT:

        - ``operator`` -- operator to realize; it is a function which
          takes ``number_of_operands`` input arguments (each out of
          ``input_alphabet``)

        - ``input_alphabet`` -- list or other iterable

        - ``number_of_operands`` -- (default: `2`) it specifies the number
          of input arguments the operator takes

        OUTPUT:

        A transducer mapping an input letter `(i_1, \dots, i_n)` to
        `\mathrm{operator}(i_1, \dots, i_n)`. Here, `n` equals
        ``number_of_operands``.

        The input alphabet of the generated transducer is the Cartesian
        product of ``number_of_operands`` copies of ``input_alphabet``.

        EXAMPLES:

        The following binary transducer realizes component-wise
        addition (this transducer is also available as :meth:`add`)::

            sage: import operator
            sage: T = transducers.operator(operator.add, [0, 1])
            sage: T.transitions()
            [Transition from 0 to 0: (0, 0)|0,
             Transition from 0 to 0: (0, 1)|1,
             Transition from 0 to 0: (1, 0)|1,
             Transition from 0 to 0: (1, 1)|2]
            sage: T.input_alphabet
            [(0, 0), (0, 1), (1, 0), (1, 1)]
            sage: T.initial_states()
            [0]
            sage: T.final_states()
            [0]
            sage: T([(0, 0), (0, 1), (1, 0), (1, 1)])
            [0, 1, 1, 2]

        Note that for a unary operator the input letters of the
        new transducer are tuples of length `1`::

            sage: T = transducers.operator(abs,
            ....:                          [-1, 0, 1],
            ....:                          number_of_operands=1)
            sage: T([-1, 1, 0])
            Traceback (most recent call last):
            ...
            ValueError: Invalid input sequence.
            sage: T([(-1,), (1,), (0,)])
            [1, 1, 0]

        Compare this with the transducer generated by :meth:`.map`::

            sage: T = transducers.map(abs,
            ....:                     [-1, 0, 1])
            sage: T([-1, 1, 0])
            [1, 1, 0]

        In fact, this transducer is also available as :meth:`.abs`::

            sage: T = transducers.abs([-1, 0, 1])
            sage: T([-1, 1, 0])
            [1, 1, 0]
        """
        from itertools import product

        def transition_function(state, operands):
            return (0, operator(*operands))
        pairs = list(product(input_alphabet, repeat=number_of_operands))
        return Transducer(transition_function,
                          input_alphabet=pairs,
                          initial_states=[0],
                          final_states=[0])

    def all(self, input_alphabet, number_of_operands=2):
        r"""
        Return a transducer which realizes logical ``and`` over the given
        input alphabet.

        INPUT:

        - ``input_alphabet`` -- list or other iterable

        - ``number_of_operands`` -- (default: `2`) specifies the number
          of input arguments for the ``and`` operation

        OUTPUT:

        A transducer mapping an input word
        `(i_{01}, \ldots, i_{0d})\ldots (i_{k1}, \ldots, i_{kd})` to the word
        `(i_{01} \land \cdots \land i_{0d})\ldots (i_{k1} \land \cdots \land i_{kd})`.

        The input alphabet of the generated transducer is the Cartesian
        product of ``number_of_operands`` copies of ``input_alphabet``.

        EXAMPLES:

        The following transducer realizes letter-wise
        logical ``and``::

            sage: T = transducers.all([False, True])
            sage: T.transitions()
            [Transition from 0 to 0: (False, False)|False,
             Transition from 0 to 0: (False, True)|False,
             Transition from 0 to 0: (True, False)|False,
             Transition from 0 to 0: (True, True)|True]
            sage: T.input_alphabet
            [(False, False), (False, True), (True, False), (True, True)]
            sage: T.initial_states()
            [0]
            sage: T.final_states()
            [0]
            sage: T([(False, False), (False, True), (True, False), (True, True)])
            [False, False, False, True]

        More than two operands and other input alphabets (with
        conversion to boolean) are also possible::

            sage: T3 = transducers.all([0, 1], number_of_operands=3)
            sage: T3([(0, 0, 0), (1, 0, 0), (1, 1, 1)])
            [False, False, True]
        """
        return self.operator(lambda *args: all(args),
                             input_alphabet, number_of_operands)

    def any(self, input_alphabet, number_of_operands=2):
        r"""
        Return a transducer which realizes logical ``or`` over the given
        input alphabet.

        INPUT:

        - ``input_alphabet`` -- list or other iterable

        - ``number_of_operands`` -- (default: `2`) specifies the number
          of input arguments for the ``or`` operation

        OUTPUT:

        A transducer mapping an input word
        `(i_{01}, \ldots, i_{0d})\ldots (i_{k1}, \ldots, i_{kd})` to the word
        `(i_{01} \lor \cdots \lor i_{0d})\ldots (i_{k1} \lor \cdots \lor i_{kd})`.

        The input alphabet of the generated transducer is the Cartesian
        product of ``number_of_operands`` copies of ``input_alphabet``.

        EXAMPLES:

        The following transducer realizes letter-wise
        logical ``or``::

            sage: T = transducers.any([False, True])
            sage: T.transitions()
            [Transition from 0 to 0: (False, False)|False,
             Transition from 0 to 0: (False, True)|True,
             Transition from 0 to 0: (True, False)|True,
             Transition from 0 to 0: (True, True)|True]
            sage: T.input_alphabet
            [(False, False), (False, True), (True, False), (True, True)]
            sage: T.initial_states()
            [0]
            sage: T.final_states()
            [0]
            sage: T([(False, False), (False, True), (True, False), (True, True)])
            [False, True, True, True]

        More than two operands and other input alphabets (with
        conversion to boolean) are also possible::

            sage: T3 = transducers.any([0, 1], number_of_operands=3)
            sage: T3([(0, 0, 0), (1, 0, 0), (1, 1, 1)])
            [False, True, True]
        """
        return self.operator(lambda *args: any(args),
                             input_alphabet, number_of_operands)

    def add(self, input_alphabet, number_of_operands=2):
        r"""
        Return a transducer which realizes addition on pairs over the
        given input alphabet.

        INPUT:

        - ``input_alphabet`` -- list or other iterable

        - ``number_of_operands`` -- (default: `2`) it specifies the number
          of input arguments the operator takes

        OUTPUT:

        A transducer mapping an input word
        `(i_{01}, \ldots, i_{0d})\ldots (i_{k1}, \ldots, i_{kd})` to the word
        `(i_{01} + \cdots + i_{0d})\ldots (i_{k1} + \cdots + i_{kd})`.

        The input alphabet of the generated transducer is the Cartesian
        product of ``number_of_operands`` copies of ``input_alphabet``.

        EXAMPLES:

        The following transducer realizes letter-wise
        addition::

            sage: T = transducers.add([0, 1])
            sage: T.transitions()
            [Transition from 0 to 0: (0, 0)|0,
             Transition from 0 to 0: (0, 1)|1,
             Transition from 0 to 0: (1, 0)|1,
             Transition from 0 to 0: (1, 1)|2]
            sage: T.input_alphabet
            [(0, 0), (0, 1), (1, 0), (1, 1)]
            sage: T.initial_states()
            [0]
            sage: T.final_states()
            [0]
            sage: T([(0, 0), (0, 1), (1, 0), (1, 1)])
            [0, 1, 1, 2]

        More than two operands can also be handled::

            sage: T3 = transducers.add([0, 1], number_of_operands=3)
            sage: T3.input_alphabet
            [(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1),
             (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1)]
            sage: T3([(0, 0, 0), (0, 1, 0), (0, 1, 1), (1, 1, 1)])
            [0, 1, 2, 3]
        """
        return self.operator(lambda *args: sum(args),
                             input_alphabet,
                             number_of_operands=number_of_operands)

    def sub(self, input_alphabet):
        r"""
        Return a transducer which realizes subtraction on pairs over
        the given input alphabet.

        INPUT:

        - ``input_alphabet`` -- list or other iterable

        OUTPUT:

        A transducer mapping an input word `(i_0, i'_0)\ldots (i_k, i'_k)`
        to the word `(i_0 - i'_0)\ldots (i_k - i'_k)`.

        The input alphabet of the generated transducer is the Cartesian
        product of two copies of ``input_alphabet``.

        EXAMPLES:

        The following transducer realizes letter-wise
        subtraction::

            sage: T = transducers.sub([0, 1])
            sage: T.transitions()
            [Transition from 0 to 0: (0, 0)|0,
             Transition from 0 to 0: (0, 1)|-1,
             Transition from 0 to 0: (1, 0)|1,
             Transition from 0 to 0: (1, 1)|0]
            sage: T.input_alphabet
            [(0, 0), (0, 1), (1, 0), (1, 1)]
            sage: T.initial_states()
            [0]
            sage: T.final_states()
            [0]
            sage: T([(0, 0), (0, 1), (1, 0), (1, 1)])
            [0, -1, 1, 0]
        """
        return self.operator(operator.sub, input_alphabet)

    def weight(self, input_alphabet, zero=0):
        r"""
        Return a transducer which realizes the Hamming weight of the input
        over the given input alphabet.

        INPUT:

        - ``input_alphabet`` -- list or other iterable

        - ``zero`` -- the zero symbol in the alphabet used

        OUTPUT:

        A transducer mapping `i_0\ldots i_k` to `(i_0\neq 0)\ldots(i_k\neq 0)`.

        The Hamming weight is defined as the number of nonzero digits in the
        input sequence over the alphabet ``input_alphabet`` (see
        :wikipedia:`Hamming_weight`). The output sequence of the transducer is
        a unary encoding of the Hamming weight. Thus the sum of the output
        sequence is the Hamming weight of the input.

        EXAMPLES::

            sage: W = transducers.weight([-1, 0, 2])
            sage: W.transitions()
            [Transition from 0 to 0: -1|1,
             Transition from 0 to 0: 0|0,
             Transition from 0 to 0: 2|1]
            sage: unary_weight = W([-1, 0, 0, 2, -1])
            sage: unary_weight
            [1, 0, 0, 1, 1]
            sage: weight = add(unary_weight)
            sage: weight
            3

        Also the joint Hamming weight can be computed::

            sage: v1 = vector([-1, 0])
            sage: v0 = vector([0, 0])
            sage: W = transducers.weight([v1, v0])
            sage: unary_weight = W([v1, v0, v1, v0])
            sage: add(unary_weight)
            2

        For the input alphabet ``[-1, 0, 1]`` the weight transducer is the
        same as the absolute value transducer
        :meth:`~TransducerGenerators.abs`::

            sage: W = transducers.weight([-1, 0, 1])
            sage: A = transducers.abs([-1, 0, 1])
            sage: W == A
            True

        For other input alphabets, we can specify the zero symbol::

            sage: W = transducers.weight(['a', 'b'], zero='a')
            sage: add(W(['a', 'b', 'b']))
            2
        """
        def weight(state, input):
            weight = int(input != zero)
            return (0, weight)
        return Transducer(weight, input_alphabet=input_alphabet,
                          initial_states=[0],
                          final_states=[0])

    def abs(self, input_alphabet):
        r"""
        Return a transducer which realizes the letter-wise
        absolute value of an input word over the given input alphabet.

        INPUT:

        - ``input_alphabet`` -- list or other iterable

        OUTPUT:

        A transducer mapping `i_0\ldots i_k`
        to `|i_0|\ldots |i_k|`.

        EXAMPLES:

        The following transducer realizes letter-wise
        absolute value::

            sage: T = transducers.abs([-1, 0, 1])
            sage: T.transitions()
            [Transition from 0 to 0: -1|1,
             Transition from 0 to 0: 0|0,
             Transition from 0 to 0: 1|1]
            sage: T.initial_states()
            [0]
            sage: T.final_states()
            [0]
            sage: T([-1, -1, 0, 1])
            [1, 1, 0, 1]
        """
        return self.map(abs, input_alphabet)

    def GrayCode(self):
        """
        Return a transducer converting the standard binary
        expansion to Gray code.

        OUTPUT: a transducer

        Cf. the :wikipedia:`Gray_code` for a description of the Gray code.

        EXAMPLES::

            sage: G = transducers.GrayCode()
            sage: G
            Transducer with 3 states
            sage: for v in srange(10):
            ....:     print("{} {}".format(v, G(v.digits(base=2))))
            0 []
            1 [1]
            2 [1, 1]
            3 [0, 1]
            4 [0, 1, 1]
            5 [1, 1, 1]
            6 [1, 0, 1]
            7 [0, 0, 1]
            8 [0, 0, 1, 1]
            9 [1, 0, 1, 1]

        In the example :ref:`Gray Code <finite_state_machine_gray_code_example>`
        in the documentation of the
        :doc:`finite_state_machine` module, the Gray code
        transducer is derived from the algorithm converting the binary
        expansion to the Gray code. The result is the same as the one
        given here.
        """
        z = ZZ(0)
        o = ZZ(1)
        return Transducer([[0, 1, z, None],
                           [0, 2, o, None],
                           [1, 1, z, z],
                           [1, 2, o, o],
                           [2, 1, z, o],
                           [2, 2, o, z]],
                          initial_states=[0],
                          final_states=[1],
                          with_final_word_out=[0])

    RecursionRule = namedtuple('RecursionRule', ['K', 'r', 'k', 's', 't'])

    def _parse_recursion_equation_(self, equation, base, function, var,
                                   word_function=None, output_rings=[ZZ, QQ]):
        """
        Parse one equation as admissible in :meth:`~.Recursion`.

        INPUT:

        - ``equation`` -- an equation of the form

          - ``f(base^K * n + r) == f(base^k * n + s) + t`` for some
            integers ``0 <= k < K``, ``r`` and some ``t``---valid for
            all ``n`` such that the arguments on both sides are
            nonnegative---

          or the form

          - ``f(r) == t`` for some integer ``r`` and some ``t``.

        - ``base`` -- see :meth:`~Recursion`

        - ``function`` -- see :meth:`~Recursion`

        - ``var`` -- see :meth:`~Recursion`

        - ``output_rings`` -- see :meth:`~Recursion`

        OUTPUT:

        A ``RecursionRule`` if the equation is of the first form
        described above and a dictionary ``{r: [t]}`` otherwise.

        EXAMPLES::

            sage: # needs sage.symbolic
            sage: var('n')
            n
            sage: function('f')
            f
            sage: transducers._parse_recursion_equation_(
            ....:     f(8*n + 7) == f(2*n + 3) + 5,
            ....:     2, f, n)
            RecursionRule(K=3, r=7, k=1, s=3, t=[5])
            sage: transducers._parse_recursion_equation_(
            ....:     f(42) == 5,
            ....:     2, f, n)
            {42: [5]}

        TESTS:

            The following tests check that the equations are well-formed::

                sage: transducers._parse_recursion_equation_(f(4*n + 1), 2, f, n)       # needs sage.symbolic
                Traceback (most recent call last):
                ...
                ValueError: f(4*n + 1) is not an equation with ==.

            ::

                sage: transducers._parse_recursion_equation_(f(n) + 1 == f(2*n),        # needs sage.symbolic
                ....:     2, f, n)
                Traceback (most recent call last):
                ...
                ValueError: f(n) + 1 is not an evaluation of f.

            ::

                sage: transducers._parse_recursion_equation_(f(2*n, 5) == 3,            # needs sage.symbolic
                ....:     2, f, n)
                Traceback (most recent call last):
                ...
                ValueError: f(2*n, 5) does not have one argument.

            ::

                sage: transducers._parse_recursion_equation_(f(1/n) == f(n) + 3,        # needs sage.symbolic
                ....:     2, f, n)
                Traceback (most recent call last):
                ...
                ValueError: 1/n is not a polynomial in n.

            ::

                sage: transducers._parse_recursion_equation_(f(n^2 + 5) == 3,           # needs sage.symbolic
                ....:     2, f, n)
                Traceback (most recent call last):
                ...
                ValueError: n^2 + 5 is not a polynomial of degree 1.

            ::

                sage: transducers._parse_recursion_equation_(f(3*n + 5) == f(n) + 7,    # needs sage.symbolic
                ....:     2, f, n)
                Traceback (most recent call last):
                ...
                ValueError: 3 is not a power of 2.

            ::

                sage: transducers._parse_recursion_equation_(f(n + 5) == f(n) + 7,      # needs sage.symbolic
                ....:     2, f, n)
                Traceback (most recent call last):
                ...
                ValueError: 1 is less than 2.

            ::

                sage: transducers._parse_recursion_equation_(                           # needs sage.symbolic
                ....:     f(2*n + 1) == f(n + 1) + f(n) + 2,
                ....:     2, f, n)
                Traceback (most recent call last):
                ...
                ValueError: f(n + 1) + f(n) + 2 does not contain
                exactly one summand which is an evaluation of f.

            ::

                sage: transducers._parse_recursion_equation_(f(2*n + 1) == sin(n) + 2,  # needs sage.symbolic
                ....:     2, f, n)
                Traceback (most recent call last):
                ...
                ValueError: sin(n) + 2 does not contain exactly one
                summand which is an evaluation of f.

            ::

                sage: transducers._parse_recursion_equation_(                           # needs sage.symbolic
                ....:     f(2*n + 1) == f(n) + n + 2,
                ....:     2, f, n)
                Traceback (most recent call last):
                ...
                ValueError: n + 2 contains n.

            ::

                sage: transducers._parse_recursion_equation_(f(2*n + 1) == sin(n),      # needs sage.symbolic
                ....:     2, f, n)
                Traceback (most recent call last):
                ...
                ValueError: sin(n) is not an evaluation of f.

            ::

                sage: transducers._parse_recursion_equation_(f(2*n + 1) == f(n, 2),     # needs sage.symbolic
                ....:     2, f, n)
                Traceback (most recent call last):
                ...
                ValueError: f(n, 2) does not have exactly one argument.

            ::

                sage: transducers._parse_recursion_equation_(f(2*n + 1) == f(1/n),      # needs sage.symbolic
                ....:     2, f, n)
                Traceback (most recent call last):
                ...
                ValueError: 1/n is not a polynomial in n.

            ::

                sage: transducers._parse_recursion_equation_(f(2*n + 1) == f(n^2 + 5),  # needs sage.symbolic
                ....:     2, f, n)
                Traceback (most recent call last):
                ...
                ValueError: n^2 + 5 is not a polynomial of degree 1.

            ::

                sage: transducers._parse_recursion_equation_(f(2*n + 1) == f(3*n + 5),  # needs sage.symbolic
                ....:     2, f, n)
                Traceback (most recent call last):
                ...
                ValueError: 3 is not a power of 2.

            ::

                sage: transducers._parse_recursion_equation_(                           # needs sage.symbolic
                ....:     f(2*n + 1) == f((1/2)*n + 5),
                ....:     QQ(2), f, n)
                Traceback (most recent call last):
                ...
                ValueError: 1/2 is less than 1.

            ::

                sage: transducers._parse_recursion_equation_(f(2*n + 1) == f(2*n + 5),  # needs sage.symbolic
                ....:     2, f, n)
                Traceback (most recent call last):
                ...
                ValueError: 2 is greater or equal than 2.
        """
        from sage.functions.log import log

        def is_scalar(expression):
            return var not in expression.variables()

        def convert_output(output):
            for ring in output_rings:
                try:
                    return ring(output)
                except (ValueError, TypeError):
                    pass
            return output

        def to_list(output):
            if output == 0:
                return []
            elif word_function is not None and output.operator() == word_function:
                return [convert_output(_) for _ in output.operands()]
            else:
                return [convert_output(output)]

        base_ring = base.parent()

        if equation.operator() != operator.eq:
            raise ValueError("%s is not an equation with ==."
                             % equation)
        assert len(equation.operands()) == 2, \
            "%s is not an equation with two operands." % equation
        left_side, right_side = equation.operands()

        if left_side.operator() != function:
            raise ValueError("%s is not an evaluation of %s."
                             % (left_side, function))
        if len(left_side.operands()) != 1:
            raise ValueError("%s does not have one argument." %
                             (left_side,))

        try:
            polynomial_left = base_ring[var](left_side.operands()[0])
        except Exception:
            raise ValueError("%s is not a polynomial "
                             "in %s." % (left_side.operands()[0], var))
        if polynomial_left in base_ring and is_scalar(right_side):
            return {polynomial_left: to_list(right_side)}

        if polynomial_left.degree() != 1:
            raise ValueError("%s is not a polynomial of degree 1."
                             % (polynomial_left,))

        [r, base_power_K] = list(polynomial_left)
        try:
            K = log(base_power_K, base=base)
        except RuntimeError:
            K = 1
        try:
            K = K.simplify()
        except AttributeError:
            pass
        if K not in ZZ:
            raise ValueError("%s is not a power of %s."
                             % (base_power_K, base))
        if K < 1:
            raise ValueError("%d is less than %d."
                             % (base_power_K, base))

        from sage.symbolic.operators import add_vararg
        if right_side.operator() == add_vararg:
            function_calls = [o for o in right_side.operands()
                              if o.operator() == function]
            other_terms = [o for o in right_side.operands()
                           if o.operator() != function]
            if len(function_calls) != 1:
                raise ValueError(
                    "%s does not contain exactly one summand which "
                    "is an evaluation of %s."
                    % (right_side, function))
            next_function = function_calls[0]
            t = sum(other_terms)
            if not is_scalar(t):
                raise ValueError("%s contains %s."
                                 % (t, var))
        else:
            next_function = right_side
            t = 0

        if next_function.operator() != function:
            raise ValueError("%s is not an evaluation of %s."
                             % (next_function, function))
        if len(next_function.operands()) != 1:
            raise ValueError("%s does not have exactly one argument."
                             % (next_function,))

        try:
            polynomial_right = base_ring[var](next_function.operands()[0])
        except Exception:
            raise ValueError("%s is not a polynomial in %s."
                             % (next_function.operands()[0], var))
        if polynomial_right.degree() != 1:
            raise ValueError("%s is not a polynomial of degree 1."
                             % (polynomial_right,))
        [s, base_power_k] = list(polynomial_right)
        k = log(base_power_k, base=base)
        try:
            k = k.simplify()
        except AttributeError:
            pass
        if k not in ZZ:
            raise ValueError("%s is not a power of %s."
                             % (base_power_k, base))
        if k < 0:
            raise ValueError("%s is less than 1."
                             % (base_power_k,))
        if k >= K:
            raise ValueError("%d is greater or equal than %d."
                             % (base_power_k, base_power_K))

        parsed_equation = function(base**K * var + r) == \
            function(base**k * var + s) + t
        assert equation == parsed_equation, \
            "Parsing of %s failed for unknown reasons." % (equation,)

        rule = self.RecursionRule(K=K, r=r, k=k, s=s, t=to_list(t))
        return rule

    def Recursion(self, recursions, base, function=None, var=None,
                  input_alphabet=None, word_function=None,
                  is_zero=None, output_rings=[ZZ, QQ]):
        r"""
        Return a transducer realizing the given recursion when reading
        the digit expansion with base ``base``.

        INPUT:

        - ``recursions`` -- list or iterable of equations. Each
          equation has either the form

          - ``f(base^K * n + r) == f(base^k * n + s) + t`` for some
            integers ``0 <= k < K``, ``r`` and some ``t``---valid for
            all ``n`` such that the arguments on both sides are
            nonnegative---

          or the form

          - ``f(r) == t`` for some integer ``r`` and some ``t``.

          Alternatively, an equation may be replaced by a
          ``transducers.RecursionRule`` with the attributes ``K``,
          ``r``, ``k``, ``s``, ``t`` as above or a tuple ``(r, t)``.
          Note that ``t`` *must* be a list in this case.

        - ``base`` -- base of the digit expansion

        - ``function`` -- symbolic function ``f`` occurring in the
          recursions

        - ``var`` -- symbolic variable

        - ``input_alphabet`` -- (default: ``None``) a list of digits
          to be used as the input alphabet. If ``None`` and the base
          is an integer, ``input_alphabet`` is chosen to be
          ``srange(base.abs())``.

        - ``word_function`` -- (default: ``None``) a symbolic function.
          If not ``None``, ``word_function(arg1, ..., argn)`` in a symbolic
          recurrence relation is interpreted as a transition with output
          ``[arg1, ..., argn]``. This could not be entered in a symbolic
          recurrence relation because lists do not coerce into the
          :class:`~sage.symbolic.ring.SymbolicRing`.

        - ``is_zero`` -- (default: ``None``) a callable. The recursion
          relations are only well-posed if there is no cycle with
          nonzero output and input consisting of zeros. This parameter
          is used to determine whether the output of such a cycle is
          nonzero. By default, the output must evaluate to ``False`` as
          a boolean.

        - ``output_rings`` -- (default: ``[ZZ, QQ]``) a list of
          rings. The output labels are converted into the first ring of
          the list in which they are contained. If they are not
          contained in any ring, they remain in whatever ring they are
          after parsing the recursions, typically the symbolic ring.

        OUTPUT: a transducer ``T``

        The transducer is constructed such that ``T(expansion) == f(n)``
        if ``expansion`` is the digit expansion of ``n`` to the base
        ``base`` with the given input alphabet as set of digits. Here,
        the ``+`` on the right hand side of the recurrence relation is
        interpreted as the concatenation of words.

        The formal equations and initial conditions in the recursion
        have to be selected such that ``f`` is uniquely defined.

        EXAMPLES:

        -   The following example computes the Hamming weight of the
            ternary expansion of integers. ::

                sage: # needs sage.symbolic
                sage: function('f')
                f
                sage: var('n')
                n
                sage: T = transducers.Recursion([
                ....:     f(3*n + 1) == f(n) + 1,
                ....:     f(3*n + 2) == f(n) + 1,
                ....:     f(3*n) == f(n),
                ....:     f(0) == 0],
                ....:     3, f, n)
                sage: T.transitions()
                [Transition from (0, 0) to (0, 0): 0|-,
                 Transition from (0, 0) to (0, 0): 1|1,
                 Transition from (0, 0) to (0, 0): 2|1]

            To illustrate what this transducer does, we consider the
            example of `n=601`::

                sage: # needs sage.symbolic
                sage: ternary_expansion = 601.digits(base=3)
                sage: ternary_expansion
                [1, 2, 0, 1, 1, 2]
                sage: weight_sequence = T(ternary_expansion)
                sage: weight_sequence
                [1, 1, 1, 1, 1]
                sage: sum(weight_sequence)
                5

            Note that the digit zero does not show up in the output because
            the equation ``f(3*n) == f(n)`` means that no output is added to
            ``f(n)``.

        -   The following example computes the Hamming weight of the
            non-adjacent form, cf. the :wikipedia:`Non-adjacent_form`. ::

                sage: # needs sage.symbolic
                sage: function('f')
                f
                sage: var('n')
                n
                sage: T = transducers.Recursion([
                ....:     f(4*n + 1) == f(n) + 1,
                ....:     f(4*n - 1) == f(n) + 1,
                ....:     f(2*n) == f(n),
                ....:     f(0) == 0],
                ....:     2, f, n)
                sage: T.transitions()
                [Transition from (0, 0) to (0, 0): 0|-,
                 Transition from (0, 0) to (1, 1): 1|-,
                 Transition from (1, 1) to (0, 0): 0|1,
                 Transition from (1, 1) to (1, 0): 1|1,
                 Transition from (1, 0) to (1, 1): 0|-,
                 Transition from (1, 0) to (1, 0): 1|-]
                sage: [(s.label(), s.final_word_out)
                ....:  for s in T.iter_final_states()]
                [((0, 0), []),
                 ((1, 1), [1]),
                 ((1, 0), [1])]

            As we are interested in the weight only, we also output `1`
            for numbers congruent to `3` mod `4`. The actual expansion
            is computed in the next example.

            Consider the example of `29=(100\bar 101)_2` (as usual,
            the digit `-1` is denoted by `\bar 1` and digits are
            written from the most significant digit at the left to the
            least significant digit at the right; for the transducer,
            we have to give the digits in the reverse order)::

                sage: NAF = [1, 0, -1, 0, 0, 1]
                sage: ZZ(NAF, base=2)
                29
                sage: binary_expansion = 29.digits(base=2)
                sage: binary_expansion
                [1, 0, 1, 1, 1]
                sage: T(binary_expansion)                                               # needs sage.symbolic
                [1, 1, 1]
                sage: sum(T(binary_expansion))                                          # needs sage.symbolic
                3

            Indeed, the given non-adjacent form has three nonzero
            digits.

        -   The following example computes the non-adjacent form from the
            binary expansion, cf. the :wikipedia:`Non-adjacent_form`. In
            contrast to the previous example, we actually compute the
            expansion, not only the weight.

            We have to write the output `0` when converting an even number.
            This cannot be encoded directly by an equation in the symbolic
            ring, because ``f(2*n) == f(n) + 0`` would be equivalent to
            ``f(2*n) == f(n)`` and an empty output would be written.
            Therefore, we wrap the output in the symbolic function ``w``
            and use the parameter ``word_function`` to announce this.

            Similarly, we use ``w(-1, 0)`` to write an output word of
            length `2` in one iteration. Finally, we write ``f(0) == w()``
            to write an empty word upon completion.

            Moreover, there is a cycle with output ``[0]`` which---from
            the point of view of this method---is a contradicting recursion.
            We override this by the parameter ``is_zero``. ::

                sage: # needs sage.symbolic
                sage: var('n')
                n
                sage: function('f w')
                (f, w)
                sage: T = transducers.Recursion([
                ....:      f(2*n) == f(n) + w(0),
                ....:      f(4*n + 1) == f(n) + w(1, 0),
                ....:      f(4*n - 1) == f(n) + w(-1, 0),
                ....:      f(0) == w()],
                ....:      2, f, n,
                ....:      word_function=w,
                ....:      is_zero=lambda x: sum(x).is_zero())
                sage: T.transitions()
                [Transition from (0, 0) to (0, 0): 0|0,
                 Transition from (0, 0) to (1, 1): 1|-,
                 Transition from (1, 1) to (0, 0): 0|1,0,
                 Transition from (1, 1) to (1, 0): 1|-1,0,
                 Transition from (1, 0) to (1, 1): 0|-,
                 Transition from (1, 0) to (1, 0): 1|0]
                sage: for s in T.iter_states():
                ....:     print("{} {}".format(s, s.final_word_out))
                (0, 0) []
                (1, 1) [1, 0]
                (1, 0) [1, 0]

            We again consider the example of `n=29`::

                sage: T(29.digits(base=2))                                              # needs sage.symbolic
                [1, 0, -1, 0, 0, 1, 0]

            The same transducer can also be entered bypassing the
            symbolic equations::

                sage: R = transducers.RecursionRule
                sage: TR = transducers.Recursion([
                ....:       R(K=1, r=0, k=0, s=0, t=[0]),
                ....:       R(K=2, r=1, k=0, s=0, t=[1, 0]),
                ....:       R(K=2, r=-1, k=0, s=0, t=[-1, 0]),
                ....:       (0, [])],
                ....:       2,
                ....:       is_zero=lambda x: sum(x).is_zero())
                sage: TR == T                                                           # needs sage.symbolic
                True

        -   Here is an artificial example where some of the `s` are
            negative::

                sage: # needs sage.symbolic
                sage: function('f')
                f
                sage: var('n')
                n
                sage: T = transducers.Recursion([
                ....:     f(2*n + 1) == f(n-1) + 1,
                ....:     f(2*n) == f(n),
                ....:     f(1) == 1,
                ....:     f(0) == 0], 2, f, n)
                sage: T.transitions()
                [Transition from (0, 0) to (0, 0): 0|-,
                 Transition from (0, 0) to (1, 1): 1|-,
                 Transition from (1, 1) to (-1, 1): 0|1,
                 Transition from (1, 1) to (0, 0): 1|1,
                 Transition from (-1, 1) to (-1, 2): 0|-,
                 Transition from (-1, 1) to (1, 2): 1|-,
                 Transition from (-1, 2) to (-1, 1): 0|1,
                 Transition from (-1, 2) to (0, 0): 1|1,
                 Transition from (1, 2) to (-1, 2): 0|1,
                 Transition from (1, 2) to (1, 2): 1|1]
                sage: [(s.label(), s.final_word_out)
                ....:  for s in T.iter_final_states()]
                [((0, 0), []),
                 ((1, 1), [1]),
                 ((-1, 1), [0]),
                 ((-1, 2), [0]),
                 ((1, 2), [1])]

        -   Abelian complexity of the paperfolding sequence
            (cf. [HKP2015]_, Example 2.8)::

                sage: # needs sage.symbolic
                sage: T = transducers.Recursion([
                ....:     f(4*n) == f(2*n),
                ....:     f(4*n+2) == f(2*n+1)+1,
                ....:     f(16*n+1) == f(8*n+1),
                ....:     f(16*n+5) == f(4*n+1)+2,
                ....:     f(16*n+11) == f(4*n+3)+2,
                ....:     f(16*n+15) == f(2*n+2)+1,
                ....:     f(1) == 2, f(0) == 0]
                ....:     + [f(16*n+jj) == f(2*n+1)+2 for jj in [3,7,9,13]],
                ....:     2, f, n)
                sage: T.transitions()
                [Transition from (0, 0) to (0, 1): 0|-,
                 Transition from (0, 0) to (1, 1): 1|-,
                 Transition from (0, 1) to (0, 1): 0|-,
                 Transition from (0, 1) to (1, 1): 1|1,
                 Transition from (1, 1) to (1, 2): 0|-,
                 Transition from (1, 1) to (3, 2): 1|-,
                 Transition from (1, 2) to (1, 3): 0|-,
                 Transition from (1, 2) to (5, 3): 1|-,
                 Transition from (3, 2) to (3, 3): 0|-,
                 Transition from (3, 2) to (7, 3): 1|-,
                 Transition from (1, 3) to (1, 3): 0|-,
                 Transition from (1, 3) to (1, 1): 1|2,
                 Transition from (5, 3) to (1, 2): 0|2,
                 Transition from (5, 3) to (1, 1): 1|2,
                 Transition from (3, 3) to (1, 1): 0|2,
                 Transition from (3, 3) to (3, 2): 1|2,
                 Transition from (7, 3) to (1, 1): 0|2,
                 Transition from (7, 3) to (2, 1): 1|1,
                 Transition from (2, 1) to (1, 1): 0|1,
                 Transition from (2, 1) to (2, 1): 1|-]
                sage: for s in T.iter_states():
                ....:     print("{} {}".format(s, s.final_word_out))
                (0, 0) []
                (0, 1) []
                (1, 1) [2]
                (1, 2) [2]
                (3, 2) [2, 2]
                (1, 3) [2]
                (5, 3) [2, 2]
                (3, 3) [2, 2]
                (7, 3) [2, 2]
                (2, 1) [1, 2]
                sage: list(sum(T(n.bits())) for n in srange(1, 21))
                [2, 3, 4, 3, 4, 5, 4, 3, 4, 5, 6, 5, 4, 5, 4, 3, 4, 5, 6, 5]

        -   We now demonstrate the use of the ``output_rings``
            parameter.  If no ``output_rings`` are specified, the
            output labels are converted into ``ZZ``::

                sage: # needs sage.symbolic
                sage: function('f')
                f
                sage: var('n')
                n
                sage: T = transducers.Recursion([
                ....:     f(2*n + 1) == f(n) + 1,
                ....:     f(2*n) == f(n),
                ....:     f(0) == 2],
                ....:     2, f, n)
                sage: for t in T.transitions():
                ....:     print([x.parent() for x in t.word_out])
                []
                [Integer Ring]
                sage: [x.parent() for x in T.states()[0].final_word_out]
                [Integer Ring]

            In contrast, if ``output_rings`` is set to the empty list, the
            results are not converted::

                sage: T = transducers.Recursion([                                       # needs sage.symbolic
                ....:     f(2*n + 1) == f(n) + 1,
                ....:     f(2*n) == f(n),
                ....:     f(0) == 2],
                ....:     2, f, n, output_rings=[])
                sage: for t in T.transitions():                                         # needs sage.symbolic
                ....:     print([x.parent() for x in t.word_out])
                []
                [Symbolic Ring]
                sage: [x.parent() for x in T.states()[0].final_word_out]                # needs sage.symbolic
                [Symbolic Ring]

            Finally, we use a somewhat questionable conversion::

                sage: T = transducers.Recursion([                                       # needs sage.rings.finite_rings sage.symbolic
                ....:     f(2*n + 1) == f(n) + 1,
                ....:     f(2*n) == f(n),
                ....:     f(0) == 0],
                ....:     2, f, n, output_rings=[GF(5)])
                sage: for t in T.transitions():                                         # needs sage.rings.finite_rings sage.symbolic
                ....:     print([x.parent() for x in t.word_out])
                []
                [Finite Field of size 5]

        .. TODO::

            Extend the method to

            - non-integral bases,

            - higher dimensions.

        ALGORITHM:

        See [HKP2015]_, Section 6. However, there are also recursion
        transitions for states of level `<\kappa` if the recursion rules
        allow such a transition. Furthermore, the intermediate step of a
        non-deterministic transducer is left out by implicitly using
        recursion transitions. The well-posedness is checked in a
        truncated version of the recursion digraph.

        TESTS:

            The following tests fail due to missing or superfluous recursions
            or initial conditions. ::

                sage: var('n')                                                          # needs sage.symbolic
                n
                sage: function('f')                                                     # needs sage.symbolic
                f
                sage: transducers.Recursion([f(2*n) == f(n)],                           # needs sage.symbolic
                ....:     2, f, n)
                Traceback (most recent call last):
                ...
                ValueError: Missing recursions for input congruent to
                [1] modulo 2.

            ::

                sage: transducers.Recursion([f(2*n + 1) == f(n),                        # needs sage.symbolic
                ....:                        f(4*n) == f(2*n) + 1,
                ....:                        f(2*n) == f(n) + 1],
                ....:                       2, f, n)
                Traceback (most recent call last):
                ...
                ValueError: Conflicting rules congruent to 0 modulo 4.

            ::

                sage: transducers.Recursion([f(2*n + 1) == f(n) + 1,                    # needs sage.symbolic
                ....:                        f(2*n) == f(n),
                ....:                        f(0) == 0,
                ....:                        f(42) == 42], 2, f, n)
                Traceback (most recent call last):
                ...
                ValueError: Superfluous initial values for [42].

            ::

                sage: transducers.Recursion([f(2*n + 1) == f(n) + 1,                    # needs sage.symbolic
                ....:                        f(2*n) == f(n - 2) + 4,
                ....:                        f(0) == 0], 2, f, n)
                Traceback (most recent call last):
                ...
                ValueError: Missing initial values for [2].

            Here is an example of a transducer with a conflicting rule
            (it cannot hold for `n = 0`)::

                sage: T = transducers.Recursion([                                       # needs sage.symbolic
                ....:     f(2*n + 1) == f(n - 1),
                ....:     f(2*n) == f(n) + 1,
                ....:     f(1) == 1,
                ....:     f(0) == 0], 2, f, n)
                Traceback (most recent call last):
                ...
                ValueError: Conflicting recursion for [0].
        """
        from sage.graphs.digraph import DiGraph
        from sage.arith.srange import srange

        if is_zero is None:
            is_zero = lambda x: not x
        RuleRight = namedtuple('Rule', ['k', 's', 't'])
        initial_values = {}
        rules = []
        if input_alphabet is None and base in ZZ:
            input_alphabet = list(srange(base.abs()))

        for equation in recursions:
            if isinstance(equation, self.RecursionRule):
                rules.append(equation)
            elif isinstance(equation, tuple) and len(equation) == 2:
                initial_values[equation[0]] = equation[1]
            else:
                parsed = self._parse_recursion_equation_(
                    equation, base, function, var, word_function, output_rings)
                if isinstance(parsed, dict):
                    initial_values.update(parsed)
                elif isinstance(parsed, self.RecursionRule):
                    rules.append(parsed)
                else:
                    assert False

        max_K = max(rule.K for rule in rules)

        residues = [[None for r in range(base**k)]
                    for k in range(max_K + 1)]

        # Aim: residues[K][R] = RuleRight(k, s, t)
        # if and only if
        # f(base^K n + R) = f(base^k n + s) + t

        for given_rule in rules:
            q, remainder = given_rule.r.quo_rem(base**given_rule.K)
            rule = self.RecursionRule(K=given_rule.K,
                                      r=remainder,
                                      k=given_rule.k,
                                      s=given_rule.s - base**given_rule.k * q,
                                      t=given_rule.t)
            for m in range(max_K - rule.K + 1):
                for ell in range(base**m):
                    R = rule.r + base**rule.K * ell
                    if residues[rule.K + m][R] is not None:
                        raise ValueError(
                            "Conflicting rules congruent to %d modulo %d."
                            % (R, base**(rule.K + m)))
                    residues[rule.K + m][R] = RuleRight(k=rule.k + m,
                                                        s=rule.s + ell * base**rule.k,
                                                        t=rule.t)

        missing_residues = [R
                            for R, rule in enumerate(residues[max_K])
                            if rule is None]
        if missing_residues:
            raise ValueError("Missing recursions for input congruent "
                             "to %s modulo %s." % (missing_residues,
                                                   base**max_K))

        required_initial_values = set()

        def recursion_transition(carry, level, force_nonnegative_target):
            """
            Compute recursion transition leaving state ``(carry, level)``.

            INPUT:

            - ``carry`` -- integer

            - ``level`` -- integer

            - ``force_nonnegative_target`` -- boolean; if ``True``, only
              recursion transitions leading to a nonnegative carry are
              returned

            OUTPUT:

            A tuple ``((new_carry, new_level), output)`` if a recursion
            transition matching the specifications exists; otherwise,
            ``None``.
            """
            if level >= max_K:
                K = max_K
            else:
                K = level
            m, r = ZZ(carry).quo_rem(base**K)
            rule = residues[K][r]
            if rule is None:
                return None
            new_carry = base**rule.k * m + rule.s
            new_level = rule.k + level - K
            if new_carry + base**new_level < 0:
                return None
            if new_carry < 0 and force_nonnegative_target:
                return None
            if new_carry < 0 and carry >= 0:
                required_initial_values.add(carry)
            return (new_carry, new_level), rule.t

        def recursion_transitions(carry, level, force_nonnegative_target):
            """
            Compute the target and output of a maximal path of
            recursion transitions starting at state ``(carry, level)``.

            INPUT:

            - ``carry`` -- integer

            - ``level`` -- integer

            - ``force_nonnegative_target`` -- boolean; if ``True``, only
              recursion transitions leading to a nonnegative carry are
              allowed

            OUTPUT:

            A tuple ``((new_carry, new_level), output)``.
            """
            (c, j) = (carry, level)
            output = []
            while True:
                transition = recursion_transition(
                    c, j, force_nonnegative_target)
                if transition is None:
                    break
                (c, j) = transition[0]
                output += transition[1]

            return ((c, j), output)

        def transition_function(states2, input):
            (state_carry, state_level) = states2
            ((carry, level), output) = recursion_transitions(
                state_carry, state_level, False)
            # no more recursion transition is possible,
            # so this is now a storing transition
            carry += input * base**level
            level += 1
            # We now may proceed along recursion transitions
            # as long as the carries stay nonnegative.
            ((carry, level), new_output) = recursion_transitions(
                carry, level, True)
            return ((carry, level), output + new_output)

        T = Transducer(transition_function,
                       initial_states=[(0, 0)],
                       input_alphabet=input_alphabet)

        def edge_recursion_digraph(n):
            r"""
            Compute the list of outgoing edges of ``n`` in the recursion digraph.

            INPUT:

            - ``n`` -- integer

            OUTPUT:

            A list ``[(A(n), label)]`` if `A(n)<\infty`; otherwise an empty list.
            """
            m, r = ZZ(n).quo_rem(base**max_K)
            rule = residues[max_K][r]
            result = base**rule.k * m + rule.s
            if result >= 0:
                return [(result, rule.t)]
            else:
                return []

        def f(n):
            """
            Compute f(n) as defined by the recursion
            """
            if n in initial_values:
                return initial_values[n]
            [(m, offset)] = edge_recursion_digraph(n)
            return offset + f(m)

        carries = set(state.label()[0] for state in T.iter_states())

        recursion_digraph = DiGraph(
            {carry: dict(edge_recursion_digraph(carry))
             for carry in carries
             if carry >= 0},
            multiedges=False)

        initial_values_set = set(initial_values)

        missing_initial_values = required_initial_values.difference(
            initial_values_set)

        if missing_initial_values:
            raise ValueError(
                "Missing initial values for %s." %
                sorted(missing_initial_values))

        for cycle in recursion_digraph.all_simple_cycles():
            assert cycle[0] is cycle[-1]
            cycle_set = set(cycle)
            intersection = cycle_set.intersection(initial_values_set)
            if not intersection:
                raise ValueError(
                    "Missing initial condition for one of %s." %
                    cycle[1:])
            if len(intersection) > 1:
                raise ValueError(
                    "Too many initial conditions, only give one of %s." %
                    cycle[1:])
            required_initial_values.update(intersection)
            output_sum = sum([e[2]
                for e in recursion_digraph.outgoing_edge_iterator(cycle[1:])],
                             [])
            if not is_zero(output_sum):
                raise ValueError(
                    "Conflicting recursion for %s." %
                    cycle[1:])

        superfluous_initial_values = initial_values_set.difference(
            required_initial_values)

        if superfluous_initial_values:
            raise ValueError(
                "Superfluous initial values for %s." %
                sorted(superfluous_initial_values))

        for state in T.iter_states():
            state.is_final = True
            if state.label()[0] >= 0:
                state.final_word_out = f(state.label()[0])
            else:
                state.final_word_out = ZZ(0)

        return T


# Easy access to the automaton and transducer generators from the command line:
automata = AutomatonGenerators()
transducers = TransducerGenerators()
