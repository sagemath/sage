# sage.doctest: needs sage.modules sage.rings.finite_rings
r"""
Channels

Given an input space and an output space, a channel takes element from the
input space (the message) and transforms it into an element of the output space
(the transmitted message).

In Sage, Channels simulate error-prone transmission over communication
channels, and we borrow the nomenclature from communication theory, such as
"transmission" and "positions" as the elements of transmitted vectors.
Transmission can be achieved with two methods:

- :meth:`Channel.transmit`. Considering a channel ``Chan`` and a message
  ``msg``, transmitting ``msg`` with ``Chan`` can be done this way::

    Chan.transmit(msg)

  It can also be written in a more convenient way::

    Chan(msg)

- :meth:`transmit_unsafe`. This does the exact same thing as
  :meth:`transmit` except that it does not check if ``msg`` belongs to the
  input space of ``Chan``::

    Chan.transmit_unsafe(msg)

This is useful in e.g. an inner-loop of a long simulation as a
lighter-weight alternative to :meth:`Channel.transmit`.

This file contains the following elements:

    - :class:`Channel`, the abstract class for Channels
    - :class:`StaticErrorRateChannel`, which creates a specific number of errors in each
      transmitted message
    - :class:`ErrorErasureChannel`, which creates a specific number of errors and a
      specific number of erasures in each transmitted message
    - :class:`StaticRankErrorChannel`, which creates an error of specific rank in each
      transmitted message. This is the rank analogue of :class:`StaticErrorRateChannel`

AUTHORS:

 - David Lucas (2015): initial version
 - Maxime Bombar (2021): Rank metric channel
"""

# ****************************************************************************
#       Copyright (C) 2015 David Lucas <david.lucas@inria.fr>
#       Copyright (C) 2021 Maxime Bombar <maxime.bombar@crans.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from copy import copy

from sage.categories.cartesian_product import cartesian_product
from sage.categories.fields import Fields
from sage.functions.other import binomial
from sage.matrix.constructor import matrix
from sage.misc.abstract_method import abstract_method
from sage.misc.prandom import randint, random, sample
from sage.modules.free_module import VectorSpace, span
from sage.modules.free_module_element import vector
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.integer import Integer
from sage.structure.sage_object import SageObject


def random_error_vector(n, F, error_positions):
    r"""
    Return a vector of length ``n`` over ``F`` filled with random nonzero coefficients
    at the positions given by ``error_positions``.

    .. NOTE::

        This is a helper function, which should only be used when implementing new channels.

    INPUT:

    - ``n`` -- the length of the vector

    - ``F`` -- the field over which the vector is defined

    - ``error_positions`` -- the nonzero positions of the vector

    OUTPUT: a vector of ``F``

    AUTHORS:

    This function is taken from codinglib (https://bitbucket.org/jsrn/codinglib/)
    and was written by Johan Nielsen.

    EXAMPLES::

        sage: from sage.coding.channel import random_error_vector
        sage: random_error_vector(5, GF(2), [1,3])
        (0, 1, 0, 1, 0)
    """
    vect = [F.zero()]*n
    for i in error_positions:
        vect[i] = F._random_nonzero_element()
    return vector(F, vect)


def format_interval(t):
    r"""
    Return a formatted string representation of ``t``.

    This method should be called by any representation function in Channel classes.

    .. NOTE::

        This is a helper function, which should only be used when implementing new channels.

    INPUT:

    - ``t`` -- list or a tuple

    OUTPUT: string

    TESTS::

        sage: from sage.coding.channel import format_interval
        sage: t = (5, 5)
        sage: format_interval(t)
        '5'

        sage: t = (2, 10)
        sage: format_interval(t)
        'between 2 and 10'
    """
    return str(t[0]) if t[0] == t[1] else 'between %s and %s' % (t[0], t[1])


class Channel(SageObject):
    r"""
    Abstract top-class for Channel objects.

    All channel objects must inherit from this class. To implement a channel subclass, one should
    do the following:

    - inherit from this class,

    - call the super constructor,

    - override :meth:`transmit_unsafe`.

    While not being mandatory, it might be useful to reimplement representation methods (``_repr_`` and
    ``_latex_``).

    This abstract class provides the following parameters:

    - ``input_space`` -- the space of the words to transmit

    - ``output_space`` -- the space of the transmitted words
    """

    def __init__(self, input_space, output_space):
        r"""
        Initialize parameters for a Channel object.

        This is a private method, which should be called by the constructor
        of every encoder, as it automatically initializes the mandatory
        parameters of a Channel object.

        INPUT:

        - ``input_space`` -- the space of the words to transmit

        - ``output_space`` -- the space of the transmitted words

        EXAMPLES:

        We first create a new Channel subclass::

            sage: from sage.coding.channel import Channel
            sage: class ChannelExample(Channel):
            ....:   def __init__(self, input_space, output_space):
            ....:       super().__init__(input_space, output_space)

        We now create a member of our newly made class::

            sage: input = VectorSpace(GF(7), 6)
            sage: output = VectorSpace(GF(7), 5)
            sage: Chan = ChannelExample(input, output)

        We can check its parameters::

            sage: Chan.input_space()
            Vector space of dimension 6 over Finite Field of size 7
            sage: Chan.output_space()
            Vector space of dimension 5 over Finite Field of size 7
        """
        self._input_space = input_space
        self._output_space = output_space

    def transmit(self, message):
        r"""
        Return ``message``, modified accordingly with the algorithm of the channel it was
        transmitted through.

        Checks if ``message`` belongs to the input space, and returns an exception if not.
        Note that ``message`` itself is never modified by the channel.

        INPUT:

        - ``message`` -- a vector

        OUTPUT: a vector of the output space of ``self``

        EXAMPLES::

            sage: F = GF(59)^6
            sage: n_err = 2
            sage: Chan = channels.StaticErrorRateChannel(F, n_err)
            sage: msg = F((4, 8, 15, 16, 23, 42))
            sage: set_random_seed(10)
            sage: Chan.transmit(msg)
            (4, 8, 4, 16, 23, 53)

        We can check that the input ``msg`` is not modified::

            sage: msg
            (4, 8, 15, 16, 23, 42)

        If we transmit a vector which is not in the input space of ``self``::

            sage: n_err = 2
            sage: Chan = channels.StaticErrorRateChannel(GF(59)^6, n_err)
            sage: msg = (4, 8, 15, 16, 23, 42)
            sage: Chan.transmit(msg)
            Traceback (most recent call last):
            ...
            TypeError: Message must be an element of the input space for the given channel

        .. NOTE::

            One can also call directly ``Chan(message)``, which does the same as ``Chan.transmit(message)``
        """
        if message in self.input_space():
            return self.transmit_unsafe(message)
        else:
            raise TypeError("Message must be an element of the input space for the given channel")

    #Alias for transmit method
    __call__ = transmit

    def input_space(self):
        r"""
        Return the input space of ``self``.

        EXAMPLES::

            sage: n_err = 2
            sage: Chan = channels.StaticErrorRateChannel(GF(59)^6, n_err)
            sage: Chan.input_space()
            Vector space of dimension 6 over Finite Field of size 59
        """
        return self._input_space

    def output_space(self):
        r"""
        Return the output space of ``self``.

        EXAMPLES::

            sage: n_err = 2
            sage: Chan = channels.StaticErrorRateChannel(GF(59)^6, n_err)
            sage: Chan.output_space()
            Vector space of dimension 6 over Finite Field of size 59
        """
        return self._output_space

    @abstract_method
    def transmit_unsafe(self, message):
        r"""
        Return ``message``, modified accordingly with the algorithm of the channel it was
        transmitted through.

        This method does not check if ``message`` belongs to the input space of ``self``.

        This is an abstract method which should be reimplemented in all the subclasses of
        Channel.

        EXAMPLES::

            sage: n_err = 2
            sage: Chan = channels.StaticErrorRateChannel(GF(59)^6, n_err)
            sage: v = Chan.input_space().random_element()
            sage: Chan.transmit_unsafe(v)  # random
            (1, 33, 46, 18, 20, 49)
        """


class StaticErrorRateChannel(Channel):
    r"""
    Channel which adds a static number of errors to each message it transmits.

    The input space and the output space of this channel are the same.

    INPUT:

    - ``space`` -- the space of both input and output

    - ``number_errors`` -- the number of errors added to each transmitted message
      It can be either an integer of a tuple. If a tuple is passed as
      argument, the number of errors will be a random integer between the
      two bounds of the tuple.

    EXAMPLES:

    We construct a :class:`StaticErrorRateChannel` which adds 2 errors
    to any transmitted message::

        sage: n_err = 2
        sage: Chan = channels.StaticErrorRateChannel(GF(59)^40, n_err)
        sage: Chan
        Static error rate channel creating 2 errors, of input and output space
        Vector space of dimension 40 over Finite Field of size 59

    We can also pass a tuple for the number of errors::

        sage: n_err = (1, 10)
        sage: Chan = channels.StaticErrorRateChannel(GF(59)^40, n_err)
        sage: Chan
        Static error rate channel creating between 1 and 10 errors,
        of input and output space Vector space of dimension 40 over Finite Field of size 59
    """

    def __init__(self, space, number_errors):
        r"""
        TESTS:

        If the number of errors exceeds the dimension of the input space,
        it will return an error::

            sage: n_err = 42
            sage: Chan = channels.StaticErrorRateChannel(GF(59)^40, n_err)
            Traceback (most recent call last):
            ...
            ValueError: There might be more errors than the dimension of the input space
        """
        if isinstance(number_errors, (Integer, int)):
            number_errors = (number_errors, number_errors)
        if not isinstance(number_errors, (tuple, list)):
            raise ValueError("number_errors must be a tuple, a list, an Integer or a Python int")
        super().__init__(space, space)
        if number_errors[1] > space.dimension():
            raise ValueError("There might be more errors than the dimension of the input space")
        self._number_errors = number_errors

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: n_err = 42
            sage: Chan = channels.StaticErrorRateChannel(GF(59)^50, n_err)
            sage: Chan
            Static error rate channel creating 42 errors, of input and output space
            Vector space of dimension 50 over Finite Field of size 59
        """
        no_err = self.number_errors()
        return "Static error rate channel creating %s errors, of input and output space %s"\
                    % (format_interval(no_err), self.input_space())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: n_err = 42
            sage: Chan = channels.StaticErrorRateChannel(GF(59)^50, n_err)
            sage: latex(Chan)
            \textnormal{Static error rate channel creating 42 errors, of
            input and output space Vector space of dimension 50 over Finite Field of size 59}
        """
        no_err = self.number_errors()
        return "\\textnormal{Static error rate channel creating %s errors, of input and output space %s}"\
                % (format_interval(no_err), self.input_space())

    def transmit_unsafe(self, message):
        r"""
        Return ``message`` with as many errors as ``self._number_errors`` in it.

        If ``self._number_errors`` was passed as a tuple for the number of errors, it will
        pick a random integer between the bounds of the tuple and use it as the number of errors.

        This method does not check if ``message`` belongs to the input space of ``self``.

        INPUT:

        - ``message`` -- a vector

        OUTPUT: a vector of the output space

        EXAMPLES::

            sage: F = GF(59)^6
            sage: n_err = 2
            sage: Chan = channels.StaticErrorRateChannel(F, n_err)
            sage: msg = F((4, 8, 15, 16, 23, 42))
            sage: set_random_seed(10)
            sage: Chan.transmit_unsafe(msg)
            (4, 8, 4, 16, 23, 53)

        This checks that :issue:`19863` is fixed::

            sage: V = VectorSpace(GF(2), 1000)
            sage: Chan = channels.StaticErrorRateChannel(V, 367)
            sage: c = V.random_element()
            sage: (c - Chan(c)).hamming_weight()
            367
        """
        w = copy(message)
        number_errors = randint(*self.number_errors())
        V = self.input_space()
        R = V.base_ring()
        for i in sample(range(V.dimension()), number_errors):
            err = R.random_element()
            while (w[i] == err):
                err = R.random_element()
            w[i] = err
        return w

    def number_errors(self):
        r"""
        Return the number of errors created by ``self``.

        EXAMPLES::

            sage: n_err = 3
            sage: Chan = channels.StaticErrorRateChannel(GF(59)^6, n_err)
            sage: Chan.number_errors()
            (3, 3)
        """
        return self._number_errors


class ErrorErasureChannel(Channel):
    r"""
    Channel which adds errors and erases several positions in any message it transmits.

    The output space of this channel is a Cartesian product between its input
    space and a VectorSpace of the same dimension over `\GF{2}`.

    INPUT:

    - ``space`` -- the input and output space

    - ``number_errors`` -- the number of errors created in each transmitted
      message. It can be either an integer of a tuple. If a tuple is passed as
      an argument, the number of errors will be a random integer between the
      two bounds of this tuple.

    - ``number_erasures`` -- the number of erasures created in each transmitted
      message. It can be either an integer of a tuple. If a tuple is passed as an
      argument, the number of erasures will be a random integer between the
      two bounds of this tuple.

    EXAMPLES:

    We construct a ErrorErasureChannel which adds 2 errors
    and 2 erasures to any transmitted message::

        sage: n_err, n_era = 2, 2
        sage: Chan = channels.ErrorErasureChannel(GF(59)^40, n_err, n_era)
        sage: Chan
        Error-and-erasure channel creating 2 errors and 2 erasures
        of input space Vector space of dimension 40 over Finite Field of size 59
        and output space The Cartesian product of (Vector space of dimension 40
        over Finite Field of size 59, Vector space of dimension 40 over Finite Field of size 2)

    We can also pass the number of errors and erasures as a couple of integers::

        sage: n_err, n_era = (1, 10), (1, 10)
        sage: Chan = channels.ErrorErasureChannel(GF(59)^40, n_err, n_era)
        sage: Chan
        Error-and-erasure channel creating between 1 and 10 errors and
        between 1 and 10 erasures of input space Vector space of dimension 40
        over Finite Field of size 59 and output space The Cartesian product of
        (Vector space of dimension 40 over Finite Field of size 59,
        Vector space of dimension 40 over Finite Field of size 2)
    """

    def __init__(self, space, number_errors, number_erasures):
        r"""
        TESTS:

        If the sum of number of errors and number of erasures
        exceeds (or may exceed, in the case of tuples) the dimension of the input space,
        it will return an error::

            sage: n_err, n_era = 21, 21
            sage: Chan = channels.ErrorErasureChannel(GF(59)^40, n_err, n_era)
            Traceback (most recent call last):
            ...
            ValueError: The total number of errors and erasures cannot exceed the dimension of the input space
        """
        if isinstance(number_errors, (Integer, int)):
            number_errors = (number_errors, number_errors)
        if not isinstance(number_errors, (tuple, list)):
            raise ValueError("number_errors must be a tuple, a list, an Integer or a Python int")

        if isinstance(number_erasures, (Integer, int)):
            number_erasures = (number_erasures, number_erasures)
        if not isinstance(number_erasures, (tuple, list)):
            raise ValueError("number_erasures must be a tuple, a list, an Integer or a Python int")

        output_space = cartesian_product([space, VectorSpace(GF(2), space.dimension())])
        super().__init__(space, output_space)
        if number_errors[1] + number_erasures[1] > space.dimension():
            raise ValueError("The total number of errors and erasures cannot exceed the dimension of the input space")
        self._number_errors = number_errors
        self._number_erasures = number_erasures

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: n_err, n_era = 21, 21
            sage: Chan = channels.ErrorErasureChannel(GF(59)^50, n_err, n_era)
            sage: Chan
            Error-and-erasure channel creating 21 errors and 21 erasures
            of input space Vector space of dimension 50 over Finite Field of size 59
            and output space The Cartesian product of (Vector space of dimension 50
            over Finite Field of size 59, Vector space of dimension 50 over Finite Field of size 2)
        """
        no_err = self.number_errors()
        no_era = self.number_erasures()
        return "Error-and-erasure channel creating %s errors and %s erasures of input space %s and output space %s"\
                % (format_interval(no_err), format_interval(no_era), self.input_space(), self.output_space())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: n_err, n_era = 21, 21
            sage: Chan = channels.ErrorErasureChannel(GF(59)^50, n_err, n_era)
            sage: latex(Chan)
            \textnormal{Error-and-erasure channel creating 21 errors and 21 erasures
            of input space Vector space of dimension 50 over Finite Field of size 59
            and output space The Cartesian product of (Vector space of dimension 50
            over Finite Field of size 59, Vector space of dimension 50 over Finite Field of size 2)}
        """
        no_err = self.number_errors()
        no_era = self.number_erasures()
        return "\\textnormal{Error-and-erasure channel creating %s errors and %s erasures of input space %s and output space %s}"\
                % (format_interval(no_err), format_interval(no_era), self.input_space(), self.output_space())

    def transmit_unsafe(self, message):
        r"""
        Return ``message`` with as many errors as ``self._number_errors``
        in it, and as many erasures as ``self._number_erasures`` in it.

        If ``self._number_errors`` was passed as a tuple for the number of errors, it will
        pick a random integer between the bounds of the tuple and use it as the number of errors.
        It does the same with ``self._number_erasures``.

        All erased positions are set to 0 in the transmitted message.
        It is guaranteed that the erasures and the errors will never overlap:
        the received message will always contains exactly as many errors and erasures
        as expected.

        This method does not check if ``message`` belongs to the input space of ``self``.

        INPUT:

        - ``message`` -- a vector

        OUTPUT: a couple of vectors, namely:

          - the transmitted message, which is ``message`` with erroneous and
            erased positions
          - the erasure vector, which contains ``1`` at the erased positions of
            the transmitted message and ``0`` elsewhere.

        EXAMPLES::

            sage: F = GF(59)^11
            sage: n_err, n_era = 2, 2
            sage: Chan = channels.ErrorErasureChannel(F, n_err, n_era)
            sage: msg = F((3, 14, 15, 9, 26, 53, 58, 9, 7, 9, 3))
            sage: set_random_seed(10)
            sage: Chan.transmit_unsafe(msg)
            ((31, 0, 15, 9, 38, 53, 58, 9, 0, 9, 3), (0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0))
        """
        number_errors = randint(*self.number_errors())
        number_erasures = randint(*self.number_erasures())
        V = self.input_space()
        n = V.dimension()
        zero = V.base_ring().zero()

        errors = sample(range(n), number_errors + number_erasures)
        error_positions = errors[:number_errors]
        erasure_positions = errors[number_errors:]

        error_vector = random_error_vector(n, V.base_ring(), error_positions)
        erasure_vector = random_error_vector(n , GF(2), erasure_positions)

        message = message + error_vector

        for i in erasure_positions:
            message[i] = zero
        return message, erasure_vector

    def number_errors(self):
        r"""
        Return the number of errors created by ``self``.

        EXAMPLES::

            sage: n_err, n_era = 3, 0
            sage: Chan = channels.ErrorErasureChannel(GF(59)^6, n_err, n_era)
            sage: Chan.number_errors()
            (3, 3)
        """
        return self._number_errors

    def number_erasures(self):
        r"""
        Return the number of erasures created by ``self``.

        EXAMPLES::

            sage: n_err, n_era = 0, 3
            sage: Chan = channels.ErrorErasureChannel(GF(59)^6, n_err, n_era)
            sage: Chan.number_erasures()
            (3, 3)
        """
        return self._number_erasures


class QarySymmetricChannel(Channel):
    r"""
    The `q`-ary symmetric, memoryless communication channel.

    Given an alphabet `\Sigma` with `|\Sigma| = q` and an error probability
    `\epsilon`, a `q`-ary symmetric channel sends an element of `\Sigma` into the
    same element with probability `1 - \epsilon`, and any one of the other `q -
    1` elements with probability `\frac{\epsilon}{q - 1}`. This implementation
    operates over vectors in `\Sigma^n`, and "transmits" each element of the
    vector independently in the above manner.

    Though `\Sigma` is usually taken to be a finite field, this implementation
    allows any structure for which Sage can represent `\Sigma^n` and for which
    `\Sigma` has a ``random_element()`` method. However, beware that if `\Sigma`
    is infinite, errors will not be uniformly distributed (since
    ``random_element()`` does not draw uniformly at random).

    The input space and the output space of this channel are the same:
    `\Sigma^n`.

    INPUT:

    - ``space`` -- the input and output space of the channel; it has to be
      `\GF{q}^n` for some finite field `\GF{q}`

    - ``epsilon`` -- the transmission error probability of the individual elements

    EXAMPLES:

    We construct a :class:`QarySymmetricChannel` which corrupts 30% of all
    transmitted symbols::

        sage: epsilon = 0.3
        sage: Chan = channels.QarySymmetricChannel(GF(59)^50, epsilon)
        sage: Chan
        q-ary symmetric channel with error probability 0.300000000000000,
         of input and output space
          Vector space of dimension 50 over Finite Field of size 59
    """

    def __init__(self, space, epsilon):
        r"""
        TESTS:

        If ``space`` is not a vector space, an error is raised::

            sage: epsilon = 0.42
            sage: Chan = channels.QarySymmetricChannel(GF(59), epsilon)
            Traceback (most recent call last):
            ...
            ValueError: space has to be of the form Sigma^n, where Sigma has a random_element() method

        If ``epsilon`` is not between 0 and 1, an error is raised::

            sage: epsilon = 42
            sage: Chan = channels.QarySymmetricChannel(GF(59)^50, epsilon)
            Traceback (most recent call last):
            ...
            ValueError: Error probability must be between 0 and 1
        """
        if epsilon >= 1 or epsilon <= 0:
            raise ValueError("Error probability must be between 0 and 1")

        super().__init__(space, space)
        self._epsilon = epsilon
        try:
            self.transmit_unsafe(space.random_element())
        except Exception:
            raise ValueError("space has to be of the form Sigma^n, where Sigma has a random_element() method")

    def __repr__(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: epsilon = 0.3
            sage: Chan = channels.QarySymmetricChannel(GF(59)^50, epsilon)
            sage: Chan
            q-ary symmetric channel with error probability 0.300000000000000,
            of input and output space Vector space of dimension 50 over Finite Field of size 59
        """
        return "q-ary symmetric channel with error probability %s, of input and output space %s"\
                    % (self.error_probability(), self.input_space())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: epsilon = 0.3
            sage: Chan = channels.QarySymmetricChannel(GF(59)^50, epsilon)
            sage: latex(Chan)
            \textnormal{q-ary symmetric channel with error probability 0.300000000000000,
            of input and output space Vector space of dimension 50 over Finite Field of size 59}
        """
        return "\\textnormal{q-ary symmetric channel with error probability %s, of input and output space %s}"\
                    % (self.error_probability(), self.input_space())

    def transmit_unsafe(self, message):
        r"""
        Return ``message`` where each of the symbols has been changed to another from the alphabet with
        probability :meth:`error_probability`.

        This method does not check if ``message`` belongs to the input space of ``self``.

        INPUT:

        - ``message`` -- a vector

        EXAMPLES::

            sage: F = GF(59)^11
            sage: epsilon = 0.3
            sage: Chan = channels.QarySymmetricChannel(F, epsilon)
            sage: msg = F((3, 14, 15, 9, 26, 53, 58, 9, 7, 9, 3))
            sage: set_random_seed(10)
            sage: Chan.transmit_unsafe(msg)
            (3, 14, 15, 53, 12, 53, 58, 9, 55, 9, 3)
        """
        epsilon = self.error_probability()
        V = self.input_space()
        F = V.base_ring()
        msg = copy(message.list())
        for i in range(len(msg)):
            if random() <= epsilon:
                a = F.random_element()
                while a == msg[i]:
                    a = F.random_element()
                msg[i] = a
        return V(msg)

    def error_probability(self):
        r"""
        Return the error probability of a single symbol transmission of
        ``self``.

        EXAMPLES::

            sage: epsilon = 0.3
            sage: Chan = channels.QarySymmetricChannel(GF(59)^50, epsilon)
            sage: Chan.error_probability()
            0.300000000000000
        """
        return self._epsilon

    def probability_of_exactly_t_errors(self, t):
        r"""
        Return the probability ``self`` has to return
        exactly ``t`` errors.

        INPUT:

        - ``t`` -- integer

        EXAMPLES::

            sage: epsilon = 0.3
            sage: Chan = channels.QarySymmetricChannel(GF(59)^50, epsilon)
            sage: Chan.probability_of_exactly_t_errors(15)
            0.122346861835401
        """
        n = self.input_space().dimension()
        epsilon = self.error_probability()
        return binomial(n, t) * epsilon**t * (1-epsilon)**(n-t)

    def probability_of_at_most_t_errors(self, t):
        r"""
        Return the probability ``self`` has to return
        at most ``t`` errors.

        INPUT:

        - ``t`` -- integer

        EXAMPLES::

            sage: epsilon = 0.3
            sage: Chan = channels.QarySymmetricChannel(GF(59)^50, epsilon)
            sage: Chan.probability_of_at_most_t_errors(20)
            0.952236164579467
        """
        return sum(self.probability_of_exactly_t_errors(i)
                for i in range(t+1))


class StaticRankErrorChannel(Channel):
    r"""
    Channel which adds an error of static rank to each message it transmits.

    The input space and the output space of this channel are the same.

    INPUT:

    - ``space`` -- the space of both input and output

    - ``rank_errors`` -- the rank of the error added to each transmitted message
      It can be either an integer of a tuple. If a tuple is passed as
      argument, the rank of the error will be a random integer between the
      two bounds of the tuple.

    - ``relative_field`` -- The field to which the extension is relative.
      If not given, it will default to the prime_subfield of the ambient space
      base_field.

    EXAMPLES:

    We construct a StaticRankErrorChannel which adds error of rank 2
    to any transmitted message::

        sage: n_err = 2
        sage: Chan = channels.StaticRankErrorChannel(GF(256)^6, n_err)
        sage: Chan
        Channel creating error of rank 2 over Finite Field of size 2,
        of input and output space Vector space of dimension 6
        over Finite Field in z8 of size 2^8

    We can also pass a tuple for the number of errors::

        sage: n_err = (1, 4)
        sage: Chan = channels.StaticRankErrorChannel(GF(256)^6, n_err)
        sage: Chan
        Channel creating error of rank between 1 and 4 over Finite Field
        of size 2, of input and output space Vector space of dimension 6
        over Finite Field in z8 of size 2^8
    """

    def __init__(self, space, rank_errors, relative_field=None):
        r"""
        TESTS:

        If the number of errors exceeds the dimension of the input space,
        or the field extension degree, then an error is raised::

            sage: n_err = 7
            sage: channels.StaticRankErrorChannel(GF(256)^6, n_err)
            Traceback (most recent call last):
            ...
            ValueError: There might be errors of rank larger than the dimension of the input space.
            sage: channels.StaticRankErrorChannel(GF(64)^8, n_err)
            Traceback (most recent call last):
            ...
            ValueError: There might be errors of rank larger than the field extension degree.

        If ``relative_field`` is specified and is not a subfield of the base field,
        it will return an error::

            sage: n_err = 2
            sage: Chan = channels.StaticRankErrorChannel(GF(16)^6, n_err, GF(8))
            Traceback (most recent call last):
            ...
            ValueError: Finite Field in z4 of size 2^4 is not an extension of
            Finite Field in z3 of size 2^3

        If ``relative_field`` is specified and is not a field,
        it will return an error::

            sage: n_err = 2
            sage: Chan = channels.StaticRankErrorChannel(GF(16)^6, n_err, GF(4)^2)
            Traceback (most recent call last):
            ...
            ValueError: relative_field must be a Field and
            Vector space of dimension 2 over Finite Field in z2
            of size 2^2 is not.
        """
        if isinstance(rank_errors, (Integer, int)):
            rank_errors = (rank_errors, rank_errors)
        if not isinstance(rank_errors, (tuple, list)):
            raise ValueError("rank_errors must be a tuple, a list, an Integer or a Python int")
        super(StaticRankErrorChannel, self).__init__(space, space)
        if rank_errors[1] > space.dimension():
            raise ValueError("There might be errors of rank larger than the dimension of the input space.")
        self._rank_errors = rank_errors
        self._base_field = space.base_field()
        if not relative_field:
            self._relative_field = self._base_field.prime_subfield()
        else:
            if relative_field not in Fields():
                raise ValueError("relative_field must be a Field and %s is not." % relative_field)
            if not relative_field.is_subring(self._base_field):
                raise ValueError("%s is not an extension of %s" % (self._base_field, relative_field))
            self._relative_field = relative_field
        V, vec_to_field, field_to_vec = self._base_field.vector_space(self._relative_field)
        self._column_space = V
        self._extension_degree = V.dimension()
        if rank_errors[1] > self._extension_degree:
            raise ValueError("There might be errors of rank larger than the field extension degree.")
        self._map_to_field = vec_to_field

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: n_err = 2
            sage: Chan = channels.StaticRankErrorChannel(GF(256)^6, n_err)
            sage: Chan
            Channel creating error of rank 2 over Finite Field of size 2,
            of input and output space Vector space of dimension 6
            over Finite Field in z8 of size 2^8
        """
        no_err = self.rank_errors()
        return "Channel creating error of rank %s over %s, of input and output space %s"\
                    % (format_interval(no_err), self._relative_field, self.input_space())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: n_err = 2
            sage: Chan = channels.StaticRankErrorChannel(GF(256)^6, n_err)
            sage: latex(Chan)
            \textnormal{Channel creating error of rank 2 over Finite Field
            of size 2, of input and output space Vector space of dimension 6
            over Finite Field in z8 of size 2^8}
        """
        no_err = self.rank_errors()
        return "\\textnormal{Channel creating error of rank %s over %s, of input and output space %s}"\
                % (format_interval(no_err), self._relative_field, self.input_space())

    def transmit_unsafe(self, message):
        r"""
        Return ``message`` with as many errors as ``self._rank_errors`` in it.

        If ``self._rank_errors`` was passed as a tuple for the number of errors, it will
        pick a random integer between the bounds of the tuple and use it as the number of errors.

        This method does not check if ``message`` belongs to the input space of``self``.

        INPUT:

        - ``message`` -- a vector

        OUTPUT:

        - a vector of the output space

        EXAMPLES::

            sage: F = GF(16)^6
            sage: n_err = 2
            sage: Chan = channels.StaticRankErrorChannel(F, n_err, GF(4))
            sage: set_random_seed(10)
            sage: msg = F.random_element()
            sage: msg
            (z4 + 1, z4, z4^3 + z4 + 1, z4^3 + z4^2 + z4 + 1, z4^2, z4^2)
            sage: c = Chan.transmit_unsafe(msg)

        TESTS::

            sage: from sage.coding.linear_rank_metric import rank_distance
            sage: rank_distance(msg, c, GF(4))
            2
        """
        rank_errors = randint(*self.rank_errors())
        Fq = self._relative_field
        V = self._column_space
        n = self.input_space().dimension()
        good = False
        while not good:
            basis = [V.random_element() for i in range(rank_errors)]
            R = span(basis)
            err = [R.random_element() for i in range(n)]
            M = matrix(Fq, err)
            good = (M.rank() == rank_errors)
        e = vector(map(self._map_to_field, M.rows()))
        w = message + e
        return w

    def rank_errors(self):
        r"""
        Return the number of rank errors created by ``self``.

        EXAMPLES::

            sage: n_err = 3
            sage: Chan = channels.StaticRankErrorChannel(GF(256)^6, n_err)
            sage: Chan.rank_errors()
            (3, 3)
        """
        return self._rank_errors

    def number_errors(self):
        r"""
        Return the number of errors created by ``self``.

        .. NOTE::

            This function is here to ease the life of people coming from the Hamming world.

        EXAMPLES::

            sage: n_err = 3
            sage: Chan = channels.StaticRankErrorChannel(GF(256)^6, n_err)
            sage: Chan.number_errors()
            (3, 3)
        """
        return self._rank_errors
