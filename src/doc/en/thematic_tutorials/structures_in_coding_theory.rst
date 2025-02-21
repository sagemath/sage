.. -*- coding: utf-8 -*-
.. _structures_in_coding_theory:

===============================================
How to write your own classes for coding theory
===============================================

.. MODULEAUTHOR:: David Lucas (2016) for initial version. Marketa Slukova (2019) for this version.

This tutorial, designed for advanced users who want to build their own classes,
will explain step by step what you need to do to write code which integrates
well in the framework of coding theory.
During this tutorial, we will cover the following parts:

- how to write a new **code family**
- how to write a new **encoder**
- how to write a new **decoder**
- how to write a new **channel**

This tutorial focuses on the example of implementation of the repetition code.
At the end of each part, we will summarize every important step of the
implementation. If one just wants a quick access to the implementation of one of
the objects cited above, one can jump directly to the end of related part,
which presents a summary of what to do.

.. contents:: Table of contents
   :depth: 2
   :class: this-will-duplicate-information-and-it-is-still-useful-here

I. Abstract classes
===================

There is a number of abstract classes representing different types of codes in
the coding module of Sage. Depending on the type of code you want to implement,
you should inherit from different abstract classes.

The most generic class is :class:`sage.coding.abstract_code:AbstractCode`. This
class makes no assumptions about linearity, metric, finiteness or the number of
alphabets. The abstract notion of "code" that is implicitly used for this class
is any enumerable subset of a cartesian product
`A_1 \times A_2 \times \ldots \times A_n` for some sets `A_i`.

If, for example, you want to create a non-linear code family, this is the class
you should inherit from.

We give a small example of creating a non-linear code family. Take the code
consisting of the following `l` words:

`\{00 \ldots 00, 10 \ldots 00, 11 \ldots 00, \ldots , 11 \ldots 10, 11 \ldots 11 \}`.

Here is how we can implement it::

    sage: from sage.coding.abstract_code import AbstractCode
    sage: class ExampleCodeFamily(AbstractCode):
    ....:     def __init__(self, length):
    ....:         super(ExampleCodeFamily, self).__init__(length)
    ....:     def __iter__(self):
    ....:         for i in range(self.length() + 1):
    ....:             yield vector([1 for j in range(i)] + [0 for k in range(i, self.length())])
    ....:     def __contains__(self, word):
    ....:         return word in list(self)
    ....:     def _repr_(self):
    ....:         return "Dummy code of length {}".format(self.length())

    sage: C = ExampleCodeFamily(4) # check that this works.
    sage: C.list()
    [(0, 0, 0, 0), (1, 0, 0, 0), (1, 1, 0, 0), (1, 1, 1, 0), (1, 1, 1, 1)]

Focusing on linear codes, the most generic representative is the class
:class:`sage.coding.linear_code_no_metric:AbstractLinearCodeNoMetric` which
contains all the methods that all linear codes share regardless of what their
metric is. If you want to implement a linear code over some metric which has
not been implemented in Sage yet, this is the class you should be inheriting
from.

We have two metric specific abstract linear code classes,
:class:`sage.coding.linear_code:AbstractLinearCode` for Hamming metric and
:class:`sage.coding.linear_rank_metric:AbstractLinearRankMetricCode` for rank
metric. If you wish to implement a linear code class over one of these metrics,
you should inherit from the given abstract class.

II. The repetition code
=======================

We want to implement in Sage the well-known repetition code.
Its definition follows:

the :math:`(n, 1)`-repetition code over :math:`\GF{q}` is the code formed
by all the vectors of :math:`\GF{q}^{n}` of the form
:math:`(i, i, i, \dots, i)` for all :math:`i \in \GF{q}`.

For example, the :math:`(3, 1)`-repetition code over :math:`\GF{2}` is:
:math:`C = \{(0, 0, 0), (1, 1, 1)\}`.

The encoding is very simple, it only consists in repeating :math:`n`
times the input symbol and pick the vector thus formed.

The decoding uses majority voting to select the right symbol
(over :math:`\GF{2}`). If we receive the word :math:`(1, 0, 1)`
(example cont'd), we deduce that the original word was :math:`(1)`.
It can correct up to :math:`\left\lceil \frac{n-1}{2} \right\rceil` errors.

Through all this tutorial, we will illustrate the implementation of the
:math:`(n, 1)`-repetition code over :math:`\GF{2}`.

III. Write a new code class
===========================

The first thing to do to write a new code class is to identify which abstract
class to inherit from. Since the repetition code is linear and we take it over
the Hamming metric, this means that we will inherit from
:class:`sage.coding.linear_code.AbstractLinearCode`.

Now we have to identify the initializing parameters of this class, which are:

- the length of the code,
- the base field of the code,
- the default encoder for the code,
- the default decoder for the code and
- any other useful argument we want to set at construction time.

For our code, we know its length, its dimension, its base field, one encoder
and one decoder.

Now we isolated the parameters of the code, we can write the
constructor of our class.
As we said, every linear code class over the Hamming metric must inherit from
:class:`sage.coding.linear_code.AbstractLinearCode`.
This class provides a lot of useful methods and, as we illustrate thereafter,
a default constructor which sets the *length*, the *base field*,
the *default encoder* and the *default decoder* as class parameters.
We also need to create the dictionary of known encoders and decoders
for the class.

Let us now write the constructor for our code class,
that we store in some file called ``repetition_code.py``::

    sage: from sage.coding.linear_code import AbstractLinearCode
    sage: from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
    sage: class BinaryRepetitionCode(AbstractLinearCode):
    ....:     _registered_encoders = {}
    ....:     _registered_decoders = {}
    ....:     def __init__(self, length):
    ....:         super(BinaryRepetitionCode, self).__init__(GF(2), length,
    ....:           "RepetitionGeneratorMatrixEncoder", "MajorityVoteDecoder")
    ....:         self._dimension = 1

As you notice, the constructor is really simple. Most of the work is indeed
managed by the topclass through the ``super`` statement.
Note that the dimension is not set by the abstract class, because for some
code families the exact dimension is hard to compute.
If the exact dimension is known, set it using ``_dimension``
as a class parameter.

We can now write representation methods for our code class::

    sage: def _repr_(self):
    ....:     return "Binary repetition code of length %s" % self.length()
    sage: def _latex_(self):
    ....:     return "\textnormal{Binary repetition code of length } %s" % self.length()

We also write a method to check equality::

    sage: def __eq__(self, other):
    ....:     return (isinstance(other, BinaryRepetitionCode)
    ....:             and self.length() == other.length()
    ....:             and self.dimension() == other.dimension())

After these examples, you probably noticed that we use two methods,
namely ``length()`` and ``dimension()`` without defining them.
That is because their implementation is provided in parent classes of
:class:`sage.coding.linear_code.AbstractLinearCode`, which are
:class:`sage.coding.linear_code_no_metric.AbstractLinearCodeNoMetric` and
:class:`sage.coding.abstract_code.AbstractCode`

They provide a default implementation of the following getter methods:

- :meth:`sage.coding.abstract_code.AbstractCode.length`,
- :meth:`sage.coding.linear_code_no_metric.AbstractLinearCodeNoMetric.dimension`,
- :meth:`sage.coding.linear_code_no_metric.AbstractLinearCodeNoMetric.base_field` and
- :meth:`sage.coding.linear_code_no_metric.AbstractLinearCodeNoMetric.ambient_space`.

They also provide several other useful methods, such as ``__contains__``. Note
that a lot of these other methods rely on the computation of a generator matrix.
It is thus highly recommended to set an encoder which knows how to compute such
a matrix as default encoder. As default encoder will be used by all these
methods which expect a generator matrix, if one provides a default encoder which
does not have a ``generator_matrix`` method, a lot of generic methods will fail.

As our code family is really simple, we do not need anything else,
and the code provided above is enough to describe properly a repetition code.

Summary of the implementation for linear codes over the Hamming metric
----------------------------------------------------------------------

1. Inherit from :class:`sage.coding.linear_code.AbstractLinearCode`.
2. Add ``_registered_encoders =  {}`` and ``_registered_decoders = {}``
   as class variables.
3. Add this line in the class' constructor::

      super(ClassName, self).__init__(base_field, length, "DefaultEncoder", "DefaultDecoder")

4. Implement representation methods (not mandatory, but highly advised)
   ``_repr_`` and ``_latex_``.
5. Implement ``__eq__``.
6. ``__ne__``, ``length`` and ``dimension`` come with the abstract class.

Please note that ``dimension`` will not work is there is no field
``_dimension`` as class parameter.

We now know how to write a new code class.
Let us see how to write a new encoder and a new decoder.


IV. Write a new encoder class
=============================

Let us continue our example. We ask the same question as before:
what do we need to describe the encoder?
For most of the cases (this one included), we only need the associated code.
In that case, writing the constructor is really straightforward
(we store the code in the same ``.py`` file as the code class)::

    sage: from sage.coding.encoder import Encoder
    sage: class BinaryRepetitionCodeGeneratorMatrixEncoder(Encoder):
    ....:     def __init__(self, code):
    ....:         super(BinaryRepetitionCodeGeneratorMatrixEncoder, self).__init__(code)

Same thing as before, as an encoder always needs to know its associated code,
the work can be done by the base class.
Remember to inherit from :class:`sage.coding.encoder.Encoder`!

We also want to override representation methods ``_repr_`` and ``_latex_``::

    sage: def _repr_(self):
    ....:     return "Binary repetition encoder for the %s" % self.code()
    sage: def _latex_(self):
    ....:     return "\textnormal{Binary repetition encoder for the } %s" % self.code()

And we want to have an equality check too::

    sage: def __eq__(self, other):
    ....:     return (isinstance(other, BinaryRepetitionCodeGeneratorMatrixEncoder)
    ....:             and self.code() == other.code())

As before, default getter method is provided by the topclass,
namely :meth:`sage.coding.encoder.Encoder.code`.

All we have to do is to implement the methods related to the encoding.
This implementation changes quite a lot whether
we have a generator matrix or not.

We have a generator matrix
--------------------------

In that case, the message space is a vector space, and it is especially easy:
the only method you need to implement is ``generator_matrix``.

Continuing our example, it will be::

    sage: def generator_matrix(self):
    ....:     n = self.code().length()
    ....:     return Matrix(GF(2), 1, n, [GF(2).one()] * n)

As the topclass provides default implementation for encode and the inverse
operation, that we call *unencode*
(see: :meth:`sage.coding.encoder.Encoder.encode` and
:meth:`sage.coding.encoder.Encoder.unencode`), alongside
with a default implementation of
:meth:`sage.coding.encoder.Encoder.message_space`, our work here is done.

.. NOTE::

    Default ``encode`` method multiplies the provide word by the generator
    matrix, while default ``unencode`` computes an information set for
    the generator matrix, inverses it and performs a matrix-vector
    multiplication to recover the original message.
    If one has a better implementation for one's specific code family,
    one should obviously override the default ``encode`` and ``unencode``.

We do not have any generator matrix
-----------------------------------

In that case, we need to override several methods, namely ``encode``,
``unencode_nocheck`` and probably ``message_space`` (in the case where
the message space is not a vector space). Note that the default
implementation of :meth:`sage.coding.encoder.Encoder.unencode` relies on
``unencode_nocheck``, so reimplementing the former is not necessary.

In our example, it is easy to create an encoder which does not need
a generator matrix to perform the encoding and the unencoding.
We propose the following implementation::

    sage: def encode(self, message):
    ....:     return vector(GF(2), [message] * self.code().length())

    sage: def unencode_nocheck(self, word):
    ....:     return word[0]

    sage: def message_space(self):
    ....:     return GF(2)

Our work here is done.

We need to do one extra thing: set this encoder in the dictionary
of known encoders for the associated code class.
To do that, just add the following line at the end of your file::

   BinaryRepetitionCode._registered_encoders["RepetitionGeneratorMatrixEncoder"] = BinaryRepetitionCodeGeneratorMatrixEncoder

.. NOTE::

    In case you are implementing a generic encoder (an encoder which works
    with any family of linear codes), please add the following statement in
    ``AbstractLinearCode``'s constructor instead:
    ``self._registered_encoders["EncName"] = MyGenericEncoder``.
    This will make it immediately available to any code
    class which inherits from `AbstractLinearCode`.

Summary of the implementation for encoders
------------------------------------------

1. Inherit from :class:`sage.coding.encoder.Encoder`.
2. Add this line in the class' constructor::

      super(ClassName, self).__init__(associated_code)

3. Implement representation methods (not mandatory) ``_repr_``
   and ``_latex_``.
4. Implement ``__eq__``
5. ``__ne__``, ``code`` come with the abstract class.
6. If a generator matrix is known, override ``generator_matrix``.
7. Else override ``encode``, ``unencode_nocheck`` and if needed
   ``message_space``.
8. Add the encoder to ``CodeClass._registered_encoders``.


V. Write a new decoder class
============================

Let us continue by writing a decoder. As before, we need to know what is
required to describe a decoder. We need of course the associated code of
the decoder. We also want to know which ``Encoder`` we should use when we
try to recover the original message from a received word containing errors.
We call this encoder ``connected_encoder``.
As different decoding algorithms do not have the same behaviour
(e.g. probabilistic decoding vs deterministic), we would like to give a few
clues about the type of a decoder. So we can store a list of keywords in the
class parameter ``_decoder_type``.
Eventually, we also need to know the input space of the decoder.
As usual, initializing these parameters can be delegated to the topclass,
and our constructor looks like that::

    sage: from sage.coding.decoder import Decoder
    sage: class BinaryRepetitionCodeMajorityVoteDecoder(Decoder):
    ....:     def __init__(self, code):
    ....:         super((BinaryRepetitionCodeMajorityVoteDecoder, self).__init__(code,
    ....:            code.ambient_space(), "RepetitionGeneratorMatrixEncoder"))

Remember to inherit from :class:`sage.coding.decoder.Decoder`!

As ``_decoder_type`` is actually a class parameter, one should set it
in the file itself, outside of any method.
For readability, we suggest to add this statement at the bottom of the file.
We'll get back to this in a moment.

We also want to override representation methods ``_repr_`` and ``_latex_``::

    sage: def _repr_(self):
    ....:     return "Majority vote-based decoder for the %s" % self.code()
    sage: def _latex_(self):
    ....:     return "\textnormal{Majority vote based-decoder for the } %s" % self.code()

And we want to have an equality check too::

    sage: def __eq__(self, other):
    ....:     return isinstance((other, BinaryRepetitionCodeMajorityVoteDecoder)
    ....:           and self.code() == other.code())

As before, default getter methods are provided by the topclass, namely
:meth:`sage.coding.decoder.Decoder.code`,
:meth:`sage.coding.decoder.Decoder.input_space`,
:meth:`sage.coding.decoder.Decoder.decoder_type` and
:meth:`sage.coding.decoder.Decoder.connected_encoder`.

All we have to do know is to implement the methods related to the decoding.

There are two methods, namely
:meth:`sage.coding.decoder.Decoder.decode_to_code`
and :meth:`sage.coding.decoder.Decoder.decode_to_message`.

By the magic of default implementation, these two are linked, as
``decode_to_message`` calls first ``decode_to_code`` and then
``unencode``, while ``decode_to_code`` calls successively
``decode_to_message`` and ``encode``.
So we only need to implement one of these two, and we choose
to override ``decode_to_code``::

    sage: def decode_to_code(self, word):
    ....:     list_word = word.list()
    ....:     count_one = list_word.count(GF(2).one())
    ....:     n = self.code().length()
    ....:     length = len(list_word)
    ....:     F = GF(2)
    ....:     if count_one > length / 2:
    ....:         return vector(F, [F.one()] * n)
    ....:     elif count_one < length / 2:
    ....:         return vector(F, [F.zero()] * n)
    ....:     else:
    ....:         raise DecodingError("impossible to find a majority")

.. NOTE::

    One notices that if default ``decode_to_code`` calls default
    ``decode_to_message`` and default ``decode_to_message`` calls default
    ``decode_to_code``, if none is overridden and one is called,
    it will end up stuck in an infinite loop. We added a trigger guard
    against this, so if none is overridden and one is called,
    an exception will be raised.

Only one method is missing: one to provide to the user the number of
errors our decoder can decode.
This is the method :meth:`sage.coding.decoder.Decoder.decoding_radius`,
which we override::

    sage: def decoding_radius(self):
    ....:     return (self.code().length()-1) // 2

As for some cases, the decoding might not be precisely known, its
implementation is not mandatory in :class:`sage.coding.decoder.Decoder`'s
subclasses.

We need to do one extra thing: set this encoder in the dictionary of
known decoders for the associated code class.
To do that, just add the following line at the end of your file::

   BinaryRepetitionCode._registered_decoders["MajorityVoteDecoder"] = BinaryRepetitionCodeMajorityVoteDecoder

Also put this line to set ``decoder_type``::

   BinaryRepetitionCode._decoder_type = {"hard-decision", "unique"}


.. NOTE::

    In case you are implementing a generic decoder (a decoder which works
    with any family of linear codes), please add the following statement in
    ``AbstractLinearCode``'s constructor instead:
    ``self._registered_decoders["DecName"] = MyGenericDecoder``.
    This will make it immediately available to any code
    class which inherits from `AbstractLinearCode`.

Summary of the implementation for decoders
------------------------------------------

1. Inherit from :class:`sage.coding.decoder.Decoder`.
2. Add this line in the class' constructor::

      super(ClassName, self).__init__(associated_code, input_space, connected_encoder_name, decoder_type)

3. Implement representation methods (not mandatory) ``_repr_`` and
   ``_latex_``.
4. Implement ``__eq__``.
5. ``__ne__``, ``code``, ``connected_encoder``, ``decoder_type`` come with
   the abstract class.
6. Override ``decode_to_code`` or ``decode_to_message`` and
   ``decoding_radius``.
7. Add the encoder to ``CodeClass._registered_decoders``.

VI. Write a new channel class
=============================

Alongside all these new structures directly related to codes, we also propose
a whole new and shiny structure to experiment on codes, and more specifically
on their decoding.

Indeed, we implemented a structure to emulate real-world communication
channels.

I'll propose here a step-by-step implementation of a dummy channel
for example's sake.

We will implement a very naive channel which works only for words over
:math:`\GF{2}` and flips as many bits as requested by the user.

As channels are not directly related to code families, but more to
vectors and words, we have a specific file, ``channel.py``
to store them.

So we will just add our new class in this file.

For starters, we ask ourselves the eternal question: What do we need to
describe a channel?
Well, we mandatorily need its ``input_space`` and its ``output_space``.
Of course, in most of the cases, the user will be able to provide some extra
information on the channel's behaviour.
In our case, it will be the number of bits to flip (aka the number of errors).

As you might have guess, there is an abstract class to take care
of the mandatory arguments!
Plus, in our case, as this channel only works for vectors
over :math:`\GF{2}`, the input and output spaces are the same.
Let us write the constructor of our new channel class::

    sage: from sage.coding.channel import Channel
    sage: class BinaryStaticErrorRateChannel(Channel):
    ....:     def __init__(self, space, number_errors):
    ....:         if space.base_ring() is not GF(2):
    ....:             raise ValueError("Provided space must be a vector space over GF(2)")
    ....:         if number_errors > space.dimension():
    ....:             raise ValueErrors("number_errors cannot be bigger than input space's dimension")
    ....:         super(BinaryStaticErrorRateChannel, self).__init__(space, space)
    ....:         self._number_errors = number_errors

Remember to inherit from :class:`sage.coding.channel.Channel`!

We also want to override representation methods ``_repr_`` and ``_latex_``::

    sage: def _repr_(self):
    ....:     return ("Binary static error rate channel creating %s errors, of input and output space %s"
    ....:             % (format_interval(no_err), self.input_space()))

    sage: def _latex_(self):
    ....:     return ("\\textnormal{Static error rate channel creating %s errors, of input and output space %s}"
    ....:             % (format_interval(no_err), self.input_space()))

We don't really see any use case for equality methods
(``__eq__`` and ``__ne__``) so do not provide any default implementation.
If one needs these, one can of course override Python's default methods.

We of course want getter methods.
There is a provided default implementation for ``input_space`` and
``output_space``, so we only need one for ``number_errors``::

    sage: def number_errors(self):
    ....:     return self._number_errors

So, now we want a method to actually add errors to words.
As it is the same thing as transmitting messages over a real-world channel,
we propose two methods, ``transmit`` and ``transmit_unsafe``.
As you can guess, ``transmit_unsafe`` tries to transmit the message
without checking if it is in the input space or not, while ``transmit`` checks
this before the transmission... Which means that ``transmit`` has a default
implementation which calls ``transmit_unsafe``.
So we only need to override ``transmit_unsafe``! Let us do it::

    sage: def transmit_unsafe(self, message):
    ....:     w = copy(message)
    ....:     number_err = self.number_errors()
    ....:     V = self.input_space()
    ....:     F = GF(2)
    ....:     for i in sample(range(V.dimension()), number_err):
    ....:         w[i] += F.one()
    ....:     return w

That is it, we now have our new channel class ready to use!

Summary of the implementation for channels
------------------------------------------

1. Inherit from :class:`sage.coding.channel.Channel`.
2. Add this line in the class' constructor:

   .. CODE-BLOCK:: python

      super(ClassName, self).__init__(input_space, output_space)

3. Implement representation methods (not mandatory) ``_repr_`` and
   ``_latex_``.
4. ``input_space`` and ``output_space`` getter methods come with the
   abstract class.
5. Override ``transmit_unsafe``.


VII. Sort our new elements
==========================

As there is many code families and channels in the coding theory library,
we do not wish to store all our classes directly in Sage's global namespace.

We propose several catalog files to store our constructions, namely:

- ``codes_catalog.py``,
- ``encoders_catalog.py``,
- ``decoders_catalog.py`` and
- ``channels_catalog.py``.

Every time one creates a new object, it should be added in the dedicated
catalog file instead of coding theory folder's ``all.py``.

Here it means the following:

- add the following in ``codes_catalog.py``:

   .. CODE-BLOCK:: python

    from sage.coding.repetition_code import BinaryRepetitionCode

- add the following in ``encoders_catalog.py``:

   .. CODE-BLOCK:: python

    from sage.coding.repetition_code import BinaryRepetitionCodeGeneratorMatrixEncoder

- add the following in ``decoders_catalog.py``:

   .. CODE-BLOCK:: python

    from sage.coding.repetition_code import BinaryRepetitionCodeMajorityVoteDecoder

- add the following in ``channels_catalog.py``:

   .. CODE-BLOCK:: python

    from sage.coding.channel import BinaryStaticErrorRateChannel

VIII. Complete code of this tutorial
====================================

If you need some base code to start from, feel free to copy-paste and
derive from the one that follows.

``repetition_code.py`` (with two encoders):

.. CODE-BLOCK:: python

    from sage.coding.linear_code import AbstractLinearCode
    from sage.coding.encoder import Encoder
    from sage.coding.decoder import Decoder
    from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF

    class BinaryRepetitionCode(AbstractLinearCode):

        _registered_encoders = {}
        _registered_decoders = {}

        def __init__(self, length):
            super(BinaryRepetitionCode, self).__init__(GF(2), length, "RepetitionGeneratorMatrixEncoder", "MajorityVoteDecoder")
            self._dimension = 1

        def _repr_(self):
            return "Binary repetition code of length %s" % self.length()

        def _latex_(self):
            return "\textnormal{Binary repetition code of length } %s" % self.length()

        def __eq__(self, other):
            return (isinstance(other, BinaryRepetitionCode)
               and self.length() == other.length()
               and self.dimension() == other.dimension())


    class BinaryRepetitionCodeGeneratorMatrixEncoder(Encoder):

        def __init__(self, code):
            super(BinaryRepetitionCodeGeneratorMatrixEncoder, self).__init__(code)

        def _repr_(self):
            return "Binary repetition encoder for the %s" % self.code()

        def _latex_(self):
            return "\textnormal{Binary repetition encoder for the } %s" % self.code()

        def __eq__(self, other):
            return (isinstance(other, BinaryRepetitionCodeGeneratorMatrixEncoder)
               and self.code() == other.code())

        def generator_matrix(self):
            n = self.code().length()
            return Matrix(GF(2), 1, n, [GF(2).one()] * n)


    class BinaryRepetitionCodeStraightforwardEncoder(Encoder):

        def __init__(self, code):
            super(BinaryRepetitionCodeStraightforwardEncoder, self).__init__(code)

        def _repr_(self):
            return "Binary repetition encoder for the %s" % self.code()

        def _latex_(self):
            return "\textnormal{Binary repetition encoder for the } %s" % self.code()

        def __eq__(self, other):
            return (isinstance(other, BinaryRepetitionCodeStraightforwardEncoder)
               and self.code() == other.code())

        def encode(self, message):
            return vector(GF(2), [message] * self.code().length())

        def unencode_nocheck(self, word):
            return word[0]

        def message_space(self):
            return GF(2)


    class BinaryRepetitionCodeMajorityVoteDecoder(Decoder):

        def __init__(self, code):
            super(BinaryRepetitionCodeMajorityVoteDecoder, self).__init__(code, code.ambient_space(),
               "RepetitionGeneratorMatrixEncoder")

        def _repr_(self):
            return "Majority vote-based decoder for the %s" % self.code()

        def _latex_(self):
            return "\textnormal{Majority vote based-decoder for the } %s" % self.code()


        def __eq__(self, other):
            return (isinstance(other, BinaryRepetitionCodeMajorityVoteDecoder)
               and self.code() == other.code())

        def decode_to_code(self, word):
            list_word = word.list()
            count_one = list_word.count(GF(2).one())
            n = self.code().length()
            length = len(list_word)
            F = GF(2)
            if count_one > length / 2:
                return vector(F, [F.one()] * n)
            elif count_one < length / 2:
               return vector(F, [F.zero()] * n)
            else:
               raise DecodingError("impossible to find a majority")

        def decoding_radius(self):
            return (self.code().length()-1) // 2


    BinaryRepetitionCode._registered_encoders["RepetitionGeneratorMatrixEncoder"] = BinaryRepetitionCodeGeneratorMatrixEncoder
    BinaryRepetitionCode._registered_encoders["RepetitionStraightforwardEncoder"] = BinaryRepetitionCodeStraightforwardEncoder
    BinaryRepetitionCode._registered_decoders["MajorityVoteDecoder"] = BinaryRepetitionCodeMajorityVoteDecoder
    BinaryRepetitionCodeMajorityVoteDecoder._decoder_type = {"hard-decision", "unique"}

``channel.py`` (continued):

.. CODE-BLOCK:: python

    class BinaryStaticErrorRateChannel(Channel):

        def __init__(self, space, number_errors):
            if space.base_ring() is not GF(2):
                raise ValueError("Provided space must be a vector space over GF(2)")
            if number_errors > space.dimension():
                raise ValueErrors("number_errors cannot be bigger than input space's dimension")
            super(BinaryStaticErrorRateChannel, self).__init__(space, space)
            self._number_errors = number_errors

        def _repr_(self):
          return ("Binary static error rate channel creating %s errors, of input and output space %s"
                  % (format_interval(no_err), self.input_space()))

        def _latex_(self):
          return ("\\textnormal{Static error rate channel creating %s errors, of input and output space %s}"
                  % (format_interval(no_err), self.input_space()))

        def number_errors(self):
          return self._number_errors

        def transmit_unsafe(self, message):
            w = copy(message)
            number_err = self.number_errors()
            V = self.input_space()
            F = GF(2)
            for i in sample(range(V.dimension()), number_err):
                w[i] += F.one()
            return w

``codes_catalog.py`` (continued):

.. CODE-BLOCK:: python

    from sage.coding.repetition_code import BinaryRepetitionCode

``encoders_catalog.py`` (continued):

.. CODE-BLOCK:: python

    from sage.coding.repetition_code import (BinaryRepetitionCodeGeneratorMatrixEncoder, BinaryRepetitionCodeStraightforwardEncoder)

``decoders_catalog.py`` (continued):

.. CODE-BLOCK:: python

    from sage.coding.repetition_code import BinaryRepetitionCodeMajorityVoteDecoder

``channels_catalog.py`` (continued):

.. CODE-BLOCK:: python

    from sage.coding.channel import (ErrorErasureChannel, StaticErrorRateChannel, BinaryStaticErrorRateChannel)
