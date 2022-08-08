.. -*- coding: utf-8 -*-
.. linkall

.. _coding_theory:

======================
Coding Theory in Sage
======================

.. MODULEAUTHOR:: David Joyner and Robert Miller (2008), edited by Ralf Stephan
                  for the initial version. David Lucas (2016) for this version.
                  Marketa Slukova (2019) for the latest version.


This tutorial, designed for beginners who want to discover how to use Sage
for their work (research, experimentation, teaching) on coding theory,
will present several key features of Sage's coding theory library
and explain how to find classes and methods you look for.

During this tutorial, we will cover the following parts:

- what can you do with **generic linear codes and associated methods**,
- what can you do with **structured code families**,
- what can you do to   **encode and recover messages, correct errors** and
- what can you do to   **easily add errors to codewords**,
- what can you do with **codes in general**.

The goal of this tutorial is to give a quick overview of what can be done
with the library and how to use the main functionalities.
It is neither a comprehensive description of all methods nor of specific
classes. If one needs some specific information on the behaviour of a
class/method related to coding theory, one should check the documentation
for this class/method.

.. contents:: Table of contents
   :depth: 3
   :class: this-will-duplicate-information-and-it-is-still-useful-here

I. Generic Linear codes and associated methods
==============================================

Let us start with the most basic code one can build: a generic linear code
without any specific structure.

To build such a code, one just need to provide
its generator matrix, as follows::

    sage: G = matrix(GF(3), [[1, 0, 0, 0, 1, 2, 1],
    ....:                    [0, 1, 0, 0, 2, 1, 0],
    ....:                    [0, 0, 1, 2, 2, 2, 2]])
    sage: C = LinearCode(G)

With these lines, you just created a linear code, congratulations!
Note that if you pass a matrix which is not full rank, Sage will turn
it into a full-rank matrix before building the code,
as illustrated in the following example::

    sage: G = matrix(GF(3), [[1, 0, 0, 0, 1, 2, 1],
    ....:                    [0, 1, 0, 0, 2, 1, 0],
    ....:                    [0, 0, 1, 2, 2, 2, 2],
    ....:                    [1, 0, 1, 2, 0, 1, 0]]) #r3 = r0 + r2
    sage: C = LinearCode(G)
    sage: C.generator_matrix()
    [1 0 0 0 1 2 1]
    [0 1 0 0 2 1 0]
    [0 0 1 2 2 2 2]

We now have a linear code... What can we do with it?
As we can a lot of things, let us start with the basic functionalities.

In the example just above, we already asked for the code's
generator matrix. It is also possible to ask the code for
its basic parameters: its *length* and *dimension* as illustrated thereafter::

    sage: C.length()
    7
    sage: C.dimension()
    3

It is also possible to ask for our code's minimum distance::

    sage: C.minimum_distance()
    3

Of course, as ``C`` is a generic linear code, an exhaustive search algorithm
is run to find the minimum distance, which will be slower and slower as the
code grows.

By just typing the name of our code, we get a sentence which briefly
describes it and gives its parameters::

    sage: C
    [7, 3] linear code over GF(3)

As the aim of this tutorial is not to give a comprehensive view of the methods,
we won't describe any other methods.

If one wants to get all methods that can be run on a linear code, one can:

- either check the manual page of the file :ref:`sage.coding.linear_code`
- or type::

      C.<tab>

  in Sage to get a list of all available methods for ``C``.
  Afterwards, typing::

      C.method?

  will show the manual page for ``method``.

.. NOTE::

    Some generic methods require the installation of the optional
    package Guava for Gap. While some work is done to always propose
    a default implementation which *does not* require an optional package,
    there exist some methods which are not up to date - yet.
    If you're receiving an error message related to Gap, please check the
    documentation of the method to verify if Guava has to be installed.


II. Structured code families and an overview of the encoding and decoding system
================================================================================

II.1 Create specific codes in Sage
----------------------------------

Now that we know how to create generic codes, we want to go deeper
and create specific code families. In Sage, all codes families can be
accessed by typing::

    codes.<tab>

Doing so, you will get the comprehensive list of all code families
Sage can build.

For the rest of this section, we will illustrate specific functionalities
of these code families by manipulating
:class:`sage.coding.grs_code.GeneralizedReedSolomonCode`.

So, for starters, we want to create a Generalized Reed-Solomon (GRS) code.

By clicking on the link provided above, or typing::

    codes.GeneralizedReedSolomonCode?

one can access the documentation page for GRS codes, find a definition
of these and learn what is needed to build one in Sage.

Here we choose to build a [12, 6] GRS code over :math:`\GF{13}`.
To do this, we need up to three elements:

- The **list of evaluation points**,
- the **dimension of the code**, and
- optionally, the **list of column multipliers**.

We build our code as follows::

    sage: F = GF(13)
    sage: length, dimension = 12, 6
    sage: evaluation_pts = F.list()[:length]
    sage: column_mults = F.list()[1:length+1]
    sage: C = codes.GeneralizedReedSolomonCode(evaluation_pts, dimension, column_mults)

Our GRS code is now created. We can ask for its parameters, as we did in the
previous section::

    sage: C.length()
    12
    sage: C.dimension()
    6
    sage: C.base_ring()
    Finite Field of size 13

It is also possible to ask for the evaluation points and
the column multipliers by calling
:meth:`sage.coding.grs_code.GeneralizedReedSolomonCode.evaluation_points` and
:meth:`sage.coding.grs_code.GeneralizedReedSolomonCode.column_multipliers`.

Now, if you know some theory for GRS codes, you know that it's especially easy
to compute their minimum distance, which is:
:math:`d = n - k + 1`, where :math:`n` is the length of the code and
:math:`k` is the dimension of the code.

Because Sage knows ``C`` is a GRS code, it will not run the exhaustive
search algorithm presented in section I to find ``C``'s minimum distance
but use the operation introduced above. And you instantly get::

    sage: C.minimum_distance()
    7

All these parameters are summarized inside the string representation
of our code::

    sage: C
    [12, 6, 7] Generalized Reed-Solomon Code over GF(13)

.. NOTE::

    Writing proper classes for code families is a work in progress.
    Some constructions under ``codes.<tab>`` might thus be functions which
    build a generic linear code, and in that case are only able to use
    generic algorithms.
    Please refer to the documentation of a construction to check if it
    is a function or a class.


II.2 Encode and decode in Sage
------------------------------

In the previous part, we learnt how to find specific code families
in Sage and create instances of these families.

In this part, we will learn how to encode and decode.

First of all, we want to generate a codeword to play with.
There is two different ways to do that:

- It is possible to just generate a random element of our code, as follows::

    sage: c = C.random_element()
    sage: c in C
    True

- Alternatively, we can create a message and then encode it into a codeword::

    sage: msg = random_vector(C.base_field(), C.dimension())
    sage: c = C.encode(msg)
    sage: c in C
    True

Either way, we obtained a codeword.
So, we might want to put some errors in it, and try to correct these
errors afterwards. We can obviously do it by changing the values at
some random positions of our codeword, but we propose here something
more general: communication channels.
:class:`sage.coding.channel.Channel` objects are meant
as abstractions for communication channels and for manipulation of
data representation. In this case, we want to emulate a communication channel
which adds some, but not too many, errors to a transmitted word::

    sage: err = 3
    sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), err)
    sage: Chan
    Static error rate channel creating 3 errors, of input and output space Vector space of dimension 12 over Finite Field of size 13
    sage: r = Chan.transmit(c)
    sage: len((c-r).nonzero_positions())
    3

If you want to learn more on Channels, please refer to section IV
of this tutorial.

Thanks to our channel, we got a "received word`, ``r``, as a codeword
with errors on it. We can try to correct the errors and recover the
original codeword::

    sage: c_dec = C.decode_to_code(r)
    sage: c_dec == c
    True

Perhaps we want the original *message* back rather than the codeword.
All we have to do then is to *unencode* it back to the message space::

    sage: m_unenc = C.unencode(c_dec)
    sage: m_unenc == msg
    True

It is also possible to perform the two previous operations
(correct the errors and recover the original message) in one line,
as illustrated below::

    sage: m_unenc2 = C.decode_to_message(r)
    sage: m_unenc2 == msg
    True


III. A deeper view of the Encoder and Decoder structure
=======================================================

In the previous section, we saw that encoding, decoding and unencoding
a vector can be easily done using methods directly on the code object.
These methods are actually shortcuts, added for usability, for when one
does not care more specifically about how encoding and decoding takes place.
At some point, however, one might need more control.

This section will thus go into details on the mechanism
of Encoders and Decoders.

At the core, the three mentioned operations are handled by
:class:`sage.coding.encoder.Encoder` and
:class:`sage.coding.decoder.Decoder`.
These objects possess their own methods to
operate on words. When one calls (as seen above)::

    C.encode(msg)

one actually calls the method :meth:`sage.coding.encoder.Encoder.encode`
on the default encoder of ``C``.
Every code object possess a list of encoders and decoders it can use.
Let us see how one can explore this::

    sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12, GF(59).list()[1:41])
    sage: C.encoders_available()
    ['EvaluationPolynomial', 'EvaluationVector', 'Systematic']
    sage: C.decoders_available()
    ['BerlekampWelch',
     'ErrorErasure',
     'Gao',
     'GuruswamiSudan',
     'InformationSet',
     'KeyEquationSyndrome',
     'NearestNeighbor',
     'Syndrome']

We got a list of the available encoders and decoders for our GRS code.
Rather than using the default ones as we did before,
we can now ask for specific encoder and decoder::

    sage: Evect = C.encoder("EvaluationVector")
    sage: Evect
    Evaluation vector-style encoder for [40, 12, 29] Generalized Reed-Solomon Code over GF(59)
    sage: type(Evect)
    <class 'sage.coding.grs_code.GRSEvaluationVectorEncoder'>
    sage: msg = random_vector(GF(59), C.dimension()) #random
    sage: c = Evect.encode(msg)
    sage: NN = C.decoder("NearestNeighbor")
    sage: NN
    Nearest neighbor decoder for [40, 12, 29] Generalized Reed-Solomon Code over GF(59)

Calling::

    C.encoder(encoder_name)

is actually a short-hand for constructing the encoder manually,
by calling the constructor for
:class:`sage.coding.grs_code.EncoderGRSEvaluationVector` yourself.
If you don't supply ``encoder_name`` to
:meth:`sage.coding.linear_code.AbstractLinearCode.encoder`
you get the default encoder for the code.
:meth:`sage.coding.linear_code.AbstractLinearCode.encoder`
also has an important side-effect: **it caches the constructed encoder**
before returning it. This means that each time one will access the same
``EvaluationVector`` encoder for ``C``, which saves construction time.

All the above things are similar for Decoders.
This reinforces that Encoders and Decoders are rarely constructed but used
many times, which allows them to perform expensive precomputation
at construction or first use, for the benefit of future use.

This gives a good idea of how the elements work internally.
Let us now go a bit more into details on specific points.

III.1 Message spaces
---------------------

The point of an Encoder is to encode messages into the code.
These messages are often just vectors over the base field of the code
and whose length matches the code's dimension.
But it could be anything: vectors over other fields, polynomials, or even
something quite different.
Therefore, each Encoder has a :meth:`sage.coding.encoder.Encoder.message_space`.
For instance, we saw earlier that our GRS code has two possible encoders;
let us investigate the one we left behind in the part just before::

    sage: Epoly = C.encoder("EvaluationPolynomial")
    sage: Epoly
    Evaluation polynomial-style encoder for [40, 12, 29] Generalized Reed-Solomon Code over GF(59)
    sage: Epoly.message_space()
    Univariate Polynomial Ring in x over Finite Field of size 59
    sage: msg_p = Epoly.message_space().random_element(degree=C.dimension()-1); msg_p #random
    31*x^11 + 49*x^10 + 56*x^9 + 31*x^8 + 36*x^6 + 58*x^5 + 9*x^4 + 17*x^3 + 29*x^2 + 50*x + 46

``Epoly`` reflects that GRS codes are often constructed
as evaluations of polynomials, and that a natural way to consider
their messages is as polynomials of degree at most :math:`k-1`,
where :math:`k` is the dimension of the code.
Notice that the message space of ``Epoly`` is all univariate polynomials:
``message_space`` is the ambient space of the messages, and sometimes an Encoder
demands that the messages are actually picked from a subspace hereof.

The default encoder of a code always has a vector space as message space,
so when we call
:meth:`sage.coding.linear_code.AbstractLinearCode.decode_to_message` or
:meth:`sage.coding.linear_code.AbstractLinearCode.unencode` on the code itself,
as illustrated on the first example, this will always return
vectors whose length is the dimension of the code.


III.2 Generator matrices
------------------------

Whenever the message space of an Encoder is a vector space
and it encodes using a linear map, the Encoder will
possess a generator matrix (note that this notion does not
make sense for other types of encoders), which specifies that linear map.

Generator matrices have been placed on Encoder objects since a code
has many generator matrices, and each of these will encode messages differently.
One will also find
:meth:`sage.coding.linear_code.AbstractLinearCode.generator_matrix`
on code objects, but this is again simply a convenience method which forwards
the query to the default encoder.

Let us see this in Sage, using the first encoder we constructed::

    sage: Evect.message_space()
    Vector space of dimension 12 over Finite Field of size 59
    sage: G = Evect.generator_matrix()
    sage: G == C.generator_matrix()
    True

III.3 Decoders and messages
---------------------------

As we saw before, any code has two generic methods for decoding, called
``decode_to_codeword`` and ``decode_to_message``.
Every Decoder also has these two methods, and the methods on the code
simply forward the calls to the default decoder of this code.

There are two reasons for having these two methods: convenience and speed.
Convenience is clear: having both methods provides a useful shortcut
depending on the user's needs.
Concerning speed, some decoders naturally decode directly to a codeword,
while others directly to a message space.
Supporting both methods therefore avoids unnecessary work
in encoding and unencoding.

However, ``decode_to_message`` implies that there is a message space
and an encoding from that space to the code behind the scenes.
A Decoder has methods ``message_space`` and ``connected_encoder``
to inform the user about this. Let us illustrate that by a long example::

    sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[1:41], 3, GF(59).list()[1:41])
    sage: c = C.random_element()
    sage: c in C
    True

    #Create two decoders: Syndrome and Gao
    sage: Syn = C.decoder("KeyEquationSyndrome")
    sage: Gao = C.decoder("Gao")

    #Check their message spaces
    sage: Syn.message_space()
    Vector space of dimension 3 over Finite Field of size 59
    sage: Gao.message_space()
    Univariate Polynomial Ring in x over Finite Field of size 59

    #and now we unencode
    sage: Syn.decode_to_message(c) #random
    (55,9,43)

    sage: Gao.decode_to_message(c) #random
    43*x^2 + 9*x + 55

IV. A deeper look at channels
=============================

In Section II, we briefly introduced the Channel objects as
a way to put errors in a word.
In this section, we will look deeper at their functionality and
introduce a second Channel.

    .. NOTE::

        Once again, we chose a specific class as a running example through
        all this section, as we do not want to make an exhaustive catalog
        of all channels.
        If one wants to get this list, one can access it by typing::

            channels.<tab>

        in Sage.

Consider again the
:meth:`sage.coding.channel.ChannelStaticErrorRate` from before.
This is a channel that places errors in the transmitted vector
but within controlled boundaries.
We can describe these boundaries in two ways:

- The first one was illustrated in Section II and consists in passing
  an integer, as shown below::

    sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
    sage: t = 14
    sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), t)
    sage: Chan
    Static error rate channel creating 14 errors, of input and output space Vector space of dimension 40 over Finite Field of size 59

- We can also pass a tuple of two integers, the first smaller than the second.
  Then each time a word is transmitted, a random number of errors
  between these two integers will be added::

    sage: t = (1, 14)
    sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), t)
    sage: Chan
    Static error rate channel creating between 1 and 14 errors, of input and output space Vector space of dimension 40 over Finite Field of size 59

We already know that a channel has a
:meth:`sage.coding.channel.Channel.transmit` method which will
perform transmission over the channel; in this case it will return
the transmitted word with some errors in it.
This method will always check if the provided word belongs to
the input space of the channel.
In a case one is absolutely certain that one's word is in the input space,
one might want to avoid this check, which is time consuming - especially
if one is simulating millions of transmissions.
For this usage there is
:meth:`sage.coding.channel.Channel.transmit_unsafe` which does
the same as
:meth:`sage.coding.channel.Channel.transmit`
but without checking the input, as illustrated thereafter::

    sage: c = C.random_element()
    sage: c in C
    True
    sage: c_trans = Chan.transmit_unsafe(c)
    sage: c_trans in C
    False

Note there exists a useful shortcut for
:meth:`sage.coding.channel.Channel.transmit` ::

    sage: r = Chan(c)
    sage: r in C
    False

A channel for errors and erasures
---------------------------------

Let us introduce a new Channel object which adds errors and erasures.
When it transmits a word, it both adds some errors
as well as it erases some positions::

    sage: Chan = channels.ErrorErasureChannel(C.ambient_space(), 3, 4)
    sage: Chan
    Error-and-erasure channel creating 3 errors and 4 erasures of input space Vector space of dimension 40 over Finite Field of size 59 and output space The Cartesian product of (Vector space of dimension 40 over Finite Field of size 59, Vector space of dimension 40 over Finite Field of size 2)


The first parameter is the input space of the channel.
The next two are (respectively) the number of errors
and the number or erasures.
Each of these can be tuples too, just as it was with
:class:`sage.coding.channel.StaticErrorRateChannel`.
As opposed to this channel though, the output of
:class:`sage.coding.channel.ErrorErasureChannel`
is not the same as its input space, i.e. the ambient space of C.
Rather, it will return two vectors: the first is the transmitted word
with the errors added and erased positions set to 0.
The second one is the erasure vector whose erased positions contain ones.
This is reflected in :meth:`sage.coding.channel.output_space`::

    sage: C = codes.random_linear_code(GF(7), 10, 5)
    sage: Chan.output_space()
    The Cartesian product of (Vector space of dimension 40 over Finite Field of size 59, Vector space of dimension 40 over Finite Field of size 2)
    sage: Chan(c) # random
    ((0, 3, 6, 4, 4, 0, 1, 0, 0, 1),
     (1, 0, 0, 0, 0, 1, 0, 0, 1, 0))

Note it is guaranteed by construction that errors and erasures
will never overlap, so when you ask for ``e`` errors and ``t`` erasures,
you will always receive a vector with ``e`` errors and ``t`` erased positions.


V. Codes in General
===================

So far we have talked only about codes which are linear and over the Hamming
metric. Sage also supports codes which are non-linear and/or over a different
metric. Since working with these usually involves creating new classes, they are
covered in detail in :ref:`structures_in_coding_theory`.

For this tutorial, we include a small example of a linear code over the rank
metric. Unlike the Hamming metric, where the distance between two words is the
number of positions in which their differ, rank metric takes distance to be the
rank of the difference of two codewords, which are represented as matrices.

We are going to create a generic linear code over the rank metric without any
additional structure. To build such a code, we only need a generator matrix::

    sage: G = Matrix(GF(4), [[1,1,0], [0,0,1]])
    sage: C = codes.LinearRankMetricCode(G)

We can do all the things that we did with our example of a linear code over the
Hamming metric. Therefore we are going to focus on the different metric.

Take a word in our code::

    sage: c = C[1]
    sage: c
    (1, 1, 0)

Over the usual Hamming metric, the weight of this word would be `2`. However,
over the rank metric, we get a different result::

    sage: C.rank_weight_of_vector(c)
    1

The weight of a word in the rank metric is simply the rank of the matrix form of
the word::

    sage: C.matrix_form_of_vector(c)
    [1 1 0]
    [0 0 0]

As we said before, the distance between two words is the rank of their
difference::

    sage: d = C[2]
    sage: d
    (z2, z2, 0)
    sage: C.rank_distance_between_vectors(c, d)
    1

Even though the words `c` and `d` differ in two positions, their distance over
the rank metric is `1`::

    sage: C.matrix_form_of_vector(c - d)
    [1 1 0]
    [1 1 0]

For more details on the linear rank metric code, see
:ref:`sage.coding.linear_rank_metric`.


VI. Conclusion - Afterword
==========================

This last section concludes our tutorial on coding theory.

After reading this, you should know enough to create
and manipulate codes in Sage!

We did not illustrate all the content of the library in this tutorial.
For instance, we did not mention how Sage manages bounds on codes.

All objects, constructions and methods related to coding theory are hidden
under the prefix ``codes`` in Sage.

For instance, it is possible to find all encoders you can build by typing::

    codes.encoders.<tab>

So, if you are looking for a specific object related to code, you should always
type::

    codes.<tab>

and check if there's a subcategory which matches your needs.

Despite all the hard work we put on it, there's always much to do!

Maybe at some point you might want to create your own codes for Sage.
If it's the case and if you don't know how to do that, don't panic!
We also wrote a tutorial for this specific case, which you can find here:
:ref:`structures_in_coding_theory`.
