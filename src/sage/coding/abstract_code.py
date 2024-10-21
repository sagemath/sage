# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.modules sage.rings.finite_rings
r"""
Codes

Class supporting methods available for any type of code (linear, non-linear) and
over any metric (Hamming, rank).

There are further abstract classes representing certain types of codes. For
linear codes,
:class:`~sage.coding.linear_code_no_metric.AbstractLinearCodeNoMetric` contains
all the methods that any linear code can use regardless of its metric.
Inheriting from this class are base classes for linear codes over specific
metrics. For example, :class:`~sage.coding.linear_code.AbstractLinearCode` is a
base class for all linear codes over the Hamming metric.

Take the class :class:`~sage.coding.hamming_code.HammingCode`. This
class inherits from :class:`~sage.coding.linear_code.AbstractLinearCode`, since
it is a linear code over the Hamming metric.
:class:`~sage.coding.linear_code.AbstractLinearCode` then inherits from
:class:`~sage.coding.linear_code_no_metric.AbstractLinearCodeNoMetric`, since it
is a linear code. Finally, this class inherits from
:class:`~sage.coding.abstract_code.AbstractCode`, since it is a code.


The following diagram shows the inheritance relationship in the coding module::

    AbstractCode
    + AbstractLinearCodeNoMetric
    | + AbstractLinearCode
    | | + ParityCheckCode
    | | + HammingCode
    | | + CyclicCode
    | | + BCHCode
    | | + GolayCode
    | | + ReedMullerCode
    | | + GeneralizedReedSolomonCode
    | | + GoppaCode
    | + AbstractLinearRankMetricCode

Any class inheriting from AbstractCode can use the encode/decode framework.

The encoder/decoder framework within the coding module offers the creation and
use of encoders/decoders independently of codes. An encoder encodes a message
into a codeword. A decoder decodes a word into a codeword or a message,
possibly with error-correction.

Instead of creating specific encoders/decoders for every code family, some
encoders/decoders can be used by multiple code families. The encoder/decoder
framework enables just that. For example,
:class:`~sage.coding.linear_code.LinearCodeGeneratorMatrixEncoder`
can be used by any code that has a generator matrix.  Similarly,
:class:`~sage.coding.linear_code.LinearCodeNearestNeighborDecoder` can be used
for any linear code with Hamming metric.

When creating a new code family, investigate the encoder/decoder catalogs,
``codes.encoders`` and ``codes.decoders``, to see if there are suitable
encoders/decoders for your code family already implemented. If this is the case,
follow the instructions in :class:`AbstractCode` to set these up.

A new encoder must have the following methods:

- ``encode`` -- method encoding a message into a codeword
- ``unencode`` -- method decoding a codeword into a message
- ``message_space`` -- ambient space of messages that can be encoded
- ``code`` -- code of the encoder

For more information about the Encoder class, see
:class:`~sage.coding.encoder.Encoder`

A new decoder must have the following methods:

- ``decode_to_code`` or ``decode_to_message`` -- method decoding a word from the
  input space into either a codeword or a message
- ``input_space`` -- ambient space of words that can be decoded
- ``code`` -- code of the decoder

For more information about the Decoder class, see
:class:`~sage.coding.decoder.Decoder`
"""

from sage.structure.parent import Parent
from sage.misc.cachefunc import cached_method
from copy import copy
from sage.rings.integer import Integer

import inspect
from sage.misc.sageinspect import sage_getargspec


def _explain_constructor(cl):
    r"""
    Internal function for use error messages when constructing encoders and decoders.

    EXAMPLES::

        sage: from sage.coding.linear_code import LinearCodeSyndromeDecoder
        sage: from sage.coding.abstract_code import _explain_constructor
        sage: cl = LinearCodeSyndromeDecoder
        sage: _explain_constructor(cl)
        "The constructor requires no arguments.\nIt takes the optional
        arguments ['maximum_error_weight'].\nSee the documentation of
        sage.coding.linear_code.LinearCodeSyndromeDecoder for more details."

        sage: from sage.coding.information_set_decoder import LinearCodeInformationSetDecoder
        sage: cl = LinearCodeInformationSetDecoder
        sage: _explain_constructor(cl)
        "The constructor requires the arguments ['number_errors'].\nIt takes the optional arguments ['algorithm'].\nIt accepts unspecified arguments as well.\nSee the documentation of sage.coding.information_set_decoder.LinearCodeInformationSetDecoder for more details."
    """
    if inspect.isclass(cl):
        argspec = sage_getargspec(cl.__init__)
        skip = 2 # skip the self and code arguments
    else:
        # Not a class, assume it's a factory function posing as a class
        argspec = sage_getargspec(cl)
        skip = 1 # skip code argument
    if argspec.defaults:
        args = argspec.args[skip:-len(argspec.defaults)]
        kwargs = argspec.args[-len(argspec.defaults):]
        opts = "It takes the optional arguments {}.".format(kwargs)
    else:
        args = argspec.args[skip:]
        opts = "It takes no optional arguments."
    if args:
        reqs = "The constructor requires the arguments {}.".format(args)
    else:
        reqs = "The constructor requires no arguments."
    if argspec.varargs or argspec.varkw:
        var = "It accepts unspecified arguments as well.\n"
    else:
        var = ""
    return ("{}\n{}\n{}See the documentation of {}.{} for more details."
            .format(reqs, opts, var, cl.__module__, cl.__name__))


class AbstractCode(Parent):
    r"""
    Abstract class for codes.

    This class contains all the methods that can be used on any code
    and on any code family. As opposed to
    :class:`sage.coding.linear_code.AbstractLinearCode`, this class makes no
    assumptions about linearity, metric, finiteness or the number of alphabets.

    The abstract notion of "code" that is implicitly used for this class is any
    enumerable subset of a cartesian product `A_1 \times A_2 \times \ldots
    \times A_n` for some sets `A_i`. Note that this class makes no attempt to
    directly represent the code in this fashion, allowing subclasses to make the
    appropriate choices. The notion of metric is also not mathematically
    enforced in any way, and is simply stored as a string value.

    Every code-related class should inherit from this abstract class.

    To implement a code, you need to:

    - inherit from :class:`AbstractCode`

    - call :class:`AbstractCode` ``__init__`` method in the subclass constructor.
      Example: ``super().__init__(length, "EncoderName",
      "DecoderName", "metric")``. "EncoderName" and "DecoderName" are set to
      ``None`` by default, a generic code class such as AbstractCode does
      not necessarily have to have general encoders/decoders. However, if you
      want to use the encoding/decoding methods, you have to add these.

    - since this class does not specify any category, it is highly recommended
      to set up the category framework in the subclass. To do this, use the
      ``Parent.__init__(self, base, facade, category)`` function in the subclass
      constructor. A good example is in
      :class:`sage.coding.linear_code.AbstractLinearCode`.

    - it is also recommended to override the ``ambient_space`` method, which is
      required by ``__call__``

    - to use the encoder/decoder framework, one has to set up the category and
      related functions ``__iter__`` and ``__contains__``. A good example is in
      :class:`sage.coding.linear_code.AbstractLinearCode`.

    - add the following two lines on the class level::

          _registered_encoders = {}
          _registered_decoders = {}


    - fill the dictionary of its encoders in ``sage.coding.__init__.py`` file.
      Example: I want to link the encoder ``MyEncoderClass`` to ``MyNewCodeClass``
      under the name ``MyEncoderName``.
      All I need to do is to write this line in the ``__init__.py`` file:
      ``MyNewCodeClass._registered_encoders["NameOfMyEncoder"] = MyEncoderClass``
      and all instances of ``MyNewCodeClass`` will be able to use instances of
      ``MyEncoderClass``.

    - fill the dictionary of its decoders in ``sage.coding.__init__`` file.
      Example: I want to link the encoder ``MyDecoderClass`` to ``MyNewCodeClass``
      under the name ``MyDecoderName``.
      All I need to do is to write this line in the ``__init__.py`` file:
      ``MyNewCodeClass._registered_decoders["NameOfMyDecoder"] = MyDecoderClass``
      and all instances of ``MyNewCodeClass`` will be able to use instances of
      ``MyDecoderClass``.


    As the class :class:`AbstractCode` is not designed to be instantiated, it
    does not have any representation methods. You should implement ``_repr_``
    and ``_latex_`` methods in the subclass.
    """

    def __init__(self, length, default_encoder_name=None,
                 default_decoder_name=None, metric='Hamming'):
        r"""
        Initialize mandatory parameters that any code shares.

        This method only exists for inheritance purposes as it initializes
        parameters that need to be known by every code. The class
        :class:`sage.coding.abstract_code.AbstractCode` should never be
        directly instantiated.

        INPUT:

        - ``length`` -- the length of ``self`` (a Python int or a Sage Integer,
          must be > 0)

        - ``default_encoder_name`` -- (default: ``None``) the name of
          the default encoder of ``self``

        - ``default_decoder_name`` -- (default: ``None``) the name of
          the default decoder of ``self``

        - ``metric`` -- (default: ``Hamming``) the name of the metric of ``self``

        EXAMPLES:

        The following example demonstrates how to use a subclass of ``AbstractCode``
        for representing a new family of codes::

            sage: from sage.coding.abstract_code import AbstractCode
            sage: class MyCodeFamily(AbstractCode):
            ....:   def __init__(self, length):
            ....:       super().__init__(length)
            ....:   def __iter__(self):
            ....:       for i in range(self.length() + 1):
            ....:            yield vector([1 for j in range(i)] + [0 for k in range(i, self.length())])
            ....:   def __contains__(self, word):
            ....:       return word in list(self)
            ....:   def _repr_(self):
            ....:       return "Dummy code of length {}".format(self.length())

        We now instantiate a member of our newly made code family::

            sage: C = MyCodeFamily(6)

        We can check its existence and parameters::

            sage: C
            Dummy code of length 6

        We can list its elements and check if an element is in the code::

            sage: list(C)
            [(0, 0, 0, 0, 0, 0),
            (1, 0, 0, 0, 0, 0),
            (1, 1, 0, 0, 0, 0),
            (1, 1, 1, 0, 0, 0),
            (1, 1, 1, 1, 0, 0),
            (1, 1, 1, 1, 1, 0),
            (1, 1, 1, 1, 1, 1)]
            sage: vector((0, 1, 0, 0, 0, 1)) in C
            False
            sage: vector((1, 1, 1, 0, 0, 0)) in C
            True

        And coming from AbstractCode code::

            sage: C.metric()
            'Hamming'

        TESTS:

        If the length field is neither a Python int nor a Sage Integer, it will
        raise a exception::

            sage: C = MyCodeFamily(10.0)
            Traceback (most recent call last):
            ...
            ValueError: length must be a Python int or a Sage Integer

        If the length of the code is not a nonzero positive integer
        (See :issue:`21326`), it will raise an exception::

            sage: C = MyCodeFamily(0)
            Traceback (most recent call last):
            ...
            ValueError: length must be a nonzero positive integer
        """

        if not isinstance(length, (int, Integer)):
            raise ValueError("length must be a Python int or a Sage Integer")
        if length <= 0:
            raise ValueError("length must be a nonzero positive integer")

        self._length = length
        self._metric = metric

        self._default_decoder_name = default_decoder_name
        self._default_encoder_name = default_encoder_name
        if not self._default_decoder_name:
            self._registered_encoders = {}
        if not self._default_encoder_name:
            self._registered_decoders = {}

    def __getstate__(self):
        """
        Used for pickling codes.

        TESTS::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: '_registered_encoders' in C.__getstate__()
            True
        """
        d = super().__getstate__()
        d['_registered_encoders'] = self._registered_encoders
        d['_registered_decoders'] = self._registered_decoders
        return d

    def __iter__(self):
        r"""
        Return an error message requiring to override ``__iter__`` in ``self``.

        As one has to implement specific category related methods (`__iter__` and
        `__contains__`) when writing a new code class which inherits from
        :class:`AbstractCode`, the generic call to `__iter__` has to fail.

        EXAMPLES:

        We create a new code class::

            sage: from sage.coding.abstract_code import AbstractCode
            sage: class MyCode(AbstractCode):
            ....:    def __init__(self):
            ....:        super().__init__(10)

        We check we get a sensible error message while asking for an
        iterator over the elements of our new class::

            sage: C = MyCode()
            sage: list(C)
            Traceback (most recent call last):
            ...
            RuntimeError: Please override __iter__ in the implementation of <class '__main__.MyCode'>
        """
        raise RuntimeError("Please override __iter__ in the implementation of {}".format(self.parent()))

    def __contains__(self, c):
        r"""
        Return an error message requiring to override ``__contains__`` in ``self``.

        As one has to implement specific category related methods (`__iter__` and
        `__contains__`) when writing a new code class which inherits from
        :class:`AbstractCode`, the generic call to `__contains__` has to fail.

        EXAMPLES:

        We create a new code class::

            sage: from sage.coding.abstract_code import AbstractCode
            sage: class MyCode(AbstractCode):
            ....:    def __init__(self, length):
            ....:        super().__init__(length)

        We check we get a sensible error message while asking if an element is
        in our new class::

            sage: C = MyCode(3)
            sage: vector((1, 0, 0, 0, 0, 1, 1)) in C
            Traceback (most recent call last):
            ...
            RuntimeError: Please override __contains__ in the implementation of <class '__main__.MyCode'>
        """
        raise RuntimeError("Please override __contains__ in the implementation of {}".format(self.parent()))

    def ambient_space(self):
        r"""
        Return an error stating ``ambient_space`` of ``self`` is not implemented.

        This method is required by :meth:`__call__`.

        EXAMPLES::

            sage: from sage.coding.abstract_code import AbstractCode
            sage: class MyCode(AbstractCode):
            ....:    def __init__(self, length):
            ....:        super().__init__(length)
            sage: C = MyCode(3)
            sage: C.ambient_space()
            Traceback (most recent call last):
            ...
            NotImplementedError: No ambient space implemented for this code.
        """
        raise NotImplementedError("No ambient space implemented for this code.")

    def __call__(self, m):
        r"""
        Return either ``m`` if it is a codeword or ``self.encode(m)``
        if it is an element of the message space of the encoder used by
        ``encode``.

        This implementation depends on :meth:`ambient_space`.

        INPUT:

        - ``m`` -- a vector whose length equals to code's length or an element
          of the message space used by ``encode``

        - ``**kwargs`` -- extra arguments are forwarded to ``encode``

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: word = vector((0, 1, 1, 0))
            sage: C(word)
            (1, 1, 0, 0, 1, 1, 0)

            sage: c = C.random_element()
            sage: C(c) == c
            True

        TESTS:

        If one passes a vector which belongs to the ambient space, it has to be a codeword.
        Otherwise, an exception is raised::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: word = vector((0, 1, 1, 0, 0, 1, 0))
            sage: C(word)
            Traceback (most recent call last):
            ...
            ValueError: If the input is a vector which belongs to the ambient space, it has to be a codeword
        """
        if m in self.ambient_space():
            if m in self:
                return m
            else:
                raise ValueError("If the input is a vector which belongs to the ambient space, it has to be a codeword")
        else:
            return self.encode(m)

    def _repr_(self):
        r"""
        Return an error message requiring to override ``_repr_`` in ``self``.

        As one has to implement specific representation methods (`_repr_` and
        `_latex_`) when writing a new code class which inherits from
        :class:`AbstractCode`, the generic call to `_repr_` has to fail.

        EXAMPLES:

        We create a new code class::

            sage: from sage.coding.abstract_code import AbstractCode
            sage: class MyCode(AbstractCode):
            ....:    def __init__(self):
            ....:        super().__init__(10)

        We check we get a sensible error message while asking for a string
        representation of an instance of our new class::

            sage: C = MyCode()
            sage: C #random
            Traceback (most recent call last):
            ...
            RuntimeError: Please override _repr_ in the implementation of <class '__main__.MyCode'>
        """
        raise RuntimeError("Please override _repr_ in the implementation of {}".format(self.parent()))

    def _latex_(self):
        r"""
        Return an error message requiring to override ``_latex_`` in ``self``.

        As one has to implement specific representation methods (`_repr_` and
        `_latex_`) when writing a new code class which inherits from
        :class:`AbstractCode`, the generic call to `_latex_` has to fail.

        EXAMPLES:

        We create a new code class::

            sage: from sage.coding.abstract_code import AbstractCode
            sage: class MyCode(AbstractCode):
            ....:    def __init__(self):
            ....:        super().__init__(10)

        We check we get a sensible error message while asking for a string
        representation of an instance of our new class::

            sage: C = MyCode()
            sage: latex(C)
            Traceback (most recent call last):
            ...
            RuntimeError: Please override _latex_ in the implementation of <class '__main__.MyCode'>
        """
        raise RuntimeError("Please override _latex_ in the implementation of {}".format(self.parent()))

    def list(self):
        r"""
        Return a list of all elements of this code.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: Clist = C.list()
            sage: Clist[5]; Clist[5] in C
            (1, 0, 1, 0, 1, 0, 1)
            True
        """
        return list(self)

    def length(self):
        r"""
        Return the length of this code.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.length()
            7
        """
        return self._length

    def metric(self):
        """
        Return the metric of ``self``.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.metric()
            'Hamming'
        """
        return self._metric

###################### Encoding-Decoding #######################################

    def add_decoder(self, name, decoder):
        r"""
        Add an decoder to the list of registered decoders of ``self``.

        .. NOTE::

            This method only adds ``decoder`` to ``self``, and not to any member of the class
            of ``self``. To know how to add an :class:`sage.coding.decoder.Decoder`, please refer
            to the documentation of :class:`AbstractCode`.

        INPUT:

        - ``name`` -- the string name for the decoder

        - ``decoder`` -- the class name of the decoder

        EXAMPLES:

        First of all, we create a (very basic) new decoder::

            sage: class MyDecoder(sage.coding.decoder.Decoder):
            ....:   def __init__(self, code):
            ....:       super().__init__(code)
            ....:   def _repr_(self):
            ....:       return "MyDecoder decoder with associated code %s" % self.code()

        We now create a new code::

            sage: C = codes.HammingCode(GF(2), 3)

        We can add our new decoder to the list of available decoders of C::

            sage: C.add_decoder("MyDecoder", MyDecoder)
            sage: sorted(C.decoders_available())
            ['InformationSet', 'MyDecoder', 'NearestNeighbor', 'Syndrome']

        We can verify that any new code will not know MyDecoder::

            sage: C2 = codes.HammingCode(GF(2), 3)
            sage: sorted(C2.decoders_available())
            ['InformationSet', 'NearestNeighbor', 'Syndrome']

        TESTS:

        It is impossible to use a name which is in the dictionary of available decoders::

            sage: C.add_decoder("Syndrome", MyDecoder)
            Traceback (most recent call last):
            ...
            ValueError: There is already a registered decoder with this name
        """
        if self._registered_decoders == self.__class__._registered_decoders:
            self._registered_decoders = copy(self._registered_decoders)
            reg_dec = self._registered_decoders
            if name in reg_dec:
                raise ValueError("There is already a registered decoder with this name")
            reg_dec[name] = decoder
        else:
            if name in self._registered_decoders:
                raise ValueError("There is already a registered decoder with this name")
            reg_dec[name] = decoder

    def add_encoder(self, name, encoder):
        r"""
        Add an encoder to the list of registered encoders of ``self``.

        .. NOTE::

            This method only adds ``encoder`` to ``self``, and not to any member of the class
            of ``self``. To know how to add an :class:`sage.coding.encoder.Encoder`, please refer
            to the documentation of :class:`AbstractCode`.

        INPUT:

        - ``name`` -- the string name for the encoder

        - ``encoder`` -- the class name of the encoder

        EXAMPLES:

        First of all, we create a (very basic) new encoder::

            sage: class MyEncoder(sage.coding.encoder.Encoder):
            ....:   def __init__(self, code):
            ....:       super().__init__(code)
            ....:   def _repr_(self):
            ....:       return "MyEncoder encoder with associated code %s" % self.code()

        We now create a new code::

            sage: C = codes.HammingCode(GF(2), 3)

        We can add our new encoder to the list of available encoders of C::

            sage: C.add_encoder("MyEncoder", MyEncoder)
            sage: sorted(C.encoders_available())
            ['MyEncoder', 'Systematic']

        We can verify that any new code will not know MyEncoder::

            sage: C2 = codes.HammingCode(GF(2), 3)
            sage: sorted(C2.encoders_available())
            ['Systematic']

        TESTS:

        It is impossible to use a name which is in the dictionary of available encoders::

            sage: C.add_encoder("Systematic", MyEncoder)
            Traceback (most recent call last):
            ...
            ValueError: There is already a registered encoder with this name
        """
        if self._registered_encoders == self.__class__._registered_encoders:
            self._registered_encoders = copy(self._registered_encoders)
            reg_enc = self._registered_encoders
            if name in reg_enc:
                raise ValueError("There is already a registered encoder with this name")
            reg_enc[name] = encoder
        else:
            if name in self._registered_encoders:
                raise ValueError("There is already a registered encoder with this name")
            reg_enc[name] = encoder

    def decode_to_code(self, word, decoder_name=None, *args, **kwargs):
        r"""
        Correct the errors in ``word`` and returns a codeword.

        INPUT:

        - ``word`` -- an element in the ambient space as ``self``

        - ``decoder_name`` -- (default: ``None``) name of the decoder which will be used
          to decode ``word``. The default decoder of ``self`` will be used if
          default value is kept.

        - ``args``, ``kwargs`` -- all additional arguments are forwarded to :meth:`decoder`

        OUTPUT: a vector of ``self``

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0], [1,0,0,1,1,0,0],
            ....:                    [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: word = vector(GF(2), (1, 1, 0, 0, 1, 1, 0))
            sage: w_err = word + vector(GF(2), (1, 0, 0, 0, 0, 0, 0))
            sage: C.decode_to_code(w_err)
            (1, 1, 0, 0, 1, 1, 0)

        It is possible to manually choose the decoder amongst the list of the available ones::

            sage: sorted(C.decoders_available())
            ['InformationSet', 'NearestNeighbor', 'Syndrome']
            sage: C.decode_to_code(w_err, 'NearestNeighbor')
            (1, 1, 0, 0, 1, 1, 0)
        """
        D = self.decoder(decoder_name, *args, **kwargs)
        return D.decode_to_code(word)

    def decode_to_message(self, word, decoder_name=None, *args, **kwargs):
        r"""
        Correct the errors in word and decodes it to the message space.

        INPUT:

        - ``word`` -- an element in the ambient space as ``self``

        - ``decoder_name`` -- (default: ``None``) name of the decoder which will be used
          to decode ``word``. The default decoder of ``self`` will be used if
          default value is kept.

        - ``args``, ``kwargs`` -- all additional arguments are forwarded to :meth:`decoder`

        OUTPUT: a vector of the message space of ``self``

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0], [1,0,0,1,1,0,0],
            ....:                    [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: word = vector(GF(2), (1, 1, 0, 0, 1, 1, 0))
            sage: C.decode_to_message(word)
            (0, 1, 1, 0)

        It is possible to manually choose the decoder amongst the list of the available ones::

            sage: sorted(C.decoders_available())
            ['InformationSet', 'NearestNeighbor', 'Syndrome']
            sage: C.decode_to_message(word, 'NearestNeighbor')
            (0, 1, 1, 0)
        """
        return self.unencode(self.decode_to_code(word, decoder_name, *args, **kwargs), **kwargs)

    @cached_method
    def decoder(self, decoder_name=None, *args, **kwargs):
        r"""
        Return a decoder of ``self``.

        INPUT:

        - ``decoder_name`` -- (default: ``None``) name of the decoder which will be
          returned. The default decoder of ``self`` will be used if
          default value is kept.

        - ``args``, ``kwargs`` -- all additional arguments will be forwarded to the constructor of the decoder
          that will be returned by this method

        OUTPUT: a decoder object

        Besides creating the decoder and returning it, this method also stores
        the decoder in a cache. With this behaviour, each decoder will be created
        at most one time for ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0], [1,0,0,1,1,0,0],
            ....:                    [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C.decoder()
            Syndrome decoder for [7, 4] linear code over GF(2) handling errors of weight up to 1

        If there is no decoder for the code, we return an error::

            sage: from sage.coding.abstract_code import AbstractCode
            sage: class MyCodeFamily(AbstractCode):
            ....:   def __init__(self, length, field):
            ....:       sage.coding.abstract_code.AbstractCode.__init__(self, length)
            ....:       Parent.__init__(self, base=field, facade=False, category=Sets())
            ....:       self._field = field
            ....:   def field(self):
            ....:       return self._field
            ....:   def _repr_(self):
            ....:       return "%d dummy code over GF(%s)" % (self.length(), self.field().cardinality())
            sage: D = MyCodeFamily(5, GF(2))
            sage: D.decoder()
            Traceback (most recent call last):
            ...
            NotImplementedError: No decoder implemented for this code.

        If the name of a decoder which is not known by ``self`` is passed,
        an exception will be raised::

            sage: sorted(C.decoders_available())
            ['InformationSet', 'NearestNeighbor', 'Syndrome']
            sage: C.decoder('Try')
            Traceback (most recent call last):
            ...
            ValueError: There is no Decoder named 'Try'.
            The known Decoders are: ['InformationSet', 'NearestNeighbor', 'Syndrome']

        Some decoders take extra arguments. If the user forgets to supply these,
        the error message attempts to be helpful::

            sage: C.decoder('InformationSet')
            Traceback (most recent call last):
            ...
            ValueError: Constructing the InformationSet decoder failed,
            possibly due to missing or incorrect parameters.
            The constructor requires the arguments ['number_errors'].
            It takes the optional arguments ['algorithm'].
            It accepts unspecified arguments as well. See the documentation of
            sage.coding.information_set_decoder.LinearCodeInformationSetDecoder
            for more details.
        """
        if not self._default_decoder_name:
            raise NotImplementedError("No decoder implemented for this code.")
        if decoder_name is None:
            decoder_name = self._default_decoder_name
        if decoder_name in self._registered_decoders:
            decClass = self._registered_decoders[decoder_name]
            try:
                return decClass(self, *args, **kwargs)
            except TypeError:
                raise ValueError(
                        "Constructing the {0} decoder failed, possibly due "
                        "to missing or incorrect parameters.\n{1}".format(
                            decoder_name, _explain_constructor(decClass)))
        else:
            raise ValueError(
                    "There is no Decoder named '{0}'. The known Decoders are: "
                    "{1}".format(decoder_name, self.decoders_available()))

    def decoders_available(self, classes=False):
        r"""
        Return a list of the available decoders' names for ``self``.

        INPUT:

        - ``classes`` -- boolean (default: ``False``); if ``classes`` is set to
          ``True``, return instead a :class:`dict` mapping available decoder
          name to the associated decoder class

        OUTPUT: list of strings, or a :class:`dict` mapping strings to classes

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0], [1,0,0,1,1,0,0],
            ....:                    [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C.decoders_available()
            ['InformationSet', 'NearestNeighbor', 'Syndrome']

            sage: dictionary = C.decoders_available(True)
            sage: sorted(dictionary.keys())
            ['InformationSet', 'NearestNeighbor', 'Syndrome']
            sage: dictionary['NearestNeighbor']
            <class 'sage.coding.linear_code.LinearCodeNearestNeighborDecoder'>
        """
        if classes:
            return copy(self._registered_decoders)

        return sorted(self._registered_decoders)

    def encode(self, word, encoder_name=None, *args, **kwargs):
        r"""
        Transform an element of a message space into a codeword.

        INPUT:

        - ``word`` -- an element of a message space of the code

        - ``encoder_name`` -- (default: ``None``) name of the encoder which
          will be used to encode ``word``. The default encoder of ``self`` will
          be used if default value is kept.

        - ``args``, ``kwargs`` -- all additional arguments are forwarded to the
          construction of the encoder that is used

        One can use the following shortcut to encode a word ::

            C(word)

        OUTPUT: a vector of ``self``

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0], [1,0,0,1,1,0,0],
            ....:                    [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: word = vector((0, 1, 1, 0))
            sage: C.encode(word)
            (1, 1, 0, 0, 1, 1, 0)
            sage: C(word)
            (1, 1, 0, 0, 1, 1, 0)

        It is possible to manually choose the encoder amongst the list of the available ones::

            sage: sorted(C.encoders_available())
            ['GeneratorMatrix', 'Systematic']
            sage: word = vector((0, 1, 1, 0))
            sage: C.encode(word, 'GeneratorMatrix')
            (1, 1, 0, 0, 1, 1, 0)
        """
        E = self.encoder(encoder_name, *args, **kwargs)
        return E.encode(word)

    @cached_method
    def encoder(self, encoder_name=None, *args, **kwargs):
        r"""
        Return an encoder of ``self``.

        The returned encoder provided by this method is cached.

        This methods creates a new instance of the encoder subclass designated by ``encoder_name``.
        While it is also possible to do the same by directly calling the subclass' constructor,
        it is strongly advised to use this method to take advantage of the caching mechanism.

        INPUT:

        - ``encoder_name`` -- (default: ``None``) name of the encoder which will be
          returned. The default encoder of ``self`` will be used if
          default value is kept.

        - ``args``, ``kwargs`` -- all additional arguments are forwarded to the
          constructor of the encoder this method will return

        OUTPUT: an Encoder object

        .. NOTE::

            The default encoder always has `F^{k}` as message space, with `k` the dimension
            of ``self`` and `F` the base ring of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0], [1,0,0,1,1,0,0],
            ....:                    [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C.encoder()
            Generator matrix-based encoder for [7, 4] linear code over GF(2)

        If there is no encoder for the code, we return an error::

            sage: from sage.coding.abstract_code import AbstractCode
            sage: class MyCodeFamily(AbstractCode):
            ....:   def __init__(self, length, field):
            ....:       sage.coding.abstract_code.AbstractCode.__init__(self, length)
            ....:       Parent.__init__(self, base=field, facade=False, category=Sets())
            ....:       self._field = field
            ....:   def field(self):
            ....:       return self._field
            ....:   def _repr_(self):
            ....:       return "%d dummy code over GF(%s)" % (self.length(),
            ....:                                             self.field().cardinality())
            sage: D = MyCodeFamily(5, GF(2))
            sage: D.encoder()
            Traceback (most recent call last):
            ...
            NotImplementedError: No encoder implemented for this code.

        We check that the returned encoder is cached::

            sage: C.encoder.is_in_cache()
            True

        If the name of an encoder which is not known by ``self`` is passed,
        an exception will be raised::

            sage: sorted(C.encoders_available())
            ['GeneratorMatrix', 'Systematic']
            sage: C.encoder('NonExistingEncoder')
            Traceback (most recent call last):
            ...
            ValueError: There is no Encoder named 'NonExistingEncoder'.
            The known Encoders are: ['GeneratorMatrix', 'Systematic']

        Some encoders take extra arguments. If the user incorrectly supplies
        these, the error message attempts to be helpful::

            sage: C.encoder('Systematic', strange_parameter=True)
            Traceback (most recent call last):
            ...
            ValueError: Constructing the Systematic encoder failed,
            possibly due to missing or incorrect parameters.
            The constructor requires no arguments. It takes the optional
            arguments ['systematic_positions']. See the documentation of
            sage.coding.linear_code_no_metric.LinearCodeSystematicEncoder
            for more details.
        """
        if not self._default_encoder_name:
            raise NotImplementedError("No encoder implemented for this code.")
        if encoder_name is None:
            encoder_name = self._default_encoder_name
        if encoder_name in self._registered_encoders:
            encClass = self._registered_encoders[encoder_name]
            try:
                return encClass(self, *args, **kwargs)
            except TypeError:
                raise ValueError(
                        "Constructing the {0} encoder failed, possibly due "
                        "to missing or incorrect parameters.\n{1}".format(
                            encoder_name, _explain_constructor(encClass)))
        else:
            raise ValueError(
                    "There is no Encoder named '{0}'. The known Encoders are: "
                    "{1}".format(encoder_name, self.encoders_available()))

    def encoders_available(self, classes=False):
        r"""
        Return a list of the available encoders' names for ``self``.

        INPUT:

        - ``classes`` -- boolean (default: ``False``); if ``classes`` is set to
          ``True``, return instead a :class:`dict` mapping available encoder
          name to the associated encoder class

        OUTPUT: list of strings, or a :class:`dict` mapping strings to classes

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0], [1,0,0,1,1,0,0],
            ....:                    [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C.encoders_available()
            ['GeneratorMatrix', 'Systematic']
            sage: dictionary = C.encoders_available(True)
            sage: sorted(dictionary.items())
            [('GeneratorMatrix', <class 'sage.coding.linear_code.LinearCodeGeneratorMatrixEncoder'>),
             ('Systematic', <class 'sage.coding.linear_code_no_metric.LinearCodeSystematicEncoder'>)]
        """
        if classes:
            return copy(self._registered_encoders)

        return sorted(self._registered_encoders)

    def unencode(self, c, encoder_name=None, nocheck=False, **kwargs):
        r"""
        Return the message corresponding to ``c``.

        This is the inverse of :meth:`encode`.

        INPUT:

        - ``c`` -- a codeword of ``self``

        - ``encoder_name`` -- (default: ``None``) name of the decoder which will be used
          to decode ``word``. The default decoder of ``self`` will be used if
          default value is kept.

        - ``nocheck`` -- boolean (default: ``False``); checks if ``c`` is in
          ``self``. You might set this to ``True`` to disable the check for
          saving computation. Note that if ``c`` is not in ``self`` and
          ``nocheck = True``, then the output of :meth:`unencode` is not
          defined (except that it will be in the message space of ``self``).

        - ``kwargs`` -- all additional arguments are forwarded to the construction of the
          encoder that is used

        OUTPUT: an element of the message space of ``encoder_name`` of ``self``

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0], [1,0,0,1,1,0,0],
            ....:                    [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: c = vector(GF(2), (1, 1, 0, 0, 1, 1, 0))
            sage: C.unencode(c)
            (0, 1, 1, 0)
        """
        E = self.encoder(encoder_name, **kwargs)
        return E.unencode(c, nocheck)

    def random_element(self, *args, **kwds):
        """
        Return a random codeword; passes other positional and keyword
        arguments to ``random_element()`` method of vector space.

        OUTPUT: random element of the vector space of this code

        EXAMPLES::

            sage: C = codes.HammingCode(GF(4,'a'), 3)
            sage: C.random_element() # random test
            (1, 0, 0, a + 1, 1, a, a, a + 1, a + 1, 1, 1, 0, a + 1, a, 0, a, a, 0, a, a, 1)

        Passes extra positional or keyword arguments through::

            sage: C.random_element(prob=.5, distribution='1/n') # random test
            (1, 0, a, 0, 0, 0, 0, a + 1, 0, 0, 0, 0, 0, 0, 0, 0, a + 1, a + 1, 1, 0, 0)

        TESTS:

        Test that the codeword returned is immutable (see :issue:`16469`)::

            sage: c = C.random_element()
            sage: c.is_immutable()
            True

        Test that codeword returned has the same parent as any non-random codeword
        (see :issue:`19653`)::

            sage: C = codes.random_linear_code(GF(16, 'a'), 10, 4)
            sage: c1 = C.random_element()
            sage: c2 = C[1]
            sage: c1.parent() == c2.parent()
            True
        """
        E = self.encoder()
        M = E.message_space()
        m = M.random_element(*args, **kwds)
        c = E.encode(m)
        c.set_immutable()
        return c
