r"""
Parity-check code

A simple way of detecting up to one error is to use the device of adding a
parity check to ensure that the sum of the digits in a transmitted word is
even.

A parity-check code of dimension `k` over `F_q` is the set:
`\{(m_1, m_2, \dots, m_k, - \Sigma_{i=1}^k m_i) \mid (m_1, m_2, \dots, m_k) \in F_q^k\}`

REFERENCE:

- [Wel1988]_
"""

# *****************************************************************************
#       Copyright (C) 2016 Tania Richmond <tania.richmond@univ-tln.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from sage.coding.linear_code import AbstractLinearCode,\
    LinearCodeGeneratorMatrixEncoder
from sage.coding.encoder import Encoder
from sage.rings.integer import Integer
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.matrix.constructor import identity_matrix
from sage.modules.free_module_element import vector
from sage.modules.free_module import VectorSpace


class ParityCheckCode(AbstractLinearCode):
    r"""
    Representation of a parity-check code.

    INPUT:

    - ``base_field`` -- the base field over which ``self`` is defined.

    - ``dimension`` -- the dimension of ``self``.

    EXAMPLES::

        sage: C = codes.ParityCheckCode(GF(5), 7)
        sage: C
        [8, 7] parity-check code over GF(5)
    """

    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, base_field=GF(2), dimension=7):
        r"""
        Initialize mandatory parameters for a parity-check code object.

        INPUT:

        - ``base_field`` -- the base field over which ``self`` is defined
                            or GF(2) if no base_field.

        - ``dimension`` -- the dimension of ``self`` or 7 if no dimension.

        EXAMPLES::

            sage: C = codes.ParityCheckCode(GF(5), 7)
            sage: C
            [8, 7] parity-check code over GF(5)
        """

        if not base_field.is_finite():
            raise ValueError("base_field has to be a finite field")
        if not isinstance(dimension, (Integer, int)):
            raise ValueError("dimension must be an integer")

        self._dimension = dimension
        super(ParityCheckCode, self).__init__(
            base_field, dimension + 1, "ParityCheckCodeGeneratorMatrixEncoder",
            "Syndrome")

    def __eq__(self, other):
        r"""
        Test equality of parity-check code objects.

        EXAMPLES::

            sage: C1 = codes.ParityCheckCode(GF(5), 7)
            sage: C2 = codes.ParityCheckCode(GF(5), 7)
            sage: C1 == C2
            True
        """
        return (isinstance(other, ParityCheckCode)
                and self.base_field() == other.base_field()
                and self.dimension() == other.dimension())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.ParityCheckCode(GF(5), 7)
            sage: C
            [8, 7] parity-check code over GF(5)
        """
        return ("[%s, %s] parity-check code over GF(%s)"
                % (self.length(), self.dimension(),
                   self.base_field().cardinality()))

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: C = codes.ParityCheckCode(GF(5), 7)
            sage: latex(C)
            [8, 7] \textnormal{parity-check code over } \Bold{F}_{5}
        """
        return "[%s, %s] \\textnormal{parity-check code over } %s"\
               % (self.length(), self.dimension(), self.base_field()._latex_())

    def minimum_distance(self):
        r"""
        Return the minimum distance of ``self``.

        It is always 2 as ``self`` is a parity-check code.

        EXAMPLES::

            sage: C = codes.ParityCheckCode(GF(5), 7)
            sage: C.minimum_distance()
            2
        """
        return 2


####################### encoders ###############################

class ParityCheckCodeGeneratorMatrixEncoder(LinearCodeGeneratorMatrixEncoder):
    r"""
    Encoder for parity-check codes which uses a generator matrix to obtain
    codewords.

    INPUT:

    - ``code`` -- the associated code of this encoder.

    EXAMPLES::

        sage: C = codes.ParityCheckCode(GF(5), 7)
        sage: E = codes.encoders.ParityCheckCodeGeneratorMatrixEncoder(C)
        sage: E
        Generator matrix-based encoder for [8, 7] parity-check code over GF(5)

    Actually, we can construct the encoder from ``C`` directly::

        sage: E = C.encoder("ParityCheckCodeGeneratorMatrixEncoder")
        sage: E
        Generator matrix-based encoder for [8, 7] parity-check code over GF(5)
    """

    def __init__(self, code):
        r"""
        TESTS:

        If ``code`` is not a parity-check code, an error is raised::

            sage: C = codes.random_linear_code(GF(11), 10, 4)
            sage: codes.encoders.ParityCheckCodeGeneratorMatrixEncoder(C)
            Traceback (most recent call last):
            ...
            ValueError: code has to be a parity-check code
        """
        if not isinstance(code, ParityCheckCode):
            raise ValueError("code has to be a parity-check code")
        super(ParityCheckCodeGeneratorMatrixEncoder, self).__init__(code)

    def generator_matrix(self):
        r"""
        Return a generator matrix of ``self``.

        EXAMPLES::

            sage: C = codes.ParityCheckCode(GF(5),7)
            sage: E = codes.encoders.ParityCheckCodeGeneratorMatrixEncoder(C)
            sage: E.generator_matrix()
            [1 0 0 0 0 0 0 4]
            [0 1 0 0 0 0 0 4]
            [0 0 1 0 0 0 0 4]
            [0 0 0 1 0 0 0 4]
            [0 0 0 0 1 0 0 4]
            [0 0 0 0 0 1 0 4]
            [0 0 0 0 0 0 1 4]
        """
        k = self.code().dimension()
        field = self.code().base_field()
        G = identity_matrix(field, k)
        G = G.augment(vector(field, [-field.one()] * k))
        G.set_immutable()
        return G


class ParityCheckCodeStraightforwardEncoder(Encoder):
    r"""
    Encoder for parity-check codes which computes the sum of message symbols
    and appends its opposite to the message to obtain codewords.

    INPUT:

    - ``code`` -- the associated code of this encoder.

    EXAMPLES::

        sage: C = codes.ParityCheckCode(GF(5), 7)
        sage: E = codes.encoders.ParityCheckCodeStraightforwardEncoder(C)
        sage: E
        Parity-check encoder for the [8, 7] parity-check code over GF(5)

    Actually, we can construct the encoder from ``C`` directly::

        sage: E = C.encoder("ParityCheckCodeStraightforwardEncoder")
        sage: E
        Parity-check encoder for the [8, 7] parity-check code over GF(5)
    """

    def __init__(self, code):
        r"""
        TESTS:

        If ``code`` is not a parity-check code, an error is raised::

            sage: C  = codes.random_linear_code(GF(11), 10, 4)
            sage: codes.encoders.ParityCheckCodeStraightforwardEncoder(C)
            Traceback (most recent call last):
            ...
            ValueError: code has to be a parity-check code
        """
        if not isinstance(code, ParityCheckCode):
            raise ValueError("code has to be a parity-check code")
        super(ParityCheckCodeStraightforwardEncoder, self).__init__(code)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.ParityCheckCode(GF(5), 7)
            sage: codes.encoders.ParityCheckCodeStraightforwardEncoder(C)
            Parity-check encoder for the [8, 7] parity-check code over GF(5)
        """
        return "Parity-check encoder for the %s" % self.code()

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: C = codes.ParityCheckCode(GF(5), 7)
            sage: E = codes.encoders.ParityCheckCodeStraightforwardEncoder(C)
            sage: latex(E)
            \textnormal{Parity-check encoder for the } [8, 7] \textnormal{parity-check code over } \Bold{F}_{5}
        """
        return ("\\textnormal{Parity-check encoder for the } %s" %
                self.code()._latex_())

    def __eq__(self, other):
        r"""
        Test equality of parity-check code objects.

        EXAMPLES::

            sage: C1 = codes.ParityCheckCode(GF(5), 7)
            sage: C2 = codes.ParityCheckCode(GF(5), 7)
            sage: C1 == C2
            True
        """
        return (isinstance(other, ParityCheckCodeStraightforwardEncoder)
                and self.code() == other.code())

    def encode(self, message):
        r"""
        Transform the vector ``message`` into a codeword of :meth:`code`.

        INPUT:

        - ``message`` -- A ``self.code().dimension()``-vector from the message
          space of ``self``.

        OUTPUT:

        - A codeword in associated code of ``self``.

        EXAMPLES::

            sage: C = codes.ParityCheckCode(GF(5),7)
            sage: message = vector(C.base_field(),[1,0,4,2,0,3,2])
            sage: C.encode(message)
            (1, 0, 4, 2, 0, 3, 2, 3)
        """
        parity=self.code().base_field().zero()
        for i in message.list():
            parity += i
        return vector(self.code().base_field(), message.list() + [-parity])

    def unencode_nocheck(self, word):
        r"""
        Return the message corresponding to the vector ``word``.

        Use this method with caution: it does not check if ``word``
        belongs to the code.

        INPUT:

        - ``word`` -- A ``self.code().length()``-vector from the ambiant space
          of ``self``.

        OUTPUT:

        - A vector corresponding to the ``self.code().dimension()``-first
          symbols in ``word``.

        EXAMPLES::

            sage: C = codes.ParityCheckCode(GF(5), 7)
            sage: word = vector(C.base_field(), [1, 0, 4, 2, 0, 3, 2, 3])
            sage: E = codes.encoders.ParityCheckCodeStraightforwardEncoder(C)
            sage: E.unencode_nocheck(word)
            (1, 0, 4, 2, 0, 3, 2)
        """
        return word[:-1]

    def message_space(self):
        r"""
        Return the message space of ``self``.

        EXAMPLES::

            sage: C = codes.ParityCheckCode(GF(5),7)
            sage: E = codes.encoders.ParityCheckCodeStraightforwardEncoder(C)
            sage: E.message_space()
            Vector space of dimension 7 over Finite Field of size 5
        """
        return VectorSpace(self.code().base_field(), self.code().dimension())

####################### registration ###############################

ParityCheckCode._registered_encoders["ParityCheckCodeGeneratorMatrixEncoder"] = ParityCheckCodeGeneratorMatrixEncoder
ParityCheckCode._registered_encoders["ParityCheckCodeStraightforwardEncoder"] = ParityCheckCodeStraightforwardEncoder
