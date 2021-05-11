r"""
Extended code

Let `C` be a linear code of length `n` over `\GF{q}`. The extended code of `C` is the code

.. MATH::

    \hat{C} = \{x_{1}x_{2}\dots x_{n+1} \in \GF{q}^{n+1} \,\vert\,  x_{1}x_{2}\dots x_{n} \in C \text{ with } x_{1} + x_{2} + \dots + x_{n+1} = 0 \}.

See [HP2003]_ (pp 15-16) for details.
"""

#*****************************************************************************
#       Copyright (C) 2016 David Lucas <david.lucas@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .linear_code import (AbstractLinearCode,
        LinearCodeSyndromeDecoder,
        LinearCodeNearestNeighborDecoder)
from .encoder import Encoder
from .decoder import Decoder
from sage.misc.cachefunc import cached_method
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from copy import copy

class ExtendedCode(AbstractLinearCode):
    r"""
    Representation of an extended code.

    INPUT:

    -  ``C`` -- A linear code

    EXAMPLES::

        sage: C = codes.random_linear_code(GF(7), 11, 5)
        sage: Ce = codes.ExtendedCode(C)
        sage: Ce
        Extension of [11, 5] linear code over GF(7)
    """

    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, C):
        r"""
        TESTS:

        ``C`` must be a linear code::

            sage: C = VectorSpace(GF(7), 11)
            sage: codes.ExtendedCode(C)
            Traceback (most recent call last):
            ...
            ValueError: Provided code must be a linear code
        """
        if not isinstance(C, AbstractLinearCode):
            raise ValueError("Provided code must be a linear code")
        super(ExtendedCode, self).__init__(C.base_ring(), C.length() + 1, "ExtendedMatrix", "OriginalDecoder")
        self._original_code = C
        self._dimension = C.dimension()

    def __eq__(self, other):
        r"""
        Tests equality between two extended codes.

        EXAMPLES::

            sage: C = codes.random_linear_code(GF(7), 11, 5)
            sage: C1 = codes.ExtendedCode(C)
            sage: C2 = codes.ExtendedCode(C)
            sage: C1 == C2
            True
        """
        return isinstance(other, ExtendedCode)\
                and self.original_code() == other.original_code()

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.random_linear_code(GF(7), 11, 5)
            sage: Ce = codes.ExtendedCode(C)
            sage: Ce
            Extension of [11, 5] linear code over GF(7)
        """
        return "Extension of %s" % self.original_code()

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: C = codes.random_linear_code(GF(7), 11, 5)
            sage: Ce = codes.ExtendedCode(C)
            sage: latex(Ce)
            \textnormal{Extension of [11, 5] linear code over GF(7)}
        """
        return "\\textnormal{Extension of %s}" % self.original_code()

    def original_code(self):
        r"""
        Returns the code which was extended to get ``self``.

        EXAMPLES::

            sage: C = codes.random_linear_code(GF(7), 11, 5)
            sage: Ce = codes.ExtendedCode(C)
            sage: Ce.original_code()
            [11, 5] linear code over GF(7)
        """
        return self._original_code

    @cached_method
    def parity_check_matrix(self):
        r"""
        Returns a parity check matrix of ``self``.

        This matrix is computed directly from :func:`original_code`.

        EXAMPLES::

            sage: C = LinearCode(matrix(GF(2),[[1,0,0,1,1],\
                                               [0,1,0,1,0],\
                                               [0,0,1,1,1]]))
            sage: C.parity_check_matrix()
            [1 0 1 0 1]
            [0 1 0 1 1]
            sage: Ce = codes.ExtendedCode(C)
            sage: Ce.parity_check_matrix()
            [1 1 1 1 1 1]
            [1 0 1 0 1 0]
            [0 1 0 1 1 0]
        """
        F = self.base_ring()
        zero = F.zero()
        one = F.one()
        H = self.original_code().parity_check_matrix()
        nr, nc = H.nrows(), H.ncols()
        Hlist = H.list()
        v = matrix(F, nr + 1, 1, [one] + [zero] * nr)
        M = matrix(F, nr + 1, nc, [one] * nc + Hlist).augment(v)
        M.set_immutable()
        return M

    def random_element(self):
        r"""
        Returns a random element of ``self``.

        This random element is computed directly from the original code,
        and does not compute a generator matrix of ``self`` in the process.

        EXAMPLES::

            sage: C = codes.random_linear_code(GF(7), 9, 5)
            sage: Ce = codes.ExtendedCode(C)
            sage: c = Ce.random_element() #random
            sage: c in Ce
            True
        """
        c = self.original_code().random_element()
        c_list = c.list()
        F = self.base_ring()
        last_element = F.zero()
        for i in c_list:
            last_element += i
        c_list.append(-last_element)
        return vector(F, c_list)










class ExtendedCodeExtendedMatrixEncoder(Encoder):
    r"""
    Encoder using original code's generator matrix to compute the extended code's one.

    INPUT:

    - ``code`` -- The associated code of ``self``.
    """

    def __init__(self, code):
        r"""
        EXAMPLES::

            sage: C = codes.random_linear_code(GF(7), 11, 5)
            sage: Ce = codes.ExtendedCode(C)
            sage: E = codes.encoders.ExtendedCodeExtendedMatrixEncoder(Ce)
            sage: E
            Matrix-based Encoder for Extension of [11, 5] linear code over GF(7)
        """
        if not isinstance(code, ExtendedCode):
            raise TypeError("code has to be an instance of ExtendedCode class")

        super(ExtendedCodeExtendedMatrixEncoder, self).__init__(code)

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.random_linear_code(GF(7), 11, 5)
            sage: Ce = codes.ExtendedCode(C)
            sage: E = codes.encoders.ExtendedCodeExtendedMatrixEncoder(Ce)
            sage: E
            Matrix-based Encoder for Extension of [11, 5] linear code over GF(7)
        """
        return "Matrix-based Encoder for %s" % self.code()

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: C = codes.random_linear_code(GF(7), 11, 5)
            sage: Ce = codes.ExtendedCode(C)
            sage: E = codes.encoders.ExtendedCodeExtendedMatrixEncoder(Ce)
            sage: latex(E)
            \textnormal{Matrix-based Encoder for }\textnormal{Extension of [11, 5] linear code over GF(7)}
        """
        return "\\textnormal{Matrix-based Encoder for }%s" % self.code()._latex_()

    def __eq__(self, other):
        r"""
        Tests equality between GRSEvaluationVectorEncoder objects.

        EXAMPLES::

            sage: C = codes.random_linear_code(GF(7), 11, 5)
            sage: Ce = codes.ExtendedCode(C)
            sage: D1 = codes.encoders.ExtendedCodeExtendedMatrixEncoder(Ce)
            sage: D2 = codes.encoders.ExtendedCodeExtendedMatrixEncoder(Ce)
            sage: D1.__eq__(D2)
            True
            sage: D1 is D2
            False
        """
        return isinstance(other, ExtendedCodeExtendedMatrixEncoder) \
                and self.code() == other.code()

    @cached_method
    def generator_matrix(self):
        r"""
        Returns a generator matrix of the associated code of ``self``.

        EXAMPLES::

            sage: C = LinearCode(matrix(GF(2),[[1,0,0,1,1],\
                                               [0,1,0,1,0],\
                                               [0,0,1,1,1]]))
            sage: Ce = codes.ExtendedCode(C)
            sage: E = codes.encoders.ExtendedCodeExtendedMatrixEncoder(Ce)
            sage: E.generator_matrix()
            [1 0 0 1 1 1]
            [0 1 0 1 0 0]
            [0 0 1 1 1 1]
        """
        C = self.code()
        F = C.base_ring()
        Cor = C.original_code()
        G = Cor.generator_matrix()
        k = C.dimension()
        extra_col = [-sum(G.rows()[i]) for i in range(k)]
        extra_col = matrix(F, k, 1, extra_col)
        M = G.augment(extra_col)
        M.set_immutable()
        return M










class ExtendedCodeOriginalCodeDecoder(Decoder):
    r"""
    Decoder which decodes through a decoder over the original code.

    INPUT:

    - ``code`` -- The associated code of this decoder

    - ``original_decoder`` -- (default: ``None``) the decoder that will be used over the original code.
      It has to be a decoder object over the original code.
      If ``original_decoder`` is set to ``None``, it will use the default decoder of the original code.

    - ``**kwargs`` -- all extra arguments are forwarded to original code's decoder

    EXAMPLES::

        sage: C = codes.GeneralizedReedSolomonCode(GF(16, 'a').list()[:15], 7)
        sage: Ce = codes.ExtendedCode(C)
        sage: D = codes.decoders.ExtendedCodeOriginalCodeDecoder(Ce)
        sage: D
        Decoder of Extension of [15, 7, 9] Reed-Solomon Code over GF(16) through Gao decoder for [15, 7, 9] Reed-Solomon Code over GF(16)
    """

    def __init__(self, code, original_decoder = None, **kwargs):
        r"""
        TESTS:

        If one tries to use a decoder whose code is not the original code, it returns an error::

            sage: C1 = codes.GeneralizedReedSolomonCode(GF(16, 'a').list()[:15], 7)
            sage: Ce = codes.ExtendedCode(C1)
            sage: C2 = codes.GeneralizedReedSolomonCode(GF(13).list()[:12], 7)
            sage: Dc2 = C2.decoder()
            sage: D = codes.decoders.ExtendedCodeOriginalCodeDecoder(Ce, original_decoder = Dc2)
            Traceback (most recent call last):
            ...
            ValueError: Original decoder must have the original code as associated code
        """
        if not isinstance(code, ExtendedCode):
            raise TypeError("code has to be an instance of ExtendedCode class")

        original_code = code.original_code()
        if original_decoder is not None and not original_decoder.code() == original_code:
            raise ValueError("Original decoder must have the original code as associated code")
        elif original_decoder is None:
            self._original_decoder = original_code.decoder()
        else:
            self._original_decoder = original_decoder
        self._decoder_type = copy(self._decoder_type)
        self._decoder_type.remove("dynamic")
        self._decoder_type = self._original_decoder.decoder_type()
        super(ExtendedCodeOriginalCodeDecoder, self).__init__(code, code.ambient_space(),\
                self._original_decoder.connected_encoder())

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(16, 'a').list()[:15], 7)
            sage: Ce = codes.ExtendedCode(C)
            sage: D = codes.decoders.ExtendedCodeOriginalCodeDecoder(Ce)
            sage: D
            Decoder of Extension of [15, 7, 9] Reed-Solomon Code over GF(16) through Gao decoder for [15, 7, 9] Reed-Solomon Code over GF(16)
        """
        return "Decoder of %s through %s" % (self.code(), self.original_decoder())

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(16, 'a').list()[:15], 7)
            sage: Ce = codes.ExtendedCode(C)
            sage: D = codes.decoders.ExtendedCodeOriginalCodeDecoder(Ce)
            sage: latex(D)
            \textnormal{Decoder of } Extension of [15, 7, 9] Reed-Solomon Code over GF(16) \textnormal{ through } Gao decoder for [15, 7, 9] Reed-Solomon Code over GF(16)
        """
        return "\\textnormal{Decoder of } %s \\textnormal{ through } %s" % (self.code(), self.original_decoder())

    def original_decoder(self):
        r"""
        Returns the decoder over the original code that will be used to decode words of
        :meth:`sage.coding.decoder.Decoder.code`.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(16, 'a').list()[:15], 7)
            sage: Ce = codes.ExtendedCode(C)
            sage: D = codes.decoders.ExtendedCodeOriginalCodeDecoder(Ce)
            sage: D.original_decoder()
            Gao decoder for [15, 7, 9] Reed-Solomon Code over GF(16)
        """
        return self._original_decoder

    def decode_to_code(self, y, **kwargs):
        r"""
        Decodes ``y`` to an element in :meth:`sage.coding.decoder.Decoder.code`.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(16, 'a').list()[:15], 7)
            sage: Ce = codes.ExtendedCode(C)
            sage: D = codes.decoders.ExtendedCodeOriginalCodeDecoder(Ce)
            sage: c = Ce.random_element()
            sage: Chan = channels.StaticErrorRateChannel(Ce.ambient_space(), D.decoding_radius())
            sage: y = Chan(c)
            sage: y in Ce
            False
            sage: D.decode_to_code(y) == c
            True

        Another example, with a list decoder::

            sage: C = codes.GeneralizedReedSolomonCode(GF(16, 'a').list()[:15], 7)
            sage: Ce = codes.ExtendedCode(C)
            sage: Dgrs = C.decoder('GuruswamiSudan', tau = 4)
            sage: D = codes.decoders.ExtendedCodeOriginalCodeDecoder(Ce, original_decoder = Dgrs)
            sage: c = Ce.random_element()
            sage: Chan = channels.StaticErrorRateChannel(Ce.ambient_space(), D.decoding_radius())
            sage: y = Chan(c)
            sage: y in Ce
            False
            sage: c in D.decode_to_code(y)
            True
        """
        D = self.original_decoder()
        C = self.code()
        F = C.base_field()
        n = C.length()
        y_original = copy(y.list())
        y_original.pop(n - 1)
        decoded = D.decode_to_code(vector(y_original), **kwargs)
        if 'list-decoder' in self.decoder_type():
            l = []
            for word in decoded:
                last_pos = F.zero()
                for i in word:
                    last_pos += i
                word_list = list(word)
                word_list.append(last_pos)
                l.append(vector(F, word_list))
            return l
        else:
            last_pos = F.zero()
            for i in decoded:
                last_pos += i
            decoded_list = list(decoded)
            decoded_list.append(last_pos)
            return vector(F, decoded_list)

    def decoding_radius(self, *args, **kwargs):
        r"""
        Returns maximal number of errors that ``self`` can decode.

        INPUT:

        - ``*args``, ``**kwargs`` -- arguments and optional arguments are
          forwarded to original decoder's ``decoding_radius`` method.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(16, 'a').list()[:15], 7)
            sage: Ce = codes.ExtendedCode(C)
            sage: D = codes.decoders.ExtendedCodeOriginalCodeDecoder(Ce)
            sage: D.decoding_radius()
            4
        """
        return self.original_decoder().decoding_radius(*args, **kwargs)


####################### registration ###############################

ExtendedCode._registered_encoders["ExtendedMatrix"] = ExtendedCodeExtendedMatrixEncoder
ExtendedCode._registered_decoders["Syndrome"] = LinearCodeSyndromeDecoder
ExtendedCode._registered_decoders["NearestNeighbor"] = LinearCodeNearestNeighborDecoder
ExtendedCode._registered_decoders["OriginalDecoder"] = ExtendedCodeOriginalCodeDecoder
ExtendedCodeOriginalCodeDecoder._decoder_type = {"dynamic"}
