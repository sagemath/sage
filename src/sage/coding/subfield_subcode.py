# sage.doctest: needs sage.modules sage.rings.finite_rings
r"""
Subfield subcode

Let `C` be a `[n, k]` code over `\GF{q^t}`.
Let `Cs = \{c \in C | \forall i, c_i \in \GF{q}\}`, `c_i` being the `i`-th
coordinate of `c`.

`Cs` is called the subfield subcode of `C` over `\GF{q}`
"""

#*****************************************************************************
#       Copyright (C) 2016 David Lucas, Inria  <david.lucas@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .linear_code import AbstractLinearCode
from sage.misc.cachefunc import cached_method
from sage.categories.homset import Hom
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from .decoder import Decoder, DecodingError
from copy import copy


class SubfieldSubcode(AbstractLinearCode):
    r"""
    Representation of a subfield subcode.

    INPUT:

    - ``original_code`` -- the code ``self`` comes from

    - ``subfield`` -- the base field of ``self``

    - ``embedding`` -- (default: ``None``) a homomorphism from ``subfield`` to
      ``original_code``'s base field. If ``None`` is provided, it will default
      to the first homomorphism of the list of homomorphisms Sage can build.

    EXAMPLES::

        sage: C = codes.random_linear_code(GF(16, 'aa'), 7, 3)
        sage: codes.SubfieldSubcode(C, GF(4, 'a'))
        Subfield subcode of [7, 3] linear code over GF(16) down to GF(4)
    """
    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, original_code, subfield, embedding=None):
        r"""
        TESTS:

        ``subfield`` has to be a finite field, otherwise an error is raised::

            sage: C = codes.random_linear_code(GF(16, 'aa'), 7, 3)
            sage: Cs = codes.SubfieldSubcode(C, RR)
            Traceback (most recent call last):
            ...
            ValueError: subfield has to be a finite field

        ``subfield`` has to be a subfield of ``original_code``'s base field,
        otherwise an error is raised::

            sage: C = codes.random_linear_code(GF(16, 'aa'), 7, 3)
            sage: Cs = codes.SubfieldSubcode(C, GF(8, 'a'))
            Traceback (most recent call last):
            ...
            ValueError: subfield has to be a subfield of the base field of the original code
        """
        if not isinstance(original_code, AbstractLinearCode):
            raise ValueError("original_code must be a linear code")
        if not subfield.is_finite():
            raise ValueError("subfield has to be a finite field")

        super().__init__(subfield, original_code.length(), "Systematic", "Syndrome")

        F = original_code.base_field()
        sm = F.degree()
        s = subfield.degree()
        if not s.divides(sm):
            raise ValueError("subfield has to be a subfield of the base field of the original code")

        H = Hom(subfield, F)
        if embedding is not None and embedding not in H:
            raise ValueError("embedding has to be an embedding from subfield to original_code's base field")
        if embedding is None:
            embedding = H[0]

        self._extension_degree = sm // s
        self._embedding = embedding
        self._original_code = original_code

    def __eq__(self, other):
        r"""
        Test equality between Subfield Subcode objects.

        EXAMPLES::

            sage: C = codes.random_linear_code(GF(16, 'aa'), 7, 3)
            sage: Cs1 = codes.SubfieldSubcode(C, GF(4, 'a'))
            sage: Cs2 = codes.SubfieldSubcode(C, GF(4, 'a'))
            sage: Cs1 == Cs2
            True
        """
        return isinstance(other, SubfieldSubcode) \
                and self.original_code() == other.original_code()\
                and self.embedding() == other.embedding()

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.random_linear_code(GF(16, 'aa'), 7, 3)
            sage: Cs = codes.SubfieldSubcode(C, GF(4, 'a'))
            sage: Cs
            Subfield subcode of [7, 3] linear code over GF(16) down to GF(4)
        """
        return "Subfield subcode of %s down to GF(%s)"\
                % (self.original_code(), self.base_field().cardinality())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: C = codes.random_linear_code(GF(16, 'aa'), 7, 3)
            sage: Cs = codes.SubfieldSubcode(C, GF(4, 'a'))
            sage: latex(Cs)
            \textnormal{Subfield subcode of }[7, 3]\textnormal{ Linear code over }\Bold{F}_{2^{4}}\textnormal{ down to }\Bold{F}_{2^{2}}
        """
        return "\\textnormal{Subfield subcode of }%s\\textnormal{ down to }%s"\
                % (self.original_code()._latex_(), self.base_field()._latex_())

    def dimension(self):
        r"""
        Return the dimension of ``self``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(16, 'aa').list()[:13], 5)
            sage: Cs = codes.SubfieldSubcode(C, GF(4, 'a'))
            sage: Cs.dimension()
            3
        """
        return self.generator_matrix().nrows()

    def dimension_upper_bound(self):
        r"""
        Return an upper bound for the dimension of ``self``.

        EXAMPLES::

            sage: C = codes.random_linear_code(GF(16, 'aa'), 7, 3)
            sage: Cs = codes.SubfieldSubcode(C, GF(4, 'a'))
            sage: Cs.dimension_upper_bound()
            3
        """
        return self.original_code().dimension()

    def dimension_lower_bound(self):
        r"""
        Return a lower bound for the dimension of ``self``.

        EXAMPLES::

            sage: C = codes.random_linear_code(GF(16, 'aa'), 7, 3)
            sage: Cs = codes.SubfieldSubcode(C, GF(4, 'a'))
            sage: Cs.dimension_lower_bound()
            -1
        """
        C = self.original_code()
        n = C.length()
        k = C.dimension()
        m = self._extension_degree
        return n - m*(n-k)

    def original_code(self):
        r"""
        Return the original code of ``self``.

        EXAMPLES::

            sage: C = codes.random_linear_code(GF(16, 'aa'), 7, 3)
            sage: Cs = codes.SubfieldSubcode(C, GF(4, 'a'))
            sage: Cs.original_code()
            [7, 3] linear code over GF(16)
        """
        return self._original_code

    def embedding(self):
        r"""
        Return the field embedding between the base field of ``self`` and
        the base field of its original code.

        EXAMPLES::

            sage: C = codes.random_linear_code(GF(16, 'aa'), 7, 3)
            sage: Cs = codes.SubfieldSubcode(C, GF(4, 'a'))
            sage: Cs.embedding()
            Ring morphism:
              From: Finite Field in a of size 2^2
              To:   Finite Field in aa of size 2^4
              Defn: a |--> aa^2 + aa
        """
        return self._embedding

    @cached_method
    def parity_check_matrix(self):
        r"""
        Return a parity check matrix of ``self``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(16, 'aa').list()[:13], 5)
            sage: Cs = codes.SubfieldSubcode(C, GF(4, 'a'))
            sage: Cs.parity_check_matrix()
            [    1     0     0     0     0     0     0     0     0     0     1 a + 1 a + 1]
            [    0     1     0     0     0     0     0     0     0     0 a + 1     0     a]
            [    0     0     1     0     0     0     0     0     0     0 a + 1     a     0]
            [    0     0     0     1     0     0     0     0     0     0     0 a + 1     a]
            [    0     0     0     0     1     0     0     0     0     0 a + 1     1 a + 1]
            [    0     0     0     0     0     1     0     0     0     0     1     1     1]
            [    0     0     0     0     0     0     1     0     0     0     a     a     1]
            [    0     0     0     0     0     0     0     1     0     0     a     1     a]
            [    0     0     0     0     0     0     0     0     1     0 a + 1 a + 1     1]
            [    0     0     0     0     0     0     0     0     0     1     a     0 a + 1]
        """
        F = self.base_field()
        n = self.length()

        C = self.original_code()
        H_original = C.parity_check_matrix()
        codimC = H_original.nrows()

        phi = self.embedding()
        E = phi.codomain()
        V, from_V, to_V = E.vector_space(phi, map=True)

        m = V.dimension()
        H = matrix(F, codimC * m, n)
        for i in range(codimC):
            for j in range(n):
                h_vec = to_V(H_original[i][j])
                for k in range(m):
                    H[i*m+k, j] = h_vec[k]

        H = H.echelon_form()
        delete = [i for i in range(H.nrows()) if H.row(i) == 0]
        M = H.delete_rows(delete)
        M.set_immutable()
        return M


class SubfieldSubcodeOriginalCodeDecoder(Decoder):
    r"""
    Decoder decoding through a decoder over the original code of ``code``.

    INPUT:

    - ``code`` -- the associated code of this decoder

    - ``original_decoder`` -- (default: ``None``) the decoder that will be used
      over the original code. It has to be a decoder object over the original
      code. If it is set to ``None``, the default decoder over the original
      code will be used.

    - ``**kwargs`` -- all extra arguments are forwarded to original code's decoder

    EXAMPLES::

        sage: C = codes.GeneralizedReedSolomonCode(GF(16, 'aa').list()[:13], 5)
        sage: Cs = codes.SubfieldSubcode(C, GF(4, 'a'))
        sage: codes.decoders.SubfieldSubcodeOriginalCodeDecoder(Cs)
        Decoder of Subfield subcode of [13, 5, 9] Reed-Solomon Code over GF(16) down to GF(4)
         through Gao decoder for [13, 5, 9] Reed-Solomon Code over GF(16)
    """

    def __init__(self, code, original_decoder=None, **kwargs):
        r"""
        TESTS:

        If the original decoder is not a decoder over ``code``'s original code, an error is
        raised::

            sage: C = codes.GeneralizedReedSolomonCode(GF(16, 'aa').list()[:13], 5)
            sage: Cs = codes.SubfieldSubcode(C, GF(4, 'a'))
            sage: Cbis = codes.GeneralizedReedSolomonCode(GF(16, 'aa').list()[:9], 5)
            sage: D = Cbis.decoder()
            sage: codes.decoders.SubfieldSubcodeOriginalCodeDecoder(Cs, original_decoder=D)
            Traceback (most recent call last):
            ...
            ValueError: original_decoder must have the original code as associated code
        """
        original_code = code.original_code()
        if original_decoder is not None and not original_decoder.code() == code.original_code():
            raise ValueError("original_decoder must have the original code as associated code")
        elif original_decoder is not None:
            self._original_decoder = original_decoder
        else:
            self._original_decoder = original_code.decoder(**kwargs)

        super().__init__(code, code.ambient_space(),
                         self._original_decoder.connected_encoder())

        self._decoder_type = copy(self._decoder_type)
        self._decoder_type.remove("dynamic")
        self._decoder_type = self._original_decoder.decoder_type()

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(16, 'aa').list()[:13], 5)
            sage: Cs = codes.SubfieldSubcode(C, GF(4, 'a'))
            sage: D = codes.decoders.SubfieldSubcodeOriginalCodeDecoder(Cs)
            sage: D
            Decoder of Subfield subcode of [13, 5, 9] Reed-Solomon Code over GF(16) down to GF(4) through Gao decoder for [13, 5, 9] Reed-Solomon Code over GF(16)
        """
        return "Decoder of %s through %s" % (self.code(), self.original_decoder())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(16, 'aa').list()[:13], 5)
            sage: Cs = codes.SubfieldSubcode(C, GF(4, 'a'))
            sage: D = codes.decoders.SubfieldSubcodeOriginalCodeDecoder(Cs)
            sage: latex(D)
            \textnormal{Decoder of Subfield subcode of [13, 5, 9] Reed-Solomon Code over GF(16) down to GF(4) through } Gao decoder for [13, 5, 9] Reed-Solomon Code over GF(16)
        """
        return "\\textnormal{Decoder of %s through } %s" % (self.code(), self.original_decoder())

    def original_decoder(self):
        r"""
        Return the decoder over the original code that will be used to decode words of
        :meth:`sage.coding.decoder.Decoder.code`.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(16, 'aa').list()[:13], 5)
            sage: Cs = codes.SubfieldSubcode(C, GF(4, 'a'))
            sage: D = codes.decoders.SubfieldSubcodeOriginalCodeDecoder(Cs)
            sage: D.original_decoder()
            Gao decoder for [13, 5, 9] Reed-Solomon Code over GF(16)
        """
        return self._original_decoder

    def decode_to_code(self, y):
        """
        Return an error-corrected codeword from ``y``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(16, 'aa').list()[:13], 5)
            sage: Cs = codes.SubfieldSubcode(C, GF(4, 'a'))
            sage: D = codes.decoders.SubfieldSubcodeOriginalCodeDecoder(Cs)
            sage: Chan = channels.StaticErrorRateChannel(Cs.ambient_space(),
            ....:                                        D.decoding_radius())
            sage: c = Cs.random_element()
            sage: y = Chan(c)
            sage: c == D.decode_to_code(y)
            True
        """
        C = self.code()
        D = self.original_decoder()
        phi = C.embedding()
        sec = phi.section()

        result = D.decode_to_code(vector([phi(i) for i in y]))

        if 'list-decoder' in self.decoder_type():
            l = []
            for cw in result:
                try:
                    l.append(vector([sec(c) for c in cw]))
                except ValueError:  # not a codeword of this code
                    pass
            return l
        else:
            try:
                cw = vector([sec(c) for c in result])
            except ValueError:  # not a codeword of this code
                raise DecodingError("Original decoder does not output a subfield codeword. "
                                    "You may have exceeded the decoding radius.")
            return cw

    def decoding_radius(self, **kwargs):
        r"""
        Return the maximal number of errors ``self`` can decode.

        INPUT:

        - ``kwargs`` -- optional arguments are forwarded to original decoder's
          :meth:`sage.coding.decoder.Decoder.decoding_radius` method

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(16, 'aa').list()[:13], 5)
            sage: Cs = codes.SubfieldSubcode(C, GF(4, 'a'))
            sage: D = codes.decoders.SubfieldSubcodeOriginalCodeDecoder(Cs)
            sage: D.decoding_radius()
            4
        """
        return self.original_decoder().decoding_radius(**kwargs)


####################### registration ###############################

SubfieldSubcode._registered_decoders["OriginalCode"] = SubfieldSubcodeOriginalCodeDecoder
SubfieldSubcodeOriginalCodeDecoder._decoder_type = {"dynamic"}
