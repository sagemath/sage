# sage.doctest: needs sage.modules sage.rings.finite_rings
r"""
Punctured code

Let `C` be a linear code. Let `C_i` be the set of all words of `C` with the
`i`-th coordinate being removed. `C_i` is the punctured code of `C`
on the `i`-th position.
"""

#*****************************************************************************
#       Copyright (C) 2015 David Lucas <david.lucas@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .linear_code import AbstractLinearCode
from .encoder import Encoder
from .decoder import Decoder, DecodingError
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer
from sage.modules.free_module import VectorSpace
from sage.modules.free_module_element import vector
from sage.rings.finite_rings.finite_field_constructor import GF
from copy import copy


def _puncture(v, points):
    r"""
    Return v punctured as the positions listed in ``points``.

    INPUT:

    - ``v`` -- a vector or a list of vectors

    - ``points`` -- set of integers, or an integer

    EXAMPLES::

        sage: v = vector(GF(7), (2,3,0,2,1,5,1,5,6,5,3))
        sage: from sage.coding.punctured_code import _puncture
        sage: _puncture(v, {4, 3})
        (2, 3, 0, 5, 1, 5, 6, 5, 3)
    """
    if not isinstance(points, (Integer, int, set)):
        raise TypeError("points must be either a Sage Integer, a Python int, or a set")
    if isinstance(v, list):
        size = len(v[0])
        S = VectorSpace(v[0].base_ring(), size - len(points))
        l = []
        for i in v:
            new_v = [i[j] for j in range(size) if j not in points]
            l.append(S(new_v))
        return l
    S = VectorSpace(v.base_ring(), len(v) - len(points))
    new_v = [v[i] for i in range(len(v)) if i not in points]
    return S(new_v)


def _insert_punctured_positions(l, punctured_points, value=None):
    r"""
    Return ``l`` with ``value`` inserted in the corresponding
    position from ``punctured_points``.

    INPUT:

    - ``l`` -- list

    - ``punctured_points`` -- set of integers

    - ``value`` -- (default: ``None``) an element to insert in every position
      given in``punctured_points``.  If it is let to ``None``, a random value
      will be chosen for each insertion.

    EXAMPLES::

        sage: from sage.coding.punctured_code import _insert_punctured_positions
        sage: _insert_punctured_positions([1,2,3,4], {2,4,5}, 1)
        [1, 2, 1, 3, 1, 1, 4]
    """
    F = l[0].base_ring()
    final = [None] * (len(l) + len(punctured_points))
    for i in punctured_points:
        if value is None:
            final[i] = F.random_element()
        else:
            final[i] = value
    index = 0
    for i in range(len(final)):
        if final[i] is None:
            final[i] = l[index]
            index += 1
    return final


class PuncturedCode(AbstractLinearCode):
    r"""
    Representation of a punctured code.

    - ``C`` -- a linear code

    - ``positions`` -- the positions where ``C`` will be punctured. It can be
      either an integer if one need to puncture only one position, a list or a
      set of positions to puncture. If the same position is passed several
      times, it will be considered only once.

    EXAMPLES::

        sage: C = codes.random_linear_code(GF(7), 11, 5)
        sage: Cp = codes.PuncturedCode(C, 3)
        sage: Cp
        Puncturing of [11, 5] linear code over GF(7) on position(s) [3]

        sage: Cp = codes.PuncturedCode(C, {3, 5})
        sage: Cp
        Puncturing of [11, 5] linear code over GF(7) on position(s) [3, 5]
    """
    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, C, positions):
        r"""
        TESTS:

        If one of the positions to puncture is bigger than the length of ``C``,
        an exception will be raised::

            sage: C = codes.random_linear_code(GF(7), 11, 5)
            sage: Cp = codes.PuncturedCode(C, {4,8,15})
            Traceback (most recent call last):
            ...
            ValueError: Positions to puncture must be positive integers smaller
            than the length of the provided code
        """
        if not isinstance(positions, (Integer, int, set, list)):
            raise TypeError("positions must be either a Sage Integer, a Python int, a set or a list")
        if isinstance(positions, (list, set)) and not all(isinstance(i, (int, Integer)) for i in positions):
            raise TypeError("if positions is a list or a set, it has to contain only Python ints or Sage Integers")
        if isinstance(positions, (Integer, int)):
            positions = {positions}
        if isinstance(positions, list):
            positions = set(positions)
        if not isinstance(C, AbstractLinearCode):
            raise ValueError("Provided code must be a linear code")
        if not all(i in range(C.length()) for i in positions):
            raise ValueError("Positions to puncture must be positive integers smaller than the length of the provided code")
        super().__init__(C.base_ring(), C.length() - len(positions),
                         "PuncturedMatrix", "OriginalCode")
        self._original_code = C
        self._positions = positions

    def __eq__(self, other):
        r"""
        Test equality between two Punctured codes.

        EXAMPLES::

            sage: C = codes.random_linear_code(GF(7), 11, 5)
            sage: Cp1 = codes.PuncturedCode(C, 2)
            sage: Cp2 = codes.PuncturedCode(C, 2)
            sage: Cp1 == Cp2
            True
        """
        return isinstance(other, PuncturedCode) \
                and self.punctured_positions() == other.punctured_positions() \
                and self.original_code() == other.original_code()

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.random_linear_code(GF(7), 11, 5)
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: Cp
            Puncturing of [11, 5] linear code over GF(7) on position(s) [3]
        """
        return "Puncturing of %s on position(s) %s"\
                % (self.original_code(), list(self.punctured_positions()))

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: C = codes.random_linear_code(GF(7), 11, 5)
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: latex(Cp)
            \textnormal{Puncturing of [11, 5] linear code over GF(7) on position(s) } [3]
        """
        return "\\textnormal{Puncturing of %s on position(s) } %s"\
                % (self.original_code(), list(self.punctured_positions()))

    def punctured_positions(self):
        r"""
        Return the list of positions which were punctured on the original code.

        EXAMPLES::

            sage: C = codes.random_linear_code(GF(7), 11, 5)
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: Cp.punctured_positions()
            {3}
        """
        return self._positions

    def original_code(self):
        r"""
        Return the linear code which was punctured to get ``self``.

        EXAMPLES::

            sage: C = codes.random_linear_code(GF(7), 11, 5)
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: Cp.original_code()
            [11, 5] linear code over GF(7)
        """
        return self._original_code

    def dimension(self):
        r"""
        Return the dimension of ``self``.

        EXAMPLES::

            sage: set_random_seed(42)
            sage: C = codes.random_linear_code(GF(7), 11, 5)
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: Cp.dimension()
            5
        """
        if hasattr(self, '_dimension'):
            return self._dimension
        self._dimension = self.generator_matrix().rank()
        return self._dimension

    def random_element(self, *args, **kwds):
        r"""
        Return a random codeword of ``self``.

        This method does not trigger the computation of
        ``self``'s :meth:`sage.coding.linear_code_no_metric.generator_matrix`.

        INPUT:

        - ``agrs``, ``kwds`` -- extra positional arguments passed to
          :meth:`sage.modules.free_module.random_element`

        EXAMPLES::

            sage: C = codes.random_linear_code(GF(7), 11, 5)
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: Cp.random_element() in Cp
            True
        """
        C_original = self.original_code()
        m = (C_original.base_ring() ** C_original.dimension()).random_element()
        c = C_original.encode(m)
        return _puncture(c, self.punctured_positions())

    def encode(self, m, original_encode=False, encoder_name=None, **kwargs):
        r"""
        Transform an element of the message space into an element of the code.

        INPUT:

        - ``m`` -- a vector of the message space of the code

        - ``original_encode`` -- boolean (default: ``False``); if this is set
          to ``True``, ``m`` will be encoded using an Encoder of ``self``'s
          :meth:`original_code`. This allow to avoid the computation of a
          generator matrix for ``self``.

        - ``encoder_name`` -- (default: ``None``) name of the encoder which will be used
          to encode ``word``. The default encoder of ``self`` will be used if
          default value is kept.

        OUTPUT: an element of ``self``

        EXAMPLES::

           sage: M = matrix(GF(7), [[1, 0, 0, 0, 3, 4, 6],
           ....:                    [0, 1, 0, 6, 1, 6, 4],
           ....:                    [0, 0, 1, 5, 2, 2, 4]])
           sage: C_original = LinearCode(M)
           sage: Cp = codes.PuncturedCode(C_original, 2)
           sage: m = vector(GF(7), [1, 3, 5])
           sage: Cp.encode(m)
            (1, 3, 5, 5, 0, 2)
        """
        if original_encode:
            c = self.original_code().encode(m, encoder_name, **kwargs)
            return _puncture(c, self.punctured_positions, self)
        return self.encoder(encoder_name, **kwargs).encode(m)

    @cached_method
    def structured_representation(self):
        r"""
        Return ``self`` as a structured code object.

        If ``self`` has a specific structured representation (e.g. a punctured GRS code is
        a GRS code too), it will return this representation, else it returns a
        :class:`sage.coding.linear_code.LinearCode`.

        EXAMPLES:

        We consider a GRS code::

            sage: C_grs = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)

        A punctured GRS code is still a GRS code::

            sage: Cp_grs = codes.PuncturedCode(C_grs, 3)
            sage: Cp_grs.structured_representation()
            [39, 12, 28] Reed-Solomon Code over GF(59)

        Another example with structureless linear codes::

            sage: set_random_seed(42)
            sage: C_lin  = codes.random_linear_code(GF(2), 10, 5)
            sage: Cp_lin = codes.PuncturedCode(C_lin, 2)
            sage: Cp_lin.structured_representation()
            [9, 5] linear code over GF(2)
        """
        C = self.original_code()
        pts = copy(self.punctured_positions())
        list_pts = list(pts)
        while isinstance(C, PuncturedCode):
            cur_pts = list(C.punctured_positions())
            list_len = len(list_pts)
            for p in cur_pts:
                for i in range(list_len):
                    if (p <= list_pts[i]):
                        list_pts[i] += 1
            list_pts += cur_pts
            C = C.original_code()
        return C._punctured_form(set(list_pts))


class PuncturedCodePuncturedMatrixEncoder(Encoder):
    r"""
    Encoder using original code generator matrix to compute the punctured code's one.

    INPUT:

    - ``code`` -- the associated code of this encoder

    EXAMPLES::

        sage: C = codes.random_linear_code(GF(7), 11, 5)
        sage: Cp = codes.PuncturedCode(C, 3)
        sage: E = codes.encoders.PuncturedCodePuncturedMatrixEncoder(Cp)
        sage: E
        Punctured matrix-based encoder for the
         Puncturing of [11, 5] linear code over GF(7) on position(s) [3]
    """

    def __init__(self, code):
        r"""
        TESTS:

        If ``code`` is not a ``PuncturedCode``, an exception is raised::

            sage: C = codes.random_linear_code(GF(7), 11, 5)
            sage: codes.encoders.PuncturedCodePuncturedMatrixEncoder(C)
            Traceback (most recent call last):
            ...
            TypeError: code has to be an instance of PuncturedCode class
        """
        if not isinstance(code, PuncturedCode):
            raise TypeError("code has to be an instance of PuncturedCode class")
        super().__init__(code)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.random_linear_code(GF(7), 11, 5)
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: E = codes.encoders.PuncturedCodePuncturedMatrixEncoder(Cp)
            sage: E
            Punctured matrix-based encoder for the Puncturing of [11, 5] linear code over GF(7) on position(s) [3]
        """
        return "Punctured matrix-based encoder for the %s" % self.code()

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: C = codes.random_linear_code(GF(7), 11, 5)
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: E = codes.encoders.PuncturedCodePuncturedMatrixEncoder(Cp)
            sage: latex(E)
            \textnormal{Punctured matrix-based encoder for the }\textnormal{Puncturing of [11, 5] linear code over GF(7) on position(s) } [3]
        """
        return "\\textnormal{Punctured matrix-based encoder for the }%s" % self.code()._latex_()

    @cached_method
    def generator_matrix(self):
        r"""
        Return a generator matrix of the associated code of ``self``.

        EXAMPLES::

            sage: set_random_seed(10)
            sage: C = codes.random_linear_code(GF(7), 11, 5)
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: E = codes.encoders.PuncturedCodePuncturedMatrixEncoder(Cp)
            sage: E.generator_matrix()
            [1 0 0 0 0 5 2 6 0 6]
            [0 1 0 0 0 5 2 2 1 1]
            [0 0 1 0 0 6 2 4 0 4]
            [0 0 0 1 0 0 6 3 3 3]
            [0 0 0 0 1 0 1 3 4 3]
        """
        C = self.code().original_code()
        pos = self.code().punctured_positions()
        M = C.generator_matrix()
        G = M.delete_columns(list(pos))
        G = G.echelon_form()
        k = G.rank()
        M = G[:k]
        M.set_immutable()
        return M


class PuncturedCodeOriginalCodeDecoder(Decoder):
    r"""
    Decoder decoding through a decoder over the original code of its punctured code.

    INPUT:

    - ``code`` -- the associated code of this encoder

    - ``strategy`` -- (default: ``None``) the strategy used to decode.
      The available strategies are:

      * ``'error-erasure'`` -- uses an error-erasure decoder over the original
        code if available, fails otherwise.

      * ``'random-values'`` -- fills the punctured positions with random elements
        in ``code``'s base field and tries to decode using
        the default decoder of the original code

      * ``'try-all'`` -- fills the punctured positions with every possible
        combination of symbols until decoding succeeds, or until every
        combination have been tried

      * ``None`` -- uses ``error-erasure`` if an error-erasure decoder is
        available, switch to ``random-values`` behaviour otherwise

    - ``original_decoder`` -- (default: ``None``) the decoder that will be used over the original code.
      It has to be a decoder object over the original code.
      This argument takes precedence over ``strategy``: if both ``original_decoder`` and ``strategy``
      are filled, ``self`` will use the ``original_decoder`` to decode over the original code.
      If ``original_decoder`` is set to ``None``, it will use the decoder picked by ``strategy``.

    - ``**kwargs`` -- all extra arguments are forwarded to original code's decoder

    EXAMPLES::

        sage: C = codes.GeneralizedReedSolomonCode(GF(16, 'a').list()[:15], 7)
        sage: Cp = codes.PuncturedCode(C, 3)
        sage: codes.decoders.PuncturedCodeOriginalCodeDecoder(Cp)
        Decoder of Puncturing of [15, 7, 9] Reed-Solomon Code over GF(16) on position(s) [3]
         through Error-Erasure decoder for [15, 7, 9] Reed-Solomon Code over GF(16)

    As seen above, if all optional are left blank, and if an error-erasure
    decoder is available, it will be chosen as the original decoder.  Now, if
    one forces ``strategy`` to ``'try-all'`` or ``'random-values'``, the default
    decoder of the original code will be chosen, even if an error-erasure is
    available::

        sage: C = codes.GeneralizedReedSolomonCode(GF(16, 'a').list()[:15], 7)
        sage: Cp = codes.PuncturedCode(C, 3)
        sage: D = codes.decoders.PuncturedCodeOriginalCodeDecoder(Cp, strategy='try-all')
        sage: "error-erasure" in D.decoder_type()
        False

    And if one fills ``original_decoder`` and ``strategy`` fields with
    contradictory elements, the ``original_decoder`` takes precedence::

        sage: C = codes.GeneralizedReedSolomonCode(GF(16, 'a').list()[:15], 7)
        sage: Cp = codes.PuncturedCode(C, 3)
        sage: Dor = C.decoder("Gao")
        sage: D = codes.decoders.PuncturedCodeOriginalCodeDecoder(Cp, original_decoder=Dor,
        ....:                                                     strategy='error-erasure')
        sage: D.original_decoder() == Dor
        True
    """

    def __init__(self, code, strategy=None, original_decoder=None, **kwargs):
        r"""
        TESTS:

        If ``code`` is not a ``PuncturedCode``, an exception is raised::

            sage: C = codes.random_linear_code(GF(7), 11, 5)
            sage: codes.decoders.PuncturedCodeOriginalCodeDecoder(C)
            Traceback (most recent call last):
            ...
            TypeError: code has to be an instance of PuncturedCode class

        If one tries to pass an original_decoder whose associated code is not the original
        code of ``self``, it returns an error::

            sage: C = codes.GeneralizedReedSolomonCode(GF(16, 'a').list()[:15], 7)
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: C2 = codes.GeneralizedReedSolomonCode(GF(7).list()[:6], 3)
            sage: D = Cp.decoder(original_decoder = C2.decoder())
            Traceback (most recent call last):
            ...
            ValueError: Original decoder must have the original code of its associated punctured code as associated code

        If one tries to use ``'error-erasure'`` strategy when the original code has no such
        decoder, it returns an error::

            sage: C = codes.random_linear_code(GF(7), 11, 5)
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: D = codes.decoders.PuncturedCodeOriginalCodeDecoder(Cp, strategy = 'error-erasure')
            Traceback (most recent call last):
            ...
            ValueError: Original code has no error-erasure decoder
        """
        if not isinstance(code, PuncturedCode):
            raise TypeError("code has to be an instance of PuncturedCode class")

        original_code = code.original_code()
        if original_decoder is not None:
            if not isinstance(original_decoder, Decoder):
                raise TypeError("original_decoder must be a decoder object")
            if not original_decoder.code() == original_code:
                raise ValueError("Original decoder must have the original code of its associated punctured code as associated code")
            if 'error-erasure' in original_decoder.decoder_type():
                strategy = 'error-erasure'
            self._original_decoder = original_decoder
        elif strategy == 'error-erasure':
            error_erasure = False
            for D in original_code._registered_decoders.values():
                if 'error-erasure' in D._decoder_type:
                    error_erasure = True
                    self._original_decoder = D(original_code, **kwargs)
                    break
            if not error_erasure:
                raise ValueError("Original code has no error-erasure decoder")
        elif strategy == 'random-values' or strategy == 'try-all':
            self._original_decoder = code.original_code().decoder(**kwargs)
        else:
            error_erasure = False
            for D in original_code._registered_decoders.values():
                if 'error-erasure' in D._decoder_type:
                    error_erasure = True
                    self._original_decoder = D(original_code, **kwargs)
                    break
            if not error_erasure:
                self._original_decoder = original_code.decoder(**kwargs)
        self._strategy = strategy
        self._decoder_type = copy(self._decoder_type)
        self._decoder_type.remove("dynamic")
        self._decoder_type = self._original_decoder.decoder_type()
        super().__init__(code, code.ambient_space(),
                         self._original_decoder.connected_encoder())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(16, 'a').list()[:15], 7)
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: D = codes.decoders.PuncturedCodeOriginalCodeDecoder(Cp)
            sage: D
            Decoder of Puncturing of [15, 7, 9] Reed-Solomon Code over GF(16) on position(s) [3] through Error-Erasure decoder for [15, 7, 9] Reed-Solomon Code over GF(16)
        """
        return "Decoder of %s through %s" % (self.code(), self.original_decoder())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(16, 'a').list()[:15], 7)
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: D = codes.decoders.PuncturedCodeOriginalCodeDecoder(Cp)
            sage: latex(D)
            \textnormal{Decoder of } Puncturing of [15, 7, 9] Reed-Solomon Code over GF(16) on position(s) [3] \textnormal{ through } Error-Erasure decoder for [15, 7, 9] Reed-Solomon Code over GF(16)
        """
        return "\\textnormal{Decoder of } %s \\textnormal{ through } %s" % (self.code(), self.original_decoder())

    def original_decoder(self):
        r"""
        Return the decoder over the original code that will be used to decode words of
        :meth:`sage.coding.decoder.Decoder.code`.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(16, 'a').list()[:15], 7)
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: D = codes.decoders.PuncturedCodeOriginalCodeDecoder(Cp)
            sage: D.original_decoder()
            Error-Erasure decoder for [15, 7, 9] Reed-Solomon Code over GF(16)
        """
        return self._original_decoder

    def decode_to_code(self, y):
        r"""
        Decode ``y`` to an element in :meth:`sage.coding.decoder.Decoder.code`.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(16, 'a').list()[:15], 7)
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: D = codes.decoders.PuncturedCodeOriginalCodeDecoder(Cp)
            sage: c = Cp.random_element()
            sage: Chan = channels.StaticErrorRateChannel(Cp.ambient_space(), 3)
            sage: y = Chan(c)
            sage: y in Cp
            False
            sage: D.decode_to_code(y) == c
            True
        """
        D = self.original_decoder()
        C = self.code()
        A = C.original_code().ambient_space()
        Cor = C.original_code()
        pts = C.punctured_positions()
        F = self.code().base_field()
        zero, one = F.zero(), F.one()
        if "error-erasure" in D.decoder_type():
            if isinstance(y, (tuple, list)):
                y, e = y[0], y[1]
                e_list = e.list()
                e_list = _insert_punctured_positions(e_list, pts, one)
            else:
                e_list = [one if i in pts else zero for i in range(Cor.length())]
            e = vector(GF(2), e_list)
            yl = y.list()
            yl = _insert_punctured_positions(yl, pts, zero)
            y = A(yl)
            return _puncture(D.decode_to_code((y, e)), pts)
        elif self._strategy == 'try-all':
            end = False
            yl = y.list()
            I = iter(VectorSpace(F, len(pts)))
            list_pts = list(pts)
            list_pts.sort()
            shift = 0
            for i in list_pts:
                yl.insert(i + shift, zero)
                shift += 1
            values = next(I)
            while not end:
                try:
                    shift = 0
                    for i in list_pts:
                        yl[i + shift] = values[shift]
                        shift += 1
                    y = A(yl)
                    values = next(I)
                    try:
                        c_or = self.original_decoder().decode_to_code(y)
                        end = True
                        break
                    except Exception:
                        pass
                except StopIteration:
                    raise DecodingError
            return _puncture(c_or, pts)
        A = Cor.ambient_space()
        yl = y.list()
        yl = _insert_punctured_positions(yl, pts)
        y = A(yl)
        return _puncture(D.decode_to_code(y), pts)

    def decoding_radius(self, number_erasures=None):
        r"""
        Return the maximal number of errors that ``self`` can decode.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(16, 'a').list()[:15], 7)
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: D = codes.decoders.PuncturedCodeOriginalCodeDecoder(Cp)
            sage: D.decoding_radius(2)
            2
        """
        punctured = len(self.code().punctured_positions())
        D = self.original_decoder()
        if self._strategy != 'try-all' and "error-erasure" not in D.decoder_type():
            if D.decoding_radius() - punctured >= 0:
                return D.decoding_radius() - punctured
            else:
                return 0
        elif "error-erasure" in D.decoder_type() and number_erasures is not None:
            diff = self.code().original_code().minimum_distance() - number_erasures - punctured - 1
            if diff <= 0:
                raise ValueError("The number of erasures exceeds decoding capability")
            return diff // 2
        elif "error-erasure" in D.decoder_type() and number_erasures is None:
            raise ValueError("You must provide the number of erasures")


####################### registration ###############################

PuncturedCode._registered_encoders["PuncturedMatrix"] = PuncturedCodePuncturedMatrixEncoder
PuncturedCode._registered_decoders["OriginalCode"] = PuncturedCodeOriginalCodeDecoder
PuncturedCodeOriginalCodeDecoder._decoder_type = {"dynamic"}
