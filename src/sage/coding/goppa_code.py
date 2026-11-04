# sage.doctest: needs sage.modules sage.rings.finite_rings
r"""
Goppa code

This module implements Goppa codes, an encoder for them, and a decoder for
binary Goppa codes using Patterson's algorithm [Pat1975]_.

EXAMPLES::

    sage: F = GF(2^6)
    sage: R.<x> = F[]
    sage: g = x^9 + 1
    sage: L = [a for a in F.list() if g(a) != 0]
    sage: C = codes.GoppaCode(g, L)
    sage: C
    [55, 16] Goppa code over GF(2)
    sage: E = codes.encoders.GoppaCodeEncoder(C)
    sage: E
    Encoder for [55, 16] Goppa code over GF(2)

A binary Goppa code with an irreducible generating polynomial can be decoded
with Patterson's algorithm::

    sage: F = GF(2^4)
    sage: R.<x> = F[]
    sage: g = R.irreducible_element(3)
    sage: L = [a for a in F.list() if g(a) != 0]
    sage: C = codes.GoppaCode(g, L)
    sage: D = C.decoder('Patterson')
    sage: c = C.random_element()
    sage: y = c + C.ambient_space().gen(0) + C.ambient_space().gen(5)
    sage: D.decode_to_code(y) == c
    True

AUTHORS:

- Filip Ion, Marketa Slukova (2019-06): initial version

- Nolan Colter (2026-10): Patterson decoder
"""
#*****************************************************************************
#       Copyright (C) 2019 Filip Ion <filip.alexandru.ion@gmail.com>,
#                          Marketa Slukova <em.slukova@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import annotations

from copy import copy
from typing import TYPE_CHECKING

from sage.coding.linear_code import AbstractLinearCode
from sage.coding.encoder import Encoder
from sage.coding.decoder import Decoder, DecodingError
from sage.modules.free_module_element import vector
import sage.coding.codes_catalog as codes

if TYPE_CHECKING:
    from sage.modules.free_module_element import FreeModuleElement
    from sage.rings.polynomial.polynomial_element import Polynomial


def _columnize(element):
    """
    Convert a finite field element to a column vector over the prime field.

    TESTS::

        sage: from sage.coding.goppa_code import _columnize
        sage: F.<a> = GF(2^6)
        sage: _columnize(a)
        [0]
        [1]
        [0]
        [0]
        [0]
        [0]
        sage: _columnize(a+1)
        [1]
        [1]
        [0]
        [0]
        [0]
        [0]
    """
    v = vector(element)
    return v.column()


def _split_even_odd(f: Polynomial) -> tuple[Polynomial, Polynomial]:
    r"""
    Return the polynomials `f_0` and `f_1` with `f = f_0^2 + x f_1^2`.

    The coefficients of ``f`` lie in a finite field of characteristic 2,
    where every element has exactly one square root. `f_0` is built from the
    square roots of the coefficients of the even powers of `x`, and `f_1`
    from those of the odd powers.

    INPUT:

    - ``f`` -- a polynomial over a finite field of characteristic 2

    OUTPUT: a pair of polynomials

    TESTS::

        sage: from sage.coding.goppa_code import _split_even_odd
        sage: F.<a> = GF(2^3)
        sage: R.<x> = F[]
        sage: f = x^5 + a*x^2 + (a + 1)*x + a^2
        sage: f0, f1 = _split_even_odd(f)
        sage: f0, f1
        ((a^2 + a)*x + a, x^2 + a^2 + a + 1)
        sage: f0^2 + x * f1^2 == f
        True
        sage: _split_even_odd(R.zero())
        (0, 0)
    """
    ring = f.parent()
    coefficients = f.list()

    # Square root of each coefficient
    # Put even powers of x into f0, odd powers into f1
    f0 = ring([c.sqrt() for c in coefficients[0::2]])
    f1 = ring([c.sqrt() for c in coefficients[1::2]])
    return f0, f1


class GoppaCode(AbstractLinearCode):
    r"""
    Implementation of Goppa codes.

    Goppa codes are a generalization of narrow-sense BCH codes.
    These codes are defined by a generating polynomial `g` over a finite field
    `\GF{p^m}`, and a defining set `L` of elements from `\GF{p^m}`, which are not roots
    of `g`. The number of defining elements determines the length of the code.

    In binary cases, the minimum distance is `2t + 1`, where `t` is the degree
    of `g`.

    INPUT:

    - ``generating_pol`` -- a monic polynomial with coefficients in a finite
      field `\GF{p^m}`, the code is defined over `\GF{p}`, `p` must be a prime number

    - ``defining_set`` -- set of elements of `\GF{p^m}` that are not roots
      of `g`, its cardinality is the length of the code

    EXAMPLES::

        sage: F = GF(2^6)
        sage: R.<x> = F[]
        sage: g = x^9 + 1
        sage: L = [a for a in F.list() if g(a) != 0]
        sage: C = codes.GoppaCode(g, L)
        sage: C
        [55, 16] Goppa code over GF(2)
    """
    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, generating_pol, defining_set):
        """
        Initialize.

        TESTS::

            sage: F = GF(2^6)
            sage: R.<x> = F[]
            sage: g = x^9 + 1
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: TestSuite(C).run()
        """
        self._field = generating_pol.base_ring().prime_subfield()
        self._length = len(defining_set)
        self._generating_pol = generating_pol
        self._defining_set = defining_set

        super().__init__(self._field, self._length, "GoppaEncoder", "Syndrome")

        if not generating_pol.is_monic():
            raise ValueError("generating polynomial must be monic")
        F = self._field
        if (not F.is_field() or not F.is_finite()):
            raise ValueError("generating polynomial must be defined over a finite field")
        for a in defining_set:
            if generating_pol(a) == 0:
                raise ValueError("defining elements cannot be roots of generating polynomial")

    def _repr_(self):
        """
        Representation of a Goppa code.

        EXAMPLES::

            sage: F  = GF(2^3)
            sage: R.<x> = F[]
            sage: g = x^2 + x+ 1
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: C
            [8, 2] Goppa code over GF(2)
        """
        return "[{}, {}] Goppa code over GF({})".format(
            self.length(), self.dimension(),
            self.base_field().cardinality())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: F  = GF(2^3)
            sage: R.<x> = F[]
            sage: g = x^2 + x+ 1
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: latex(C)
            [8, 2]\text{ Goppa code over }\Bold{F}_{2}
        """
        return r"[{}, {}]\text{{ Goppa code over }}{}".format(self.length(), self.dimension(),
                                                              self.base_field()._latex_())

    def __eq__(self, other):
        """
        Test equality with ``other``.

        Two Goppa codes are considered the same when defining sets and
        generating polynomials are the same.

        EXAMPLES::

            sage: F  = GF(2^3)
            sage: R.<x> = F[]
            sage: g = x^2 + x+ 1
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: D = codes.GoppaCode(g, L)
            sage: C == D
            True

        Note that equality check will be false if ``other`` represents the same
        linear code as ``self`` but not constructed as a Goppa code::

            sage: E = LinearCode(C.generator_matrix())
            sage: C == E
            False
        """
        return (isinstance(other, GoppaCode)
           and self.length() == other.length()
           and self._generating_pol == other._generating_pol
           and self._defining_set == other._defining_set)

    def parity_check_matrix(self):
        r"""
        Return a parity check matrix for ``self``.

        The element in row `t`, column `i` is `h[i](D[i]^t)`, where:

        - `h[i]` -- is the inverse of `g(D[i])`
        - `D[i]` -- is the `i`-th element of the defining set

        In the resulting `d \times n` matrix we interpret each entry as an
        `m`-column vector and return a `dm \times n` matrix.

        EXAMPLES::

            sage: F  = GF(2^3)
            sage: R.<x> = F[]
            sage: g = x^2 + x+ 1
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: C
            [8, 2] Goppa code over GF(2)
            sage: C.parity_check_matrix()
            [1 0 0 0 0 0 0 1]
            [0 0 1 0 1 1 1 0]
            [0 1 1 1 0 0 1 0]
            [0 1 1 1 1 1 1 1]
            [0 1 0 1 1 0 1 0]
            [0 0 1 1 1 1 0 0]
        """
        g = self._generating_pol
        F = g.base_ring()
        n = self._length
        d = g.degree()
        alpha = F.primitive_element()

        D = self._defining_set
        h = [(g(D[i]).inverse_of_unit()) for i in range(n)]

        #assemble top row
        M = _columnize(alpha)
        for i in range(n):
            v = _columnize(h[i])
            M = M.augment(v)
        M = M.delete_columns([0])
        old = M

        for t in range(1,d):
            #assemble row
            M = _columnize(alpha)
            for i in range(n):
                v = _columnize(h[i]*(D[i]**t))
                M = M.augment(v)
            M = M.delete_columns([0])
            new = M
            old = old.stack(new)

        return old

    def _parity_check_matrix_vandermonde(self):
        """
        Return a parity check matrix for ``self`` using Vandermonde matrix.

        EXAMPLES::

            sage: F  = GF(2^3)
            sage: R.<x> = F[]
            sage: g = x^2 + x+ 1
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: C
            [8, 2] Goppa code over GF(2)
            sage: C._parity_check_matrix_vandermonde()
            [1 0 0 0 0 0 0 1]
            [0 0 1 0 1 1 1 0]
            [0 1 1 1 0 0 1 0]
            [---------------]
            [0 1 1 1 1 1 1 1]
            [0 1 0 1 1 0 1 0]
            [0 0 1 1 1 1 0 0]
        """
        L = self._defining_set
        g = self._generating_pol
        t = g.degree()

        from sage.matrix.constructor import matrix, diagonal_matrix, block_matrix

        V = matrix.vandermonde(L)
        V = V.transpose()
        GL = [g(i) for i in L]
        GLI = [j.inverse_of_unit() for j in GL]
        D = diagonal_matrix(GLI)
        VF = matrix([V.row(i) for i in range(t)])
        H = VF*D

        matrices = [matrix([vector(i) for i in H.row(j)]) for j in range(t)]
        matrices = [m.transpose() for m in matrices]

        m = block_matrix(t, 1, matrices)
        return m

    def distance_bound(self):
        """
        Return a lower bound for the minimum distance of the code.

        Computed using the degree of the generating polynomial of ``self``.
        The minimum distance is guaranteed to be bigger than or equal to this bound.

        EXAMPLES::

            sage: F = GF(2^3)
            sage: R.<x> = F[]
            sage: g = x^2 + x + 1
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: C
            [8, 2] Goppa code over GF(2)
            sage: C.distance_bound()
            3
            sage: C.minimum_distance()                                                  # needs sage.libs.gap
            5
        """
        return 1 + (self._generating_pol).degree()


class GoppaCodeEncoder(Encoder):
    r"""
    Encoder for Goppa codes.

    Encodes words represented as vectors of length `k`, where `k` is
    the dimension of ``self``, with entries from `\GF{p}`, the prime field of
    the base field of the generating polynomial of ``self``, into codewords
    of length `n`, with entries from `\GF{p}`.

    EXAMPLES::

        sage: F = GF(2^3)
        sage: R.<x> = F[]
        sage: g = x^2 + x + 1
        sage: L = [a for a in F.list() if g(a) != 0]
        sage: C = codes.GoppaCode(g, L)
        sage: C
        [8, 2] Goppa code over GF(2)
        sage: E = codes.encoders.GoppaCodeEncoder(C)
        sage: E
        Encoder for [8, 2] Goppa code over GF(2)
        sage: word = vector(GF(2), (0, 1))
        sage: c = E.encode(word)
        sage: c
        (0, 1, 1, 1, 1, 1, 1, 0)
        sage: c in C
        True
    """
    def __init__(self, code):
        """
        Initialize.

        TESTS::

            sage: F = GF(2^3)
            sage: R.<x> = F[]
            sage: g = x^2 + x + 1
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: E = codes.encoders.GoppaCodeEncoder(C)
            sage: TestSuite(E).run()
        """
        super().__init__(code)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: F = GF(2^3)
            sage: R.<x> = F[]
            sage: g = (x^2 + x + 1)^2
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: E = codes.encoders.GoppaCodeEncoder(C)
            sage: E
            Encoder for [8, 2] Goppa code over GF(2)
        """
        return "Encoder for {}".format(self.code())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: F = GF(2^3)
            sage: R.<x> = F[]
            sage: g = (x^2 + x + 1)^2
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: E = codes.encoders.GoppaCodeEncoder(C)
            sage: latex(E)
            \text{Encoder for }[8, 2]\text{ Goppa code over }\Bold{F}_{2}
        """
        return r"\text{{Encoder for }}{}".format(self.code()._latex_())

    def __eq__(self, other):
        """
        Test equality with ``other``.

        EXAMPLES::

            sage: F = GF(2^3)
            sage: R.<x> = F[]
            sage: g = (x^2 + x + 1)^2
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: E1 = codes.encoders.GoppaCodeEncoder(C)
            sage: E2 = codes.encoders.GoppaCodeEncoder(C)
            sage: E1 == E2
            True
        """
        return (isinstance(other, GoppaCodeEncoder)
           and self.code() == other.code())

    def generator_matrix(self):
        r"""
        A generator matrix for ``self``.

        Dimension of resulting matrix is `k \times n`, where `k` is
        the dimension of ``self`` and `n` is the length of ``self``.

        EXAMPLES::

            sage: F = GF(2^3)
            sage: R.<x> = F[]
            sage: g = (x^2 + x + 1)^2
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: C
            [8, 2] Goppa code over GF(2)
            sage: C.generator_matrix()
            [1 0 0 1 0 1 1 1]
            [0 1 1 1 1 1 1 0]
        """
        c = self.code()
        pmat = c.parity_check_matrix()
        aux = codes.from_parity_check_matrix(pmat)
        return aux.generator_matrix()


class GoppaCodePattersonDecoder(Decoder):
    r"""
    Decoder for binary Goppa codes using Patterson's algorithm [Pat1975]_.

    It corrects up to `t` errors, where `t` is the degree of the generating
    polynomial `g`. The code must be binary, `g` must be squarefree, and every
    irreducible factor of `g` must have degree above `(t - 1)/2`. That means
    `g` is irreducible, or `t` is even and `g` is a product of two different
    irreducible polynomials of degree `t/2`. Other Goppa codes are refused,
    since Patterson's algorithm can fail to decode some of their words with at
    most `t` errors.

    Square roots modulo `g` are taken as in [Hub1996]_.

    INPUT:

    - ``code`` -- a Goppa code over `\GF{2}` with a squarefree generating
      polynomial, every irreducible factor of degree above `(t - 1)/2`

    - ``check`` -- boolean (default: ``True``); if ``False``, don't refuse a
      squarefree generating polynomial because of its factor degrees. Some
      words with at most `t` errors can then raise
      :exc:`~sage.coding.decoder.DecodingError`, whose message still says the
      word has more than `t` errors. So ``decoding_radius`` and the decoder type
      ``"always-succeed"`` are not guarantees for such a code.

    EXAMPLES::

        sage: F = GF(2^3)
        sage: R.<x> = F[]
        sage: g = x^2 + x + 1
        sage: L = [a for a in F.list() if g(a) != 0]
        sage: C = codes.GoppaCode(g, L)
        sage: D = codes.decoders.GoppaCodePattersonDecoder(C)
        sage: D
        Patterson decoder for [8, 2] Goppa code over GF(2)
        sage: D.decoding_radius()
        2

    To build a code the decoder accepts, choose `g` irreducible::

        sage: F = GF(2^4)
        sage: R.<x> = F[]
        sage: g = R.irreducible_element(3)
        sage: L = [a for a in F.list() if g(a) != 0]
        sage: C = codes.GoppaCode(g, L)
        sage: D = codes.decoders.GoppaCodePattersonDecoder(C)
        sage: D.decoding_radius()
        3

    With ``check=False``, the decoder accepts a `g` with small factors. Most
    words still decode, and this one with 5 errors doesn't::

        sage: F = GF(2^5)
        sage: R.<x> = F[]
        sage: g = (x^2 + x + 1) * (x^3 + x + 1)
        sage: L = [a for a in F.list() if g(a) != 0]
        sage: C = codes.GoppaCode(g, L)
        sage: D = C.decoder("Patterson", check=False)
        sage: def word_with_errors(positions):
        ....:     return vector(GF(2), [int(i in positions) for i in range(C.length())])
        sage: D.decode_to_code(word_with_errors([0, 1, 2, 3, 4])).is_zero()
        True
        sage: D.decode_to_code(word_with_errors([0, 1, 6, 26, 29]))
        Traceback (most recent call last):
        ...
        DecodingError: Patterson decoding failed, the word has more than 5 errors

    A generating polynomial that isn't squarefree is still refused::

        sage: F = GF(2^3)
        sage: R.<x> = F[]
        sage: g = (x^2 + x + 1)^2
        sage: L = [a for a in F.list() if g(a) != 0]
        sage: C = codes.GoppaCode(g, L)
        sage: codes.decoders.GoppaCodePattersonDecoder(C, check=False)
        Traceback (most recent call last):
        ...
        ValueError: Patterson decoding needs a squarefree generating polynomial
    """
    def __init__(self, code: GoppaCode, check: bool = True) -> None:
        r"""
        Initialize.

        TESTS::

            sage: F = GF(2^3)
            sage: R.<x> = F[]
            sage: g = x^2 + x + 1
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: D = codes.decoders.GoppaCodePattersonDecoder(C)
            sage: TestSuite(D).run()

        If ``code`` is not a Goppa code, an error is raised::

            sage: C = codes.random_linear_code(GF(2), 10, 4)
            sage: codes.decoders.GoppaCodePattersonDecoder(C)
            Traceback (most recent call last):
            ...
            ValueError: Patterson decoding needs a Goppa code

        An error is also raised if ``code`` is not binary::

            sage: F = GF(3^2)
            sage: R.<x> = F[]
            sage: g = x^2 + 1
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: C
            [7, 3] Goppa code over GF(3)
            sage: codes.decoders.GoppaCodePattersonDecoder(C)
            Traceback (most recent call last):
            ...
            ValueError: Patterson decoding needs a binary Goppa code, not one over Finite Field of size 3

        A squarefree polynomial with large enough factors is accepted::

            sage: F = GF(2^5)
            sage: R.<x> = F[]
            sage: g = (x^3 + x + 1) * (x^3 + x^2 + 1)
            sage: C = codes.GoppaCode(g, F.list())
            sage: D = codes.decoders.GoppaCodePattersonDecoder(C)
            sage: D.decode_to_code(vector(GF(2), [1]*6 + [0]*26)).is_zero()
            True

        An error is raised if the generating polynomial is not squarefree, even
        when it has no roots in the defining set::

            sage: F = GF(2^3)
            sage: R.<x> = F[]
            sage: g = (x^2 + x + 1)^2
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: codes.decoders.GoppaCodePattersonDecoder(C)
            Traceback (most recent call last):
            ...
            ValueError: Patterson decoding needs a squarefree generating polynomial with every factor of degree above 1

        It's also raised if a factor has degree at most `(t - 1)/2`::

            sage: F = GF(2^5)
            sage: R.<x> = F[]
            sage: g = (x^2 + x + 1) * (x^3 + x + 1)
            sage: g.is_squarefree(), g.is_irreducible()
            (True, False)
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: codes.decoders.GoppaCodePattersonDecoder(C)
            Traceback (most recent call last):
            ...
            ValueError: Patterson decoding needs a squarefree generating polynomial with every factor of degree above 2

        So the code in this module's first example is refused, since `x^9 + 1`
        splits into linear factors::

            sage: F = GF(2^6)
            sage: R.<x> = F[]
            sage: g = x^9 + 1
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: codes.decoders.GoppaCodePattersonDecoder(C)
            Traceback (most recent call last):
            ...
            ValueError: Patterson decoding needs a squarefree generating polynomial with every factor of degree above 4

        A constant polynomial is refused too::

            sage: F = GF(2^4)
            sage: R.<x> = F[]
            sage: C = codes.GoppaCode(R.one(), F.list())
            sage: codes.decoders.GoppaCodePattersonDecoder(C)
            Traceback (most recent call last):
            ...
            ValueError: Patterson decoding needs a nonconstant generating polynomial
        """
        # Refuse every code Patterson can't decode correctly
        if not isinstance(code, GoppaCode):
            raise ValueError("Patterson decoding needs a Goppa code")

        # Binary only
        # The square root step needs characteristic 2
        field = code.base_field()
        if field.cardinality() != 2:
            raise ValueError("Patterson decoding needs a binary Goppa code, "
                             f"not one over {field}")

        # factor() of a constant is empty, so the check below wouldn't refuse it
        t = code._generating_pol.degree()
        if t < 1:
            raise ValueError("Patterson decoding needs a nonconstant "
                             "generating polynomial")

        # Squarefree always, since square roots mod g need it
        if not check:
            if not code._generating_pol.is_squarefree():
                raise ValueError("Patterson decoding needs a squarefree "
                                 "generating polynomial")
        # No factor small enough to divide b in sigma = a^2 + x b^2,
        # otherwise a correctable word's syndrome can have no inverse
        elif any(multiplicity > 1 or factor.degree() <= (t - 1) // 2
                 for factor, multiplicity in code._generating_pol.factor()):
            raise ValueError("Patterson decoding needs a squarefree generating "
                             "polynomial with every factor of degree above "
                             f"{(t - 1) // 2}")

        super().__init__(code, code.ambient_space(), "GoppaEncoder")

        # Store 1/(x - a) mod g for every point a of the defining set
        # Inverted once here, since every decode sums about n/2 of them
        generating_pol = code._generating_pol
        x = generating_pol.parent().gen()
        self._syndrome_terms = [(x - point).inverse_mod(generating_pol)
                                for point in code._defining_set]

        # Store sqrt(x) mod g, which every square root is built from
        # g = g0^2 + x g1^2, so x = (g0 / g1)^2 mod g
        # (g1 has an inverse because g is squarefree)
        g0, g1 = _split_even_odd(generating_pol)
        self._sqrt_x = (g0 * g1.inverse_mod(generating_pol)) % generating_pol

        # Map each point of the defining set to its position in the word
        # (turns the error locator's roots into the positions to flip)
        self._positions = {point: i
                           for i, point in enumerate(code._defining_set)}

    def _square_root(self, f: Polynomial) -> Polynomial:
        r"""
        Return the square root of ``f`` modulo `g`.

        Every polynomial has exactly one square root modulo `g`. Writing
        `f = f_0^2 + x f_1^2`, it is `f_0 + \sqrt{x} f_1`. The decoder computes
        `\sqrt{x}` once, as `g_0 g_1^{-1}` from `g = g_0^2 + x g_1^2`
        [Hub1996]_.

        INPUT:

        - ``f`` -- a polynomial over the field of `g`

        OUTPUT: a polynomial of degree less than the degree of `g`

        TESTS:

        The stored `\sqrt{x}` squares to `x`::

            sage: F = GF(2^3)
            sage: R.<x> = F[]
            sage: g = x^2 + x + 1
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: D = codes.decoders.GoppaCodePattersonDecoder(C)
            sage: (D._sqrt_x^2) % g == x
            True

        A square root squares back to the polynomial, for any polynomial::

            sage: f = R.random_element(degree=5)
            sage: (D._square_root(f)^2 - f) % g == 0
            True

        It agrees with squaring `m t - 1` times, which gives the square root
        since `f^{2^{m t}} = f` modulo `g`::

            sage: f = R.random_element(degree=1)
            sage: root = f
            sage: for _ in range(3 * 2 - 1):
            ....:     root = root^2 % g
            sage: D._square_root(f) == root
            True

        The same holds on a larger code::

            sage: F = GF(2^5)
            sage: R.<x> = F[]
            sage: g = R.irreducible_element(6)
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: D = codes.decoders.GoppaCodePattersonDecoder(C)
            sage: f = R.random_element(degree=11)
            sage: (D._square_root(f)^2 - f) % g == 0
            True
        """
        # Split f into f0^2 + x f1^2, then sqrt(f) = f0 + sqrt(x) f1 mod g
        f0, f1 = _split_even_odd(f)
        return (f0 + self._sqrt_x * f1) % self.code()._generating_pol

    def _syndrome(self, word: FreeModuleElement) -> Polynomial:
        r"""
        Return the syndrome polynomial of ``word``.

        This is the sum of `1/(x - \alpha_i)` modulo `g`, over the positions
        `i` where ``word`` has a 1. Here `\alpha_i` is the `i`-th element of
        the defining set.

        A codeword's syndrome is zero, so the syndrome of a received word only
        depends on its errors.

        INPUT:

        - ``word`` -- a vector in the ambient space of the code

        OUTPUT: a polynomial of degree less than the degree of `g`

        TESTS:

        A codeword has syndrome zero::

            sage: F = GF(2^3)
            sage: R.<x> = F[]
            sage: g = x^2 + x + 1
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: D = codes.decoders.GoppaCodePattersonDecoder(C)
            sage: c = vector(GF(2), (0, 1, 1, 1, 1, 1, 1, 0))
            sage: c in C
            True
            sage: D._syndrome(c)
            0

        A single error at position `i` has syndrome `1/(x - \alpha_i)`::

            sage: e = C.ambient_space().gen(3)
            sage: e
            (0, 0, 0, 1, 0, 0, 0, 0)
            sage: (D._syndrome(e) * (x - L[3])) % g
            1

        It agrees with the syndrome built from the field-valued syndromes of a
        word `y`, `s_j = \sum_i y_i \alpha_i^j / g(\alpha_i)`, the entries of a
        parity check matrix before they are written out in bits::

            sage: t = g.degree()
            sage: y = c + e
            sage: s = [sum(y[i] * L[i]^j / g(L[i]) for i in range(len(L)))
            ....:      for j in range(t)]
            sage: S = sum(s[j] * (g // x^(j + 1)) for j in range(t))
            sage: D._syndrome(y) == S
            True

        And it is zero exactly when the word is in the code::

            sage: H = C.parity_check_matrix()
            sage: all((D._syndrome(w) == 0) == (H * w == 0) for w in [c, y, e])
            True
        """
        # Add up the stored terms where the word has a 1
        terms = (self._syndrome_terms[i] for i, bit in enumerate(word) if bit)
        return sum(terms, self.code()._generating_pol.parent().zero())

    def decode_to_code(self, word: FreeModuleElement) -> FreeModuleElement:
        r"""
        Correct the errors in ``word`` and return a codeword.

        If ``word`` has more than `t` errors, either a codeword at distance at
        most `t` from ``word`` is returned, or a
        :class:`~sage.coding.decoder.DecodingError` is raised. A word that is
        not in the ambient space of the code is refused.

        INPUT:

        - ``word`` -- a vector in the ambient space of the code

        OUTPUT: a vector of the code

        EXAMPLES:

        The code knows the decoder by name::

            sage: F = GF(2^4)
            sage: R.<x> = F[]
            sage: g = R.irreducible_element(3)
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: 'Patterson' in C.decoders_available()
            True
            sage: D = C.decoder('Patterson')
            sage: D
            Patterson decoder for [16, 4] Goppa code over GF(2)

        Any `t` errors are corrected::

            sage: t = D.decoding_radius()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), t)
            sage: c = C.random_element()
            sage: D.decode_to_code(Chan(c)) == c
            True

        So a message comes back too, through the code's encoder::

            sage: m = random_vector(GF(2), C.dimension())
            sage: D.decode_to_message(Chan(C.encode(m))) == m
            True

        A codeword comes back unchanged::

            sage: D.decode_to_code(c) == c
            True

        TESTS:

        With more than `t` errors, the result is a codeword within distance
        `t`, or an error::

            sage: from sage.coding.decoder import DecodingError
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(),
            ....:                                        t + 1)
            sage: y = Chan(c)
            sage: try:
            ....:     z = D.decode_to_code(y)
            ....:     print(z in C and (z - y).hamming_weight() <= t)
            ....: except DecodingError:
            ....:     print(True)
            True

        A word too far from every codeword is refused::

            sage: F = GF(2^3)
            sage: R.<x> = F[]
            sage: g = x^2 + x + 1
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: D = codes.GoppaCode(g, L).decoder('Patterson')
            sage: D.decode_to_code(vector(GF(2), (0, 0, 0, 0, 1, 0, 1, 1)))
            Traceback (most recent call last):
            ...
            DecodingError: Patterson decoding failed, the word has more than 2 errors

        One error at the point 0 of the defining set is corrected, wherever 0
        sits in the defining set::

            sage: F = GF(2^3)
            sage: R.<x> = F[]
            sage: g = x^2 + x + 1
            sage: L = [a for a in F.list() if a != 0] + [F.zero()]
            sage: C = codes.GoppaCode(g, L)
            sage: D = C.decoder('Patterson')
            sage: c = C.random_element()
            sage: D.decode_to_code(c + C.ambient_space().gen(7)) == c
            True

        The word passed in is left as it was::

            sage: y = c + C.ambient_space().gen(7)
            sage: D.decode_to_code(y) == c
            True
            sage: y == c + C.ambient_space().gen(7)
            True

        A codeword passed in comes back as a copy, not as the same vector::

            sage: D.decode_to_code(c) is c
            False

        A word Sage can convert into the ambient space is decoded as that
        converted word::

            sage: w = vector(ZZ, [2, 0, 0, 0, 0, 0, 0, 0])
            sage: z = D.decode_to_code(w)
            sage: z, z.parent() == C.ambient_space()
            ((0, 0, 0, 0, 0, 0, 0, 0), True)
            sage: w = vector(GF(4), c + C.ambient_space().gen(2))
            sage: z = D.decode_to_code(w)
            sage: z == c, z.parent() == C.ambient_space()
            (True, True)

        A word of the wrong length, or over the wrong field, is refused::

            sage: D.decode_to_code(vector(GF(2), 5))
            Traceback (most recent call last):
            ...
            ValueError: Patterson decoding needs a word in Vector space of dimension 8 over Finite Field of size 2
            sage: D.decode_to_code(vector(GF(3), 8))
            Traceback (most recent call last):
            ...
            ValueError: Patterson decoding needs a word in Vector space of dimension 8 over Finite Field of size 2

        A word past `t` errors whose syndrome has no inverse is refused too,
        which only happens for a reducible `g`::

            sage: F = GF(2^4)
            sage: R.<x> = F[]
            sage: g = x * (x + 1)
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: D = codes.decoders.GoppaCodePattersonDecoder(C)
            sage: D.decode_to_code(vector(GF(2), [1]*7 + [0]*7))
            Traceback (most recent call last):
            ...
            DecodingError: Patterson decoding failed, the word has more than 2 errors
        """
        # Refuse a word of the wrong length or over the wrong field
        if word not in self.input_space():
            raise ValueError("Patterson decoding needs a word in "
                             f"{self.input_space()}")
        # Work on the word as a vector over GF(2)
        # (a ZZ or GF(4) word can pass the check, and a 2 there is a 0)
        word = self.input_space()(word)

        generating_pol = self.code()._generating_pol
        t = generating_pol.degree()
        x = generating_pol.parent().gen()

        # Get the word's syndrome (zero means it's already a codeword)
        syndrome = self._syndrome(word)
        if syndrome.is_zero():
            return copy(word)

        # Build tau = sqrt(1/S + x) mod g
        # (with at most t errors, a nonzero syndrome always has an inverse mod g)
        try:
            syndrome_inverse = syndrome.inverse_mod(generating_pol)
        except ArithmeticError:
            raise DecodingError("Patterson decoding failed, the word has more "
                                f"than {t} errors")
        tau = self._square_root(syndrome_inverse + x)

        # Solve the key equation a = b tau mod g, deg a <= t/2, deg b <= (t-1)/2
        # Run Euclid on g and tau, stopping once the remainder is small enough
        # (remainder = coefficient * tau mod g holds the whole way)
        remainder_before, remainder = generating_pol, tau
        coefficient_before, coefficient = x.parent().zero(), x.parent().one()
        while remainder.degree() > t // 2:
            quotient, rest = remainder_before.quo_rem(remainder)
            remainder_before, remainder = remainder, rest
            coefficient_before, coefficient = (
                coefficient, coefficient_before - quotient * coefficient)
        a, b = remainder, coefficient

        # Error locator sigma = a^2 + x b^2
        # Its roots are the points of the defining set where the errors are
        error_locator = a**2 + x * b**2
        roots = error_locator.roots(multiplicities=False)
        positions = [self._positions.get(root) for root in roots]

        # Refuse unless sigma's roots are distinct and all in the defining set
        # The above is true if there are at most t errors
        # With more it can still be true, then a different codeword is returned
        # It's still within t of the word, since S sigma = sigma' = b^2 mod g
        if len(positions) != error_locator.degree() or None in positions:
            raise DecodingError("Patterson decoding failed, the word has more "
                                f"than {t} errors")

        # Flip the bits at the error positions
        corrected = copy(word)
        for i in positions:
            corrected[i] += 1
        return corrected

    def _repr_(self) -> str:
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: F = GF(2^3)
            sage: R.<x> = F[]
            sage: g = x^2 + x + 1
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: D = codes.decoders.GoppaCodePattersonDecoder(C)
            sage: D
            Patterson decoder for [8, 2] Goppa code over GF(2)
        """
        return "Patterson decoder for {}".format(self.code())

    def _latex_(self) -> str:
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: F = GF(2^3)
            sage: R.<x> = F[]
            sage: g = x^2 + x + 1
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: D = codes.decoders.GoppaCodePattersonDecoder(C)
            sage: latex(D)
            \text{Patterson decoder for }[8, 2]\text{ Goppa code over }\Bold{F}_{2}
        """
        return r"\text{{Patterson decoder for }}{}".format(
            self.code()._latex_())

    def __eq__(self, other: object) -> bool:
        r"""
        Test equality with ``other``.

        Two Patterson decoders are equal when their codes are.

        EXAMPLES::

            sage: F = GF(2^3)
            sage: R.<x> = F[]
            sage: g = x^2 + x + 1
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: D1 = codes.decoders.GoppaCodePattersonDecoder(C)
            sage: D2 = codes.decoders.GoppaCodePattersonDecoder(C)
            sage: D1 == D2
            True
            sage: D1 is D2
            False
            sage: g2 = x^2 + x + F.gen()^3
            sage: L2 = [a for a in F.list() if g2(a) != 0]
            sage: C3 = codes.GoppaCode(g2, L2)
            sage: D3 = codes.decoders.GoppaCodePattersonDecoder(C3)
            sage: D1 == D3
            False

        A decoder is never equal to something that is not a decoder::

            sage: D1 == 1
            False

        Equal decoders hash the same, so a set keeps one of them::

            sage: len({D1, D2, D3})
            2
        """
        return (isinstance(other, GoppaCodePattersonDecoder)
                and self.code() == other.code())

    def __hash__(self) -> int:
        r"""
        Return the hash value of ``self``.

        Equal decoders have equal hashes, as both depend only on the
        generating polynomial and the defining set of the code.

        EXAMPLES::

            sage: F = GF(2^3)
            sage: R.<x> = F[]
            sage: g = x^2 + x + 1
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: D1 = codes.decoders.GoppaCodePattersonDecoder(C)
            sage: D2 = codes.decoders.GoppaCodePattersonDecoder(C)
            sage: hash(D1) == hash(D2)
            True
            sage: len({D1, D2})
            1
        """
        # The code itself isn't hashable, so hash what its equality compares
        code = self.code()
        return hash((code._generating_pol, tuple(code._defining_set)))

    def decoding_radius(self) -> int:
        r"""
        Return the maximal number of errors that ``self`` can decode.

        This is `t`, the degree of the generating polynomial.

        EXAMPLES::

            sage: F = GF(2^3)
            sage: R.<x> = F[]
            sage: g = x^2 + x + 1
            sage: L = [a for a in F.list() if g(a) != 0]
            sage: C = codes.GoppaCode(g, L)
            sage: D = codes.decoders.GoppaCodePattersonDecoder(C)
            sage: D.decoding_radius()
            2

        The types of ``self``::

            sage: sorted(D.decoder_type())
            ['always-succeed', 'hard-decision']
        """
        # Degree of the generating polynomial, not the minimum distance
        # (computing the minimum distance takes exponential time)
        return self.code()._generating_pol.degree()


GoppaCode._registered_encoders["GoppaEncoder"] = GoppaCodeEncoder
GoppaCode._registered_decoders["Patterson"] = GoppaCodePattersonDecoder
GoppaCodePattersonDecoder._decoder_type = {"hard-decision", "always-succeed"}
