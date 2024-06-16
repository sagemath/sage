# sage.doctest: needs sage.modules sage.rings.finite_rings
r"""
Goppa code

This module implements Goppa codes and an encoder for them.

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

AUTHORS:

- Filip Ion, Marketa Slukova (2019-06): initial version
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
from sage.coding.linear_code import AbstractLinearCode
from sage.coding.encoder import Encoder
from sage.modules.free_module_element import vector
from sage.coding.all import codes

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


GoppaCode._registered_encoders["GoppaEncoder"] = GoppaCodeEncoder
