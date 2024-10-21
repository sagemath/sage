# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.modules sage.rings.finite_rings
r"""
Kasami code

This module implements a construction for the extended Kasami codes.
The "regular" Kasami codes are obtained from truncating the extended version.

The extended Kasami code with parameters `(s,t)` is defined as

.. MATH::

    \{ v \in \GF{2}^s \mid
    \sum_{a \in \GF{s}} v_a =
    \sum_{a \in \GF{s}} a v_a =
    \sum_{a \in \GF{s}} a^{t+1} v_a = 0 \}


It follows that these are subfield subcodes of the code having those three
equations as parity checks.  The only valid parameters `s,t` are given by the
below, where `q` is a power of 2

* `s = q^{2j+1}`, `t = q^m` with `m \leq j` and `\gcd(m,2j+1) = 1`
* `s = q^2`, `t=q`

The coset graphs of the Kasami codes are distance-regular.  In particular, the
extended Kasami codes result in distance-regular graphs with intersection arrays

* `[q^{2j+1}, q^{2j+1} - 1, q^{2j+1} - q, q^{2j+1} - q^{2j} + 1;`
  `1, q, q^{2j} -1, q^{2j+1}]`
* `[q^2, q^2 - 1, q^2 - q, 1; 1, q, q^2 - 1, q^2]`

The Kasami codes result in distance-regular graphs with intersection arrays

* `[q^{2j+1} - 1, q^{2j+1} - q, q^{2j+1} - q^{2j} + 1; 1, q, q^{2j} -1]`
* `[q^2 - 1, q^2 - q, 1; 1, q, q^2 - 1]`

REFERENCES:

- [BCN1989]_ p. 358 for a definition.

- [Kas1966a]_

- [Kas1966b]_

- [Kas1971]_

AUTHORS:

- Ivo Maffei (2020-07-09): initial version
"""

# ****************************************************************************
#       Copyright (C) 2020 Ivo Maffei <ivomaffei@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.finite_rings.finite_field_constructor import GF
from sage.matrix.constructor import matrix
from sage.coding.linear_code import (AbstractLinearCode,
                                     LinearCodeGeneratorMatrixEncoder)
from sage.arith.misc import is_prime_power, gcd


class KasamiCode(AbstractLinearCode):
    r"""
    Representation of a Kasami Code.

    The extended Kasami code with parameters `(s,t)` is defined as

    .. MATH::

        \{ v \in \GF{2}^s \mid
        \sum_{a \in \GF{s}} v_a =
        \sum_{a \in \GF{s}} a v_a =
        \sum_{a \in \GF{s}} a^{t+1} v_a = 0 \}

    The only valid parameters `s,t` are given by the below,
    where `q` is a power of 2:

    * `s = q^{2j+1}`, `t = q^m` with `m \leq j` and `\gcd(m,2j+1) = 1`
    * `s = q^2`, `t=q`

    The Kasami code `(s,t)` is obtained from the extended
    Kasami code `(s,t)`, via truncation of all words.

    INPUT:

    - ``s``, ``t`` -- integer; the parameters of the Kasami code

    - ``extended`` -- boolean (default: ``True``); if set to ``True``,
      creates an extended Kasami code

    EXAMPLES::

        sage: codes.KasamiCode(16,4)
        [16, 9] Extended (16, 4)-Kasami code
        sage: _.minimum_distance()                                                      # needs sage.libs.gap
        4

        sage: codes.KasamiCode(8, 2, extended=False)
        [7, 1] (8, 2)-Kasami code

        sage: codes.KasamiCode(8,4)
        Traceback (most recent call last):
        ...
        ValueError: The parameters(=8,4) are invalid. Check the documentation

    The extended Kasami code is the extension of the Kasami code::

        sage: C = codes.KasamiCode(16, 4, extended=False)
        sage: Cext = C.extended_code()
        sage: D = codes.KasamiCode(16, 4, extended=True)
        sage: D.generator_matrix() == Cext.generator_matrix()
        True

    .. SEEALSO::

        :mod:`sage.coding.linear_code`.

    REFERENCES:

    For more information on Kasami codes and their use see [BCN1989]_
    or [Kas1966a]_, [Kas1966b]_, [Kas1971]_

    TESTS::

        sage: C1 = codes.KasamiCode(16, 4)
        sage: C2 = codes.KasamiCode(16, 4, extended=False)
        sage: C1.parameters() == C2.parameters()
        True
        sage: C1 == C2
        False
        sage: C1.minimum_distance() == C2.minimum_distance() + 1                        # needs sage.libs.gap
        True
        sage: C = codes.KasamiCode(4,2)
        sage: C.dimension()
        0
        sage: C.generator_matrix()
        []
    """

    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, s, t, extended=True):
        r"""
        Constructor for the ``KasamiCode`` class.

        TESTS::

            sage: codes.KasamiCode(64,8)
            [64, 54] Extended (64, 8)-Kasami code

            sage: codes.KasamiCode(64,8, extended=False)
            [63, 54] (64, 8)-Kasami code

            sage: codes.KasamiCode(3,5)
            Traceback (most recent call last):
            ...
            ValueError: The parameter t(=5) must be a power of 2
        """
        # Check validity of s and t
        (p,i) = is_prime_power(t,get_data=True)
        if p != 2:
            raise ValueError(f"The parameter t(={t}) must be a power of 2")

        if s != t*t:
            # then we must have s=q^{2j+1} and t = q^m
            (p,k) = is_prime_power(s,get_data=True)
            if p != 2:
                raise ValueError(f"The parameter s(={s}) must be a power of 2")

            # q= 2^l here l = gcd(k,i)
            l = gcd(i,k)
            q = 2**l
            m = i // l

            if (k//l) % 2 == 0:
                raise ValueError(
                    f"The parameter s(={s}) is invalid. Check the documentation")

            j = ((k//l) - 1) // 2

            # gcd(m,2*j+1) = gcd( i/l, k/l) = 1
            if m > j:
                raise ValueError(
                    (f"The parameters(={s},{t}) are invalid. "
                     "Check the documentation"))

        # s and t are valid!!!
        self._s = s
        self._t = t
        self._extended = extended

        length = s - 1
        if extended:
            length += 1

        super().__init__(GF(2), length, "GeneratorMatrix", "Syndrome")

    def parameters(self):
        r"""
        Return the parameters `s,t` of ``self``.

        EXAMPLES::

            sage: C = codes.KasamiCode(16, 4, extended=True)
            sage: C.parameters()
            (16, 4)
            sage: D = codes.KasamiCode(16, 4, extended=False)
            sage: D.parameters()
            (16, 4)
            sage: C = codes.KasamiCode(8,2)
            sage: C.parameters()
            (8, 2)
        """
        return (self._s,self._t)

    def __eq__(self, other):
        r"""
        Test equality between Kasami Code objects.

        Two Kasami codes are the same if they
        have the same `s,t` values and are both
        extended or both regular.

        EXAMPLES::

            sage: C1 = codes.KasamiCode(8,2)
            sage: C2 = codes.KasamiCode(8,2)
            sage: C1.__eq__(C2)
            True
        """
        # Check that s, t, extended values of both
        # objects are the same
        return isinstance(other, KasamiCode) \
                and self.parameters() == other.parameters() \
                and self._extended == other._extended

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: codes.KasamiCode(4,2,extended=True)
            [4, 0] Extended (4, 2)-Kasami code
        """
        ext = ""
        if self._extended:
            ext = " Extended"
        return "[%s, %s]%s (%s, %s)-Kasami code"\
                % (self.length(),self.dimension(), ext, self._s, self._t)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: C = codes.KasamiCode(16,4)
            sage: latex(C)
            [16, 9]\textnormal{ Extended} (16, 4)\textnormal{-Kasami code}
        """
        ext = ""
        if self._extended:
            ext = " Extended"
        return "[%s, %s]\\textnormal{%s} (%s, %s)\\textnormal{-Kasami code}"\
                % (self.length(), self.dimension(), ext, self._s, self._t)

    def generator_matrix(self):
        r"""
        Return a generator matrix of ``self``.

        EXAMPLES::

            sage: C = codes.KasamiCode(16, 4, extended=False)
            sage: C.generator_matrix()
            [1 0 0 0 0 0 0 0 0 1 0 0 1 1 1]
            [0 1 0 0 0 0 0 0 0 1 1 0 1 0 0]
            [0 0 1 0 0 0 0 0 0 0 1 1 0 1 0]
            [0 0 0 1 0 0 0 0 0 0 0 1 1 0 1]
            [0 0 0 0 1 0 0 0 0 0 0 0 1 1 0]
            [0 0 0 0 0 1 0 0 0 1 1 0 1 1 1]
            [0 0 0 0 0 0 1 0 0 0 1 1 0 1 1]
            [0 0 0 0 0 0 0 1 0 1 1 1 0 0 1]
            [0 0 0 0 0 0 0 0 1 1 0 1 0 0 0]

        ALGORITHM:

        We build the parity check matrix given by the three equations that
        the codewords must satisfy. Then we generate the parity check matrix
        over `\GF{2}` and from this the obtain the generator matrix for the
        extended Kasami codes.

        For the Kasami codes, we truncate the last column.

        TESTS::

            sage: C = codes.KasamiCode(4,2)
            sage: C.generator_matrix()
            []
            sage: C = codes.KasamiCode(8,2)
            sage: C.generator_matrix()
            [1 1 1 1 1 1 1 1]
            sage: C.minimum_distance()
            8
            sage: C = codes.KasamiCode(8, 2, extended=False)
            sage: C.generator_matrix()
            [1 1 1 1 1 1 1]
            sage: C.minimum_distance()
            7
            sage: C = codes.KasamiCode(4, 2, extended=False)
            sage: C.generator_matrix()
            []
            sage: C = codes.KasamiCode(16, 4, extended=False)
            sage: C.minimum_distance()
            3
        """
        from sage.misc.functional import log

        m = log(self._s, 2)
        F = GF(self._s)

        def exp(row):
            return matrix(F,
                          [x + [0]*(m - len(x)) for x in
                            [a.polynomial().list() for a in row]]).transpose()

        # Parity check matrix over GF(s)
        Hs = matrix(F, [[1]*self._s,
                        F.list(),
                        [a**(self._t + 1) for a in F]])

        # Parity check matrix over GF(2)
        H = matrix.block(GF(2), 3, 1, [exp(row) for row in Hs.rows()])

        # Generator matrix
        G = H.right_kernel().basis_matrix()

        if self._extended:
            return G

        newG = [v[:-1] for v in G]
        return matrix(GF(2), newG)


KasamiCode._registered_encoders[
    "GeneratorMatrix"] = LinearCodeGeneratorMatrixEncoder
