# sage.doctest: needs sage.combinat sage.rings.finite_rings
"""
Stream Ciphers
"""
#*****************************************************************************
#       Copyright (C) 2007 David Kohel <kohel@maths.usyd.edu.au>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .lfsr import lfsr_sequence
from .cipher import SymmetricKeyCipher
from sage.monoids.string_monoid_element import StringMonoidElement

class LFSRCipher(SymmetricKeyCipher):
    def __init__(self, parent, poly, IS):
        """
        Create a linear feedback shift register (LFSR) cipher.

        INPUT:

        - ``parent`` -- parent

        - ``poly`` -- connection polynomial

        - ``IS`` -- initial state

        EXAMPLES::

            sage: FF = FiniteField(2)
            sage: P.<x> = PolynomialRing(FF)
            sage: E = LFSRCryptosystem(FF)
            sage: E
            LFSR cryptosystem over Finite Field of size 2
            sage: IS = [ FF(a) for a in [0,1,1,1,0,1,1] ]
            sage: g = x^7 + x + 1
            sage: e = E((g,IS))
            sage: B = BinaryStrings()
            sage: m = B.encoding("THECATINTHEHAT")
            sage: e(m)
            0010001101111010111010101010001100000000110100010101011100001011110010010000011111100100100011001101101000001111
            sage: FF = FiniteField(2)
            sage: P.<x> = PolynomialRing(FF)
            sage: LFSR = LFSRCryptosystem(FF)
            sage: e = LFSR((x^2+x+1,[FF(0),FF(1)]))
            sage: B = e.domain()
            sage: m = B.encoding("The cat in the hat.")
            sage: e(m)
            00111001110111101011111001001101110101011011101000011001100101101011001000000011100101101010111100000101110100111111101100000101110101111010111101000011
            sage: m == e(e(m))
            True

        TESTS::

            sage: FF = FiniteField(2)
            sage: P.<x> = PolynomialRing(FF)
            sage: E = LFSRCryptosystem(FF)
            sage: E == loads(dumps(E))
            True
        """
        SymmetricKeyCipher.__init__(self, parent, key=(poly, IS))

    def __call__(self, M, mode='ECB'):
        r"""
        Generate key stream from the binary string ``M``.

        INPUT:

        - ``M`` -- a StringMonoidElement

        - ``mode`` -- ignored (default: ``'ECB'``)

        EXAMPLES::

            sage: k = GF(2)
            sage: P.<x> = PolynomialRing( k )
            sage: LFSR = LFSRCryptosystem( k )
            sage: e = LFSR((x^2+x+1,[k(0), k(1)]))
            sage: B = e.domain()
            sage: m = B.encoding('The cat in the hat.')
            sage: e(m)
            00111001110111101011111001001101110101011011101000011001100101101011001000000011100101101010111100000101110100111111101100000101110101111010111101000011
        """
        B = self.domain() # = plaintext_space = ciphertext_space
        if not isinstance(M, StringMonoidElement) and M.parent() == B:
            raise TypeError("Argument M (= %s) must be a string in the plaintext space." % M)
        (poly, IS) = self.key()
        n = B.ngens() # two for binary strings
        N = len(M)
        Melt = M._element_list
        Kelt = lfsr_sequence(poly.list(), IS, N)
        return B([ (Melt[i]+int(Kelt[i])) % n for i in range(N) ])

    def _repr_(self):
        r"""
        Return the string representation of this LFSR cipher.

        EXAMPLES::

            sage: FF = FiniteField(2)
            sage: P.<x> = PolynomialRing(FF)
            sage: LFSR = LFSRCryptosystem(FF)
            sage: IS_1 = [ FF(a) for a in [0,1,0,1,0,0,0] ]
            sage: e1 = LFSR((x^7 + x + 1,IS_1))
            sage: IS_2 = [ FF(a) for a in [0,0,1,0,0,0,1,0,1] ]
            sage: e2 = LFSR((x^9 + x^3 + 1,IS_2))
            sage: E = ShrinkingGeneratorCryptosystem()
            sage: e = E((e1,e2))
            sage: e.keystream_cipher()
            LFSR cipher on Free binary string monoid
        """
        return "LFSR cipher on %s" % self.domain()

    def connection_polynomial(self):
        """
        The connection polynomial defining the LFSR of the cipher.

        EXAMPLES::

            sage: k = GF(2)
            sage: P.<x> = PolynomialRing( k )
            sage: LFSR = LFSRCryptosystem( k )
            sage: e = LFSR((x^2+x+1,[k(0), k(1)]))
            sage: e.connection_polynomial()
            x^2 + x + 1
        """
        return self.key()[0]

    def initial_state(self):
        """
        The initial state of the LFSR cipher.

        EXAMPLES::

            sage: k = GF(2)
            sage: P.<x> = PolynomialRing( k )
            sage: LFSR = LFSRCryptosystem( k )
            sage: e = LFSR((x^2+x+1,[k(0), k(1)]))
            sage: e.initial_state()
            [0, 1]
        """
        return self.key()[1]

class ShrinkingGeneratorCipher(SymmetricKeyCipher):
    def __init__(self, parent, e1, e2):
        """
        Create a shrinking generator cipher.

        INPUT:

        - ``parent`` -- parent

        - ``poly`` -- connection polynomial

        - ``IS`` -- initial state

        EXAMPLES::

            sage: FF = FiniteField(2)
            sage: P.<x> = PolynomialRing(FF)
            sage: LFSR = LFSRCryptosystem(FF)
            sage: IS_1 = [ FF(a) for a in [0,1,0,1,0,0,0] ]
            sage: e1 = LFSR((x^7 + x + 1,IS_1))
            sage: IS_2 = [ FF(a) for a in [0,0,1,0,0,0,1,0,1] ]
            sage: e2 = LFSR((x^9 + x^3 + 1,IS_2))
            sage: E = ShrinkingGeneratorCryptosystem()
            sage: e = E((e1,e2))
            sage: e
            Shrinking generator cipher on Free binary string monoid
        """
        if not isinstance(e1, LFSRCipher):
            raise TypeError("Argument e1 (= %s) must be a LFSR cipher." % e1)
        if not isinstance(e2, LFSRCipher):
            raise TypeError("Argument e2 (= %s) must be a LFSR cipher." % e2)
        SymmetricKeyCipher.__init__(self, parent, key=(e1, e2))

    def keystream_cipher(self):
        """
        The LFSR cipher generating the output key stream.

        EXAMPLES::

            sage: FF = FiniteField(2)
            sage: P.<x> = PolynomialRing(FF)
            sage: LFSR = LFSRCryptosystem(FF)
            sage: IS_1 = [ FF(a) for a in [0,1,0,1,0,0,0] ]
            sage: e1 = LFSR((x^7 + x + 1,IS_1))
            sage: IS_2 = [ FF(a) for a in [0,0,1,0,0,0,1,0,1] ]
            sage: e2 = LFSR((x^9 + x^3 + 1,IS_2))
            sage: E = ShrinkingGeneratorCryptosystem()
            sage: e = E((e1,e2))
            sage: e.keystream_cipher()
            LFSR cipher on Free binary string monoid
        """
        return self.key()[0]

    def decimating_cipher(self):
        """
        The LFSR cipher generating the decimating key stream.

        EXAMPLES::

            sage: FF = FiniteField(2)
            sage: P.<x> = PolynomialRing(FF)
            sage: LFSR = LFSRCryptosystem(FF)
            sage: IS_1 = [ FF(a) for a in [0,1,0,1,0,0,0] ]
            sage: e1 = LFSR((x^7 + x + 1,IS_1))
            sage: IS_2 = [ FF(a) for a in [0,0,1,0,0,0,1,0,1] ]
            sage: e2 = LFSR((x^9 + x^3 + 1,IS_2))
            sage: E = ShrinkingGeneratorCryptosystem()
            sage: e = E((e1,e2))
            sage: e.decimating_cipher()
            LFSR cipher on Free binary string monoid
        """
        return self.key()[1]

    def __call__(self, M, mode='ECB'):
        r"""
        INPUT:

        - ``M`` -- a StringMonoidElement

        - ``mode`` -- ignored (default: ``'ECB'``)

        EXAMPLES::

            sage: FF = FiniteField(2)
            sage: P.<x> = PolynomialRing(FF)
            sage: LFSR = LFSRCryptosystem(FF)
            sage: IS_1 = [ FF(a) for a in [0,1,0,1,0,0,0] ]
            sage: e1 = LFSR((x^7 + x + 1,IS_1))
            sage: IS_2 = [ FF(a) for a in [0,0,1,0,0,0,1,0,1] ]
            sage: e2 = LFSR((x^9 + x^3 + 1,IS_2))
            sage: E = ShrinkingGeneratorCryptosystem()
            sage: e = E((e1,e2))
            sage: B = BinaryStrings()
            sage: m = B.encoding("THECATINTHEHAT")
            sage: c = e(m)
            sage: c.decoding()
            "t\xb6\xc1'\x83\x17\xae\xc9ZO\x84V\x7fX"
            sage: e(e(m)) == m
            True
            sage: m.decoding()
            'THECATINTHEHAT'
        """
        B = self.domain() # = plaintext_space = ciphertext_space
        if not isinstance(M, StringMonoidElement) and M.parent() == B:
            raise TypeError("Argument M (= %s) must be a string in the plaintext space." % M)
        (e1, e2) = self.key()
        MStream = M._element_list
        g1 = e1.connection_polynomial()
        n1 = g1.degree()
        IS_1 = e1.initial_state()
        g2 = e2.connection_polynomial()
        n2 = g2.degree()
        IS_2 = e2.initial_state()
        k = 0
        N = len(M)
        n = max(n1, n2)
        CStream = []
        while k < N:
            r = max(N-k,2*n)
            KStream = lfsr_sequence(g1.list(), IS_1, r)
            DStream = lfsr_sequence(g2.list(), IS_2, r)
            for i in range(r - n):
                if DStream[i] != 0:
                    CStream.append(int(MStream[k] + KStream[i]))
                    k += 1
                if k == N:
                    break
            IS_1 = KStream[r-n-1:r-n+n1]
            IS_2 = DStream[r-n-1:r-n+n2]
        return B(CStream)

    def _repr_(self):
        r"""
        Return the string representation of this shrinking generator cipher.

        EXAMPLES::

            sage: FF = FiniteField(2)
            sage: P.<x> = PolynomialRing(FF)
            sage: LFSR = LFSRCryptosystem(FF)
            sage: IS_1 = [ FF(a) for a in [0,1,0,1,0,0,0] ]
            sage: e1 = LFSR((x^7 + x + 1,IS_1))
            sage: IS_2 = [ FF(a) for a in [0,0,1,0,0,0,1,0,1] ]
            sage: e2 = LFSR((x^9 + x^3 + 1,IS_2))
            sage: E = ShrinkingGeneratorCryptosystem()
            sage: e = E((e1,e2)); e
            Shrinking generator cipher on Free binary string monoid
        """
        return "Shrinking generator cipher on %s" % self.domain()
