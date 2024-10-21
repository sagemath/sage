# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.combinat
"""
Classical Ciphers
"""
# ****************************************************************************
#       Copyright (C) 2007 David Kohel <kohel@maths.usyd.edu.au>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.misc.lazy_import import lazy_import
from sage.monoids.string_monoid_element import StringMonoidElement

lazy_import('sage.modules.free_module', 'FreeModule')

from .cipher import SymmetricKeyCipher


class AffineCipher(SymmetricKeyCipher):
    r"""
    Affine cipher class. This is the class that does the actual work of
    encryption and decryption. Users should not directly instantiate or
    create objects of this class. Instead, functionalities of this class
    should be accessed via
    :class:`AffineCryptosystem <sage.crypto.classical.AffineCryptosystem>`
    as the latter provides a convenient user interface.
    """

    def __init__(self, parent, key):
        r"""
        Create an affine cipher.

        INPUT:

        - ``parent`` -- an ``AffineCryptosystem`` object

        - ``key`` -- a secret key; let `N` be the size of the cipher domain.
          A key of this affine cipher is an ordered pair
          `(a, b) \in \ZZ_N \times \ZZ_N` such that `\gcd(a, N) = 1`.

        EXAMPLES:

        Testing of dumping and loading object::

            sage: A = AffineCryptosystem(AlphabeticStrings())
            sage: AC = A(3, 5)
            sage: AC == loads(dumps(AC))
            True
        """
        SymmetricKeyCipher.__init__(self, parent, key)

    def __eq__(self, other):
        r"""
        Comparing this ``AffineCipher`` with ``other``. Two ``AffineCipher``
        objects are the same if they are of the same type, have the same
        parent, and share the same secret key.

        INPUT:

        - ``other`` -- another object to compare with

        OUTPUT: ``True`` if ``self`` and ``other`` are the same
        ``AffineCipher`` object; ``False`` otherwise

        EXAMPLES::

            sage: aff1 = AffineCryptosystem(AlphabeticStrings())
            sage: aff2 = AffineCryptosystem(AlphabeticStrings())
            sage: aff1 == aff2
            True
            sage: aff1(1, 2) == aff2(1, 2)
            True
        """
        return type(self) is type(other) and self.parent() == other.parent() and self.key() == other.key()

    def __call__(self, M):
        r"""
        Return the ciphertext (respectively, plaintext) corresponding to
        ``M``. This is the main place where encryption and decryption takes
        place.

        INPUT:

        - ``M`` -- a message to be encrypted or decrypted. This message must
          be encoded using the plaintext or ciphertext alphabet. The current
          behaviour is that the plaintext and ciphertext alphabets are the
          same alphabet.

        - ``algorithm`` -- (default: ``'encrypt'``) whether to use the
          encryption or decryption algorithm on ``M``. The flag ``'encrypt'``
          signifies using the encryption algorithm, while ``'decrypt'``
          signifies using the decryption algorithm. The only acceptable
          values for ``algorithm`` are: ``'encrypt'`` and ``'decrypt'``.

        OUTPUT: the ciphertext or plaintext corresponding to ``M``

        EXAMPLES::

            sage: A = AffineCryptosystem(AlphabeticStrings()); A
            Affine cryptosystem on Free alphabetic string monoid on A-Z
            sage: P = A.encoding("The affine cryptosystem generalizes the shift cipher.")
            sage: P
            THEAFFINECRYPTOSYSTEMGENERALIZESTHESHIFTCIPHER
            sage: a, b = (9, 13)
            sage: E = A(a, b); E
            Affine cipher on Free alphabetic string monoid on A-Z
            sage: C = E(P); C
            CYXNGGHAXFKVSCJTVTCXRPXAXKNIHEXTCYXTYHGCFHSYXK
            sage: aInv, bInv = A.inverse_key(a, b)
            sage: D = A(aInv, bInv); D
            Affine cipher on Free alphabetic string monoid on A-Z
            sage: D(C)
            THEAFFINECRYPTOSYSTEMGENERALIZESTHESHIFTCIPHER
            sage: D(C) == P
            True
            sage: D(C) == P == D(E(P))
            True
        """
        # sanity check
        D = self.domain()  # = plaintext_space = ciphertext_space
        if M.parent() != D:
            raise TypeError("Argument M must be a string in the plaintext/ciphertext space.")

        from sage.rings.finite_rings.integer_mod import Mod
        A = list(D.alphabet())     # plaintext/ciphertext alphabet as a list
        N = self.domain().ngens()  # number of elements in this alphabet
        a, b = self.key()          # encryption/decryption key (a,b)
        # Let I be the index list of M. That is, the i-th element of M has
        # index k in the cipher domain D. We store this cipher domain index
        # as the i-th element of I.
        I = [A.index(str(e)) for e in M]

        # Encrypt the plaintext M. For each element i in I, the ciphertext
        # corresponding to i is ai + b (mod N). This can also be used for
        # decryption, in which case (a, b) is the inverse key corresponding
        # to a secret key.
        return D([ A.index(A[Mod(a*i + b, N).lift()]) for i in I ])

    def _repr_(self):
        r"""
        Return the string representation of this affine cipher.

        EXAMPLES::

            sage: A = AffineCryptosystem(AlphabeticStrings())
            sage: A(1, 2)
            Affine cipher on Free alphabetic string monoid on A-Z
        """
        # The affine cipher has the plaintext and ciphertext spaces defined
        # over the same non-empty alphabet. The cipher domain is the same
        # as the alphabet used for the plaintext and ciphertext spaces.
        return "Affine cipher on %s" % self.parent().cipher_domain()


class HillCipher(SymmetricKeyCipher):
    """
    Hill cipher class
    """
    def __init__(self, parent, key):
        """
        Create a Hill cipher.

        EXAMPLES::

            sage: # needs sage.modules
            sage: S = AlphabeticStrings()
            sage: E = HillCryptosystem(S,3); E
            Hill cryptosystem on Free alphabetic string monoid on A-Z of block length 3
            sage: M = E.key_space()
            sage: A = M([[1,0,1],[0,1,1],[2,2,3]]); A
            [1 0 1]
            [0 1 1]
            [2 2 3]
            sage: e = E(A); e
            Hill cipher on Free alphabetic string monoid on A-Z of block length 3
            sage: e(S("LAMAISONBLANCHE"))
            JYVKSKQPELAYKPV

        TESTS::

            sage: S = AlphabeticStrings()
            sage: E = HillCryptosystem(S,3)                                             # needs sage.modules
            sage: E == loads(dumps(E))                                                  # needs sage.modules
            True
        """
        # TODO: some type checking that the key is an invertible matrix?
        SymmetricKeyCipher.__init__(self, parent, key)

    def __eq__(self, right):
        return type(self) is type(right) and self.parent() == right.parent() and self.key() == right.key()

    def __call__(self, M):
        S = self.domain() # = plaintext_space = ciphertext_space
        if not isinstance(M, StringMonoidElement) and M.parent() == S:
            raise TypeError("Argument M (= %s) must be a string in the plaintext space." % M)
        m = self.parent().block_length()
        if len(M) % m != 0:
            raise TypeError("The length of M (= %s) must be a multiple of %s." % (M, m ))
        Alph = list(S.alphabet())
        A = self.key() # A is an m x m matrix
        R = A.parent().base_ring()
        V = FreeModule(R,m)
        Mstr = str(M)
        C = []
        for i in range(len(M)//m):
            v = V([ Alph.index(Mstr[m*i+j]) for j in range(m) ])
            C += (v * A).list()
        return S([ k.lift() for k in C ])

    def _repr_(self):
        r"""
        Return the string representation of this Hill cipher.

        EXAMPLES::

            sage: # needs sage.modules
            sage: H = HillCryptosystem(AlphabeticStrings(), 3)
            sage: M = MatrixSpace(IntegerModRing(26), 3, 3)
            sage: A = M([[1,0,1], [0,1,1], [2,2,3]])
            sage: e = H(A); e
            Hill cipher on Free alphabetic string monoid on A-Z of block length 3
        """
        return "Hill cipher on %s of block length %s" % (
            self.parent().cipher_domain(), self.parent().block_length() )

    def inverse(self):
        E = self.parent()
        try:
            B = E.inverse_key(self.key())
        except Exception:
            raise ValueError("Argument\n\n%s\n\nmust be an invertible cipher." % self)
        return E(B)


class ShiftCipher(SymmetricKeyCipher):
    r"""
    Shift cipher class. This is the class that does the actual work of
    encryption and decryption. Users should not directly instantiate or
    create objects of this class. Instead, functionalities of this class
    should be accessed via
    :class:`ShiftCryptosystem <sage.crypto.classical.ShiftCryptosystem>`
    as the latter provides a convenient user interface.
    """

    def __init__(self, parent, key):
        r"""
        Create a shift cipher.

        INPUT:

        - ``parent`` -- a ``ShiftCryptosystem`` object

        - ``key`` -- a secret key

        EXAMPLES::

            sage: S = ShiftCryptosystem(AlphabeticStrings()); S
            Shift cryptosystem on Free alphabetic string monoid on A-Z
            sage: P = S.encoding("The shift cryptosystem generalizes the Caesar cipher.")
            sage: P
            THESHIFTCRYPTOSYSTEMGENERALIZESTHECAESARCIPHER
            sage: K = 7
            sage: C = S.enciphering(K, P); C
            AOLZOPMAJYFWAVZFZALTNLULYHSPGLZAOLJHLZHYJPWOLY
            sage: S.deciphering(K, C)
            THESHIFTCRYPTOSYSTEMGENERALIZESTHECAESARCIPHER
            sage: S.deciphering(K, C) == P
            True
        """
        SymmetricKeyCipher.__init__(self, parent, key)

    def __eq__(self, other):
        r"""
        Comparing this ``ShiftCipher`` with ``other``. Two ``ShiftCipher``
        objects are the same if they are of the same type, have the same
        parent, and share the same secret key.

        INPUT:

        - ``other`` -- another object to compare with

        OUTPUT: ``True`` if ``self`` and ``other`` are the same ``ShiftCipher``
        object; ``False`` otherwise.

        EXAMPLES::

            sage: shift1 = ShiftCryptosystem(AlphabeticStrings())
            sage: shift2 = ShiftCryptosystem(AlphabeticStrings())
            sage: shift1 == shift2
            True
            sage: shift1 = ShiftCryptosystem(HexadecimalStrings())
            sage: shift2 = ShiftCryptosystem(BinaryStrings())
            sage: shift1 == shift2
            False
        """
        return type(self) is type(other) and self.parent() == other.parent() and self.key() == other.key()

    def __call__(self, M):
        r"""
        Return the ciphertext (respectively, plaintext) corresponding to
        ``M``. This is the main place where encryption and decryption takes
        place.

        INPUT:

        - ``M`` -- a message to be encrypted or decrypted. This message must
          be encoded using the plaintext or ciphertext alphabet. The current
          behaviour is that the plaintext and ciphertext alphabets are the
          same alphabet.

        OUTPUT: the ciphertext or plaintext corresponding to ``M``

        EXAMPLES:

        These are indirect doctests. Functionalities of this class are
        usually invoked from
        :class:`ShiftCryptosystem <sage.crypto.classical.ShiftCryptosystem>`.

        ::

            sage: S = ShiftCryptosystem(AlphabeticStrings())
            sage: S.enciphering(12, S.encoding("Stop shifting me."))
            EFABETURFUZSYQ
            sage: S = ShiftCryptosystem(HexadecimalStrings())
            sage: S.enciphering(7, S.encoding("Shift me now."))
            cadfd0ddeb97d4dc97d5d6ee95
            sage: S = ShiftCryptosystem(BinaryStrings())
            sage: S.enciphering(1, S.encoding("OK, enough shifting."))
            1011000010110100110100111101111110011010100100011001000010001010100110001001011111011111100011001001011110010110100110011000101110010110100100011001100011010001
        """
        dom = self.domain()  # = plaintext_space = ciphertext_space
        if not isinstance(M, StringMonoidElement) and M.parent() == dom:
            raise TypeError("Argument M (= %s) must be a string in the plaintext/ciphertext space." % M)
        from sage.rings.finite_rings.integer_mod import Mod
        A = list(dom.alphabet())   # plaintext/ciphertext alphabet as a list
        N = self.domain().ngens()  # number of elements in this alphabet
        K = self.key()             # encryption/decryption key
        # Here, M is a message encoded within the ciphertext/plaintext
        # alphabet of this shift cryptosystem. The list A above is a list of
        # all elements of this alphabet, each element being associated with
        # an index 0 <= i < n, where n is the size of A. Now get the index
        # of each character in the message M, storing all these indices
        # in the index list I.
        # TODO: the following two lines are coded with clarity and
        # readability in mind. It is possible to remove the overhead of
        # doing a list comprehension to get the index associated with each
        # element of M. For example, the next two lines can be crammed into
        # one list comprehension as follows:
        # return dom([ A.index(A[Mod(A.index(str(i)) + K, N).lift()]) for i in M ])
        # But it performs badly compared to the following two lines.
        I = [A.index(str(e)) for e in M]
        # Perform encryption/decryption on the whole message M, returning
        # the result as a string encoded in the alphabet A.
        return dom([ A.index(A[Mod(i + K, N).lift()]) for i in I ])

    def _repr_(self):
        r"""
        Return the string representation of this shift cipher.

        EXAMPLES::

            sage: S = ShiftCryptosystem(AlphabeticStrings())
            sage: S(3)
            Shift cipher on Free alphabetic string monoid on A-Z
            sage: S = ShiftCryptosystem(HexadecimalStrings())
            sage: S(5)
            Shift cipher on Free hexadecimal string monoid
            sage: S = ShiftCryptosystem(BinaryStrings())
            sage: S(1)
            Shift cipher on Free binary string monoid
        """
        # The shift cipher has the plaintext and ciphertext spaces defined
        # over the same non-empty alphabet. The cipher domain is the same
        # as the alphabet used for the plaintext and ciphertext spaces.
        return "Shift cipher on %s" % self.parent().cipher_domain()


class SubstitutionCipher(SymmetricKeyCipher):
    """
    Substitution cipher class
    """
    def __init__(self, parent, key):
        """
        Create a substitution cipher.

        EXAMPLES::

            sage: S = AlphabeticStrings()
            sage: E = SubstitutionCryptosystem(S)
            sage: E
            Substitution cryptosystem on Free alphabetic string monoid on A-Z
            sage: K = S([ 25-i for i in range(26) ])
            sage: K
            ZYXWVUTSRQPONMLKJIHGFEDCBA
            sage: e = E(K)
            sage: m = S("THECATINTHEHAT")
            sage: e(m)
            GSVXZGRMGSVSZG

        TESTS::

            sage: S = AlphabeticStrings()
            sage: E = SubstitutionCryptosystem(S)
            sage: E == loads(dumps(E))
            True
        """
        SymmetricKeyCipher.__init__(self, parent, key)

    def __eq__(self, right):
        return type(self) is type(right) and self.parent() == right.parent() and self.key() == right.key()

    def __call__(self, M):
        S = self.domain() # = plaintext_space = ciphertext_space
        if not isinstance(M, StringMonoidElement) and M.parent() == S:
            raise TypeError("Argument M (= %s) must be a string in the plaintext space." % M)
        A = list(S.alphabet())
        K = str(self.key()) # K is a string, while we want the indices:
        I = [ A.index(K[i]) for i in range(len(K)) ]
        Mstr = str(M)
        return S([ I[A.index(Mstr[i])] for i in range(len(Mstr)) ])

    def _repr_(self):
        r"""
        Return the string representation of this substitution cipher.

        EXAMPLES::

            sage: A = HexadecimalStrings(); S = SubstitutionCryptosystem(A)
            sage: K = A([16-i for i in range(16)]); S(K)
            Substitution cipher on Free hexadecimal string monoid
            sage: A = AlphabeticStrings(); S = SubstitutionCryptosystem(A)
            sage: K = A([26-i for i in range(26)]); S(K)
            Substitution cipher on Free alphabetic string monoid on A-Z
            sage: A = OctalStrings(); S = SubstitutionCryptosystem(A)
            sage: K = A([8-i for i in range(8)]); S(K)
            Substitution cipher on Free octal string monoid
            sage: A = BinaryStrings(); S = SubstitutionCryptosystem(A)
            sage: K = A([2-i for i in range(2)]); S(K)
            Substitution cipher on Free binary string monoid
            sage: A = Radix64Strings(); S = SubstitutionCryptosystem(A)
            sage: K = A([64-i for i in range(64)]); S(K)
            Substitution cipher on Free radix 64 string monoid
        """
        return "Substitution cipher on %s" % self.parent().cipher_domain()

    def inverse(self):
        E = self.parent()
        K = E.inverse_key(self.key())
        return E(K)


class TranspositionCipher(SymmetricKeyCipher):
    """
    Transition cipher class
    """
    def __init__(self, parent, key):
        """
        Create a transposition cipher.

        EXAMPLES::

            sage: # needs sage.groups
            sage: S = AlphabeticStrings()
            sage: E = TranspositionCryptosystem(S,14); E
            Transposition cryptosystem on
             Free alphabetic string monoid on A-Z of block length 14
            sage: K = [ 14-i for i in range(14) ]; K
            [14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
            sage: e = E(K)
            sage: m = S("THECATINTHEHAT")
            sage: e(m)
            TAHEHTNITACEHT

        EXAMPLES::

            sage: # needs sage.groups
            sage: S = AlphabeticStrings()
            sage: E = TranspositionCryptosystem(S,15)
            sage: m = S("THECATANDTHEHAT")
            sage: G = E.key_space(); G
            Symmetric group of order 15! as a permutation group
            sage: g = G([ 3, 2, 1, 6, 5, 4, 9, 8, 7, 12, 11, 10, 15, 14, 13 ])
            sage: e = E(g)
            sage: e(m)
            EHTTACDNAEHTTAH

        TESTS::

            sage: S = AlphabeticStrings()
            sage: E = TranspositionCryptosystem(S,14)                                   # needs sage.groups
            sage: E == loads(dumps(E))                                                  # needs sage.groups
            True
        """
        n = parent.block_length()
        if isinstance(key, list) and not len(key) == n:
            raise ValueError("key (= %s) must have block length %s" % (key, n))
        SymmetricKeyCipher.__init__(self, parent, key)

    def __call__(self, M, mode='ECB'):
        S = self.domain() # = plaintext_space = ciphertext_space
        if not isinstance(M, StringMonoidElement) and M.parent() == S:
            raise TypeError("Argument M (= %s) must be a string in the plaintext space." % M)
        if not mode == "ECB":
            raise NotImplementedError("Enciphering not implemented for mode (= %s) other than 'ECB'." % mode)
        g = self.key()
        N = len(M)
        m = self.parent().block_length()
        if not N % m == 0:
            raise TypeError("Argument M (= %s) must be a string of length k*%s." % (M, m))
        Melt = M._element_list # this uses the internal structure of string monoids
        # Caution: this is parsed as an outer loop in k and an inner loop in i:
        #     for k in range(N//m):
        #         for i in range(m):
        #             S([ Melt[g(i+1)-1+k*m]
        return S([ Melt[g(i+1)-1+k*m] for k in range(N//m) for i in range(m) ])

    def inverse(self):
        E = self.parent()
        K = E.inverse_key(self.key())
        return E(K)


class VigenereCipher(SymmetricKeyCipher):
    """
    Vigenere cipher class
    """
    def __init__(self, parent, key):
        """
        Create a Vigenere cipher.

        EXAMPLES::

            sage: S = AlphabeticStrings()
            sage: E = VigenereCryptosystem(S,11)
            sage: K = S("SHAKESPEARE")
            sage: e = E(K)
            sage: m = S("THECATINTHEHAT")
            sage: e(m)
            LOEMELXRTYIZHT

        TESTS::

            sage: S = AlphabeticStrings()
            sage: E = VigenereCryptosystem(S,11)
            sage: E == loads(dumps(E))
            True
        """
        SymmetricKeyCipher.__init__(self, parent, key)

    def __call__(self, M, mode='ECB'):
        S = self.domain() # = plaintext_space = ciphertext_space
        if not isinstance(M, StringMonoidElement) and M.parent() == S:
            raise TypeError("Argument M (= %s) must be a string in the plaintext space." % M)
        if not mode == "ECB":
            raise NotImplementedError("Enciphering not implemented for mode (= %s) other than 'ECB'." % mode)
        K = self.key()
        m = self.parent().period()
        n = S.ngens()
        # This uses the internal structure of string monoids
        Melt = M._element_list
        Kelt = K._element_list
        return S([ (Melt[i]+Kelt[i % m]) % n for i in range(len(M)) ])

    def inverse(self):
        E = self.parent()
        K = E.inverse_key(self.key())
        return E(K)
