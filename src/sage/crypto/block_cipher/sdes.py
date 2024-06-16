# sage.doctest: needs sage.combinat sage.rings.finite_rings
r"""
Simplified DES

A simplified variant of the Data Encryption Standard (DES). Note that
Simplified DES or S-DES is for educational purposes only. It is a
small-scale version of the DES designed to help beginners understand the
basic structure of DES.

AUTHORS:

- Minh Van Nguyen (2009-06): initial version
"""

###########################################################################
# Copyright (c) 2009 Minh Van Nguyen <nguyenminh2@gmail.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# http://www.gnu.org/licenses/
###########################################################################

from sage.monoids.string_monoid import BinaryStrings
from sage.structure.sage_object import SageObject

class SimplifiedDES(SageObject):
    r"""
    This class implements the Simplified Data Encryption Standard (S-DES)
    described in [Sch1996]_. Schaefer's S-DES is for educational purposes
    only and is not secure for practical purposes. S-DES is a version of
    the DES with all parameters significantly reduced, but at the same time
    preserving the structure of DES. The goal of S-DES is to allow a
    beginner to understand the structure of DES, thus laying a foundation
    for a thorough study of DES. Its goal is as a teaching tool in the same
    spirit as Phan's
    :mod:`Mini-AES <sage.crypto.block_cipher.miniaes>` [Pha2002]_.

    EXAMPLES:

    Encrypt a random block of 8-bit plaintext using a random key, decrypt
    the ciphertext, and compare the result with the original plaintext::

        sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
        sage: sdes = SimplifiedDES(); sdes
        Simplified DES block cipher with 10-bit keys
        sage: bin = BinaryStrings()
        sage: P = [bin(str(randint(0, 1))) for i in range(8)]
        sage: K = sdes.random_key()
        sage: C = sdes.encrypt(P, K)
        sage: plaintxt = sdes.decrypt(C, K)
        sage: plaintxt == P
        True

    We can also encrypt binary strings that are larger than 8 bits in length.
    However, the number of bits in that binary string must be positive
    and a multiple of 8::

        sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
        sage: sdes = SimplifiedDES()
        sage: bin = BinaryStrings()
        sage: P = bin.encoding("Encrypt this using S-DES!")
        sage: Mod(len(P), 8) == 0
        True
        sage: K = sdes.list_to_string(sdes.random_key())
        sage: C = sdes(P, K, algorithm='encrypt')
        sage: plaintxt = sdes(C, K, algorithm='decrypt')
        sage: plaintxt == P
        True
    """

    def __init__(self):
        r"""
        A simplified variant of the Data Encryption Standard (DES).

        EXAMPLES::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES(); sdes
            Simplified DES block cipher with 10-bit keys
            sage: B = BinaryStrings()
            sage: P = [B(str(randint(0, 1))) for i in range(8)]
            sage: K = sdes.random_key()
            sage: C = sdes.encrypt(P, K)
            sage: plaintxt = sdes.decrypt(C, K)
            sage: plaintxt == P
            True
        """
        from sage.crypto.sbox import SBox
        # the number of bits in a secret key
        self._key_size = 10
        # the S-box S_0
        self._sbox0 = SBox(1, 0, 3, 2, 3, 2, 1, 0, 0, 2, 1, 3, 3, 1, 3, 2)
        # the S-box S_1
        self._sbox1 = SBox(0, 1, 2, 3, 2, 0, 1, 3, 3, 0, 1, 0, 2, 1, 0, 3)

    def __call__(self, B, K, algorithm='encrypt'):
        r"""
        Apply S-DES encryption or decryption on the binary string ``B``
        using the key ``K``.  The flag ``algorithm`` controls what action is
        to be performed on ``B``.

        INPUT:

        - ``B`` -- a binary string, where the number of bits is positive and
          a multiple of 8

        - ``K`` -- a secret key; this must be a 10-bit binary string

        - ``algorithm`` -- (default: ``'encrypt'``) a string; a flag to signify
          whether encryption or decryption is to be applied to the binary
          string ``B``. The encryption flag is ``'encrypt'`` and the decryption
          flag is ``'decrypt'``.

        OUTPUT:

        - The ciphertext (respectively plaintext) corresponding to the
          binary string ``B``.

        EXAMPLES:

        Encrypt a plaintext, decrypt the ciphertext, and compare the
        result with the original plaintext::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: bin = BinaryStrings()
            sage: P = bin.encoding("Encrypt this using DES!")
            sage: K = sdes.random_key()
            sage: K = sdes.list_to_string(K)
            sage: C = sdes(P, K, algorithm='encrypt')
            sage: plaintxt = sdes(C, K, algorithm='decrypt')
            sage: plaintxt == P
            True

        TESTS:

        The binary string ``B`` must be non-empty and the number of bits must
        be a multiple of 8::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: sdes("B", "K")
            Traceback (most recent call last):
            ...
            TypeError: input B must be a non-empty binary string with number of bits a multiple of 8
            sage: bin = BinaryStrings()
            sage: B = bin("101")
            sage: sdes(B, "K")
            Traceback (most recent call last):
            ...
            ValueError: the number of bits in the binary string B must be positive and a multiple of 8

        The secret key ``K`` must be a block of 10 bits::

            sage: B = bin.encoding("abc")
            sage: sdes(B, "K")
            Traceback (most recent call last):
            ...
            TypeError: secret key must be a 10-bit binary string
            sage: K = bin("1010")
            sage: sdes(B, K)
            Traceback (most recent call last):
            ...
            ValueError: secret key must be a 10-bit binary string

        The value for ``algorithm`` must be either ``'encrypt'`` or
        ``'decrypt'``::

            sage: B = bin.encoding("abc")
            sage: K = sdes.list_to_string(sdes.random_key())
            sage: sdes(B, K, algorithm='e')
            Traceback (most recent call last):
            ...
            ValueError: algorithm must be either 'encrypt' or 'decrypt'
            sage: sdes(B, K, algorithm='d')
            Traceback (most recent call last):
            ...
            ValueError: algorithm must be either 'encrypt' or 'decrypt'
            sage: sdes(B, K, algorithm='abc')
            Traceback (most recent call last):
            ...
            ValueError: algorithm must be either 'encrypt' or 'decrypt'
        """
        from sage.monoids.string_monoid_element import StringMonoidElement
        from sage.rings.finite_rings.integer_mod import Mod
        # S-DES operates on 8-bit ciphertext/plaintext blocks
        Blength = 8

        if not isinstance(B, StringMonoidElement):
            raise TypeError("input B must be a non-empty binary string with number of bits a multiple of 8")
        if (len(B) == 0) or (Mod(len(B), Blength).lift() != 0):
            raise ValueError("the number of bits in the binary string B must be positive and a multiple of 8")
        if not isinstance(K, StringMonoidElement):
            raise TypeError("secret key must be a 10-bit binary string")
        if len(K) != self._key_size:
            raise ValueError("secret key must be a 10-bit binary string")

        N = len(B) // Blength  # the number of 8-bit blocks
        S = ""
        bin = BinaryStrings()
        # encrypt each 8-bit block in succession
        if algorithm == "encrypt":
            for i in range(N):
                # get an 8-bit block
                block = B[i*Blength : (i+1)*Blength]
                block = self.string_to_list(str(block))
                key = self.string_to_list(str(K))
                # encrypt the block using key
                C = self.encrypt(block, key)
                C = self.list_to_string(C)
                # append encrypted block to ciphertext string
                S = "".join([S, str(C)])
            return bin(S)
        # decrypt each 8-bit block in succession
        elif algorithm == "decrypt":
            for i in range(N):
                # get an 8-bit block
                block = B[i*Blength : (i+1)*Blength]
                block = self.string_to_list(str(block))
                key = self.string_to_list(str(K))
                # decrypt the block using key
                P = self.decrypt(block, key)
                P = self.list_to_string(P)
                # append decrypted block to plaintext string
                S = "".join([S, str(P)])
            return bin(S)
        # invalid value for algorithm option
        else:
            raise ValueError("algorithm must be either 'encrypt' or 'decrypt'")

    def __eq__(self, other):
        r"""
        Compare ``self`` with ``other``.

        Simplified DES objects are the same if they have the same key size
        and S-boxes.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: s = SimplifiedDES()
            sage: s == loads(dumps(s))
            True
        """
        return ( (self._key_size == other._key_size) and
                 (self._sbox0 == other._sbox0) and
                 (self._sbox1 == other._sbox1) )

    def __repr__(self):
        r"""
        A string representation of this Simplified DES.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: SimplifiedDES()
            Simplified DES block cipher with 10-bit keys
        """
        return "Simplified DES block cipher with 10-bit keys"

    def block_length(self):
        r"""
        Return the block length of Schaefer's S-DES block cipher. A key in
        Schaefer's S-DES is a block of 10 bits.

        OUTPUT:

        - The block (or key) length in number of bits.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: sdes.block_length()
            10
        """
        return self._key_size

    def decrypt(self, C, K):
        r"""
        Return an 8-bit plaintext corresponding to the ciphertext ``C``,
        using S-DES decryption with key ``K``. The decryption process of
        S-DES is as follows. Let `P` be the initial permutation function,
        `P^{-1}` the corresponding inverse permutation, `\Pi_F` the
        permutation/substitution function, and `\sigma` the switch function.
        The ciphertext block ``C`` first goes through `P`, the output of
        which goes through `\Pi_F` using the second subkey. Then we apply
        the switch function to the output of the last function, and the
        result is then fed into `\Pi_F` using the first subkey. Finally,
        run the output through `P^{-1}` to get the plaintext.

        INPUT:

        - ``C`` -- an 8-bit ciphertext; a block of 8 bits

        - ``K`` -- a 10-bit key; a block of 10 bits

        OUTPUT:

        The 8-bit plaintext corresponding to ``C``, obtained using the
        key ``K``.

        EXAMPLES:

        Decrypt an 8-bit ciphertext block::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: C = [0, 1, 0, 1, 0, 1, 0, 1]
            sage: K = [1, 0, 1, 0, 0, 0, 0, 0, 1, 0]
            sage: sdes.decrypt(C, K)
            [0, 0, 0, 1, 0, 1, 0, 1]

        We can also work with strings of bits::

            sage: C = "01010101"
            sage: K = "1010000010"
            sage: sdes.decrypt(sdes.string_to_list(C), sdes.string_to_list(K))
            [0, 0, 0, 1, 0, 1, 0, 1]

        TESTS:

        The ciphertext must be a block of 8 bits::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: sdes.decrypt("C", "K")
            Traceback (most recent call last):
            ...
            TypeError: ciphertext must be a list of 8 bits
            sage: sdes.decrypt([], "K")
            Traceback (most recent call last):
            ...
            ValueError: ciphertext must be a list of 8 bits
            sage: sdes.decrypt([1, 2, 3, 4], "K")
            Traceback (most recent call last):
            ...
            ValueError: ciphertext must be a list of 8 bits

        The key must be a block of 10 bits::

            sage: sdes.decrypt([1, 0, 1, 0, 1, 1, 0, 1], "K")
            Traceback (most recent call last):
            ...
            TypeError: the key must be a list of 10 bits
            sage: sdes.decrypt([1, 0, 1, 0, 1, 1, 0, 1], [])
            Traceback (most recent call last):
            ...
            TypeError: the key must be a list of 10 bits
            sage: sdes.decrypt([1, 0, 1, 0, 1, 1, 0, 1], [1, 2, 3, 4, 5])
            Traceback (most recent call last):
            ...
            TypeError: the key must be a list of 10 bits

        The value of each element of ``C`` or ``K`` must be either 0 or 1::

            sage: C = [1, 2, 3, 4, 5, 6, 7, 8]
            sage: K = [11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
            sage: sdes.decrypt(C, K)
            Traceback (most recent call last):
            ...
            TypeError: Argument x (= 2) is not a valid string.
            sage: C = [0, 1, 0, 0, 1, 1, 1, 0]
            sage: K = [11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
            sage: sdes.decrypt(C, K)
            Traceback (most recent call last):
            ...
            TypeError: Argument x (= 13) is not a valid string.
        """
        # sanity check
        if not isinstance(C, list):
            raise TypeError("ciphertext must be a list of 8 bits")
        if len(C) != 8:
            raise ValueError("ciphertext must be a list of 8 bits")
        if not isinstance(K, list):
            raise TypeError("the key must be a list of 10 bits")
        if len(K) != 10:
            raise TypeError("the key must be a list of 10 bits")

        # run through initial permutation
        P = self.initial_permutation(C, inverse=False)
        # run through Pi_F with subkey 2
        P = self.permute_substitute(P, self.subkey(K, n=2))
        # run through switch function
        P = self.switch(P)
        # run through Pi_F with subkey 1
        P = self.permute_substitute(P, self.subkey(K, n=1))
        # run through inverse permutation
        P = self.initial_permutation(P, inverse=True)
        # output the plaintext
        return P

    def encrypt(self, P, K):
        r"""
        Return an 8-bit ciphertext corresponding to the plaintext ``P``,
        using S-DES encryption with key ``K``. The encryption process of
        S-DES is as follows. Let `P` be the initial permutation function,
        `P^{-1}` the corresponding inverse permutation, `\Pi_F` the
        permutation/substitution function, and `\sigma` the switch function.
        The plaintext block ``P`` first goes through `P`, the output of
        which goes through `\Pi_F` using the first subkey. Then we apply
        the switch function to the output of the last function, and the
        result is then fed into `\Pi_F` using the second subkey. Finally,
        run the output through `P^{-1}` to get the ciphertext.

        INPUT:

        - ``P`` -- an 8-bit plaintext; a block of 8 bits

        - ``K`` -- a 10-bit key; a block of 10 bits

        OUTPUT:

        The 8-bit ciphertext corresponding to ``P``, obtained using the
        key ``K``.

        EXAMPLES:

        Encrypt an 8-bit plaintext block::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: P = [0, 1, 0, 1, 0, 1, 0, 1]
            sage: K = [1, 0, 1, 0, 0, 0, 0, 0, 1, 0]
            sage: sdes.encrypt(P, K)
            [1, 1, 0, 0, 0, 0, 0, 1]

        We can also work with strings of bits::

            sage: P = "01010101"
            sage: K = "1010000010"
            sage: sdes.encrypt(sdes.string_to_list(P), sdes.string_to_list(K))
            [1, 1, 0, 0, 0, 0, 0, 1]

        TESTS:

        The plaintext must be a block of 8 bits::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: sdes.encrypt("P", "K")
            Traceback (most recent call last):
            ...
            TypeError: plaintext must be a list of 8 bits
            sage: sdes.encrypt([], "K")
            Traceback (most recent call last):
            ...
            ValueError: plaintext must be a list of 8 bits
            sage: sdes.encrypt([1, 2, 3, 4], "K")
            Traceback (most recent call last):
            ...
            ValueError: plaintext must be a list of 8 bits

        The key must be a block of 10 bits::

            sage: sdes.encrypt([1, 0, 1, 0, 1, 1, 0, 1], "K")
            Traceback (most recent call last):
            ...
            TypeError: the key must be a list of 10 bits
            sage: sdes.encrypt([1, 0, 1, 0, 1, 1, 0, 1], [])
            Traceback (most recent call last):
            ...
            TypeError: the key must be a list of 10 bits
            sage: sdes.encrypt([1, 0, 1, 0, 1, 1, 0, 1], [1, 2, 3, 4, 5])
            Traceback (most recent call last):
            ...
            TypeError: the key must be a list of 10 bits

        The value of each element of ``P`` or ``K`` must be either 0 or 1::

            sage: P = [1, 2, 3, 4, 5, 6, 7, 8]
            sage: K = [11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
            sage: sdes.encrypt(P, K)
            Traceback (most recent call last):
            ...
            TypeError: Argument x (= 2) is not a valid string.
            sage: P = [0, 1, 0, 0, 1, 1, 1, 0]
            sage: K = [11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
            sage: sdes.encrypt(P, K)
            Traceback (most recent call last):
            ...
            TypeError: Argument x (= 13) is not a valid string.
        """
        # sanity check
        if not isinstance(P, list):
            raise TypeError("plaintext must be a list of 8 bits")
        if len(P) != 8:
            raise ValueError("plaintext must be a list of 8 bits")
        if not isinstance(K, list):
            raise TypeError("the key must be a list of 10 bits")
        if len(K) != 10:
            raise TypeError("the key must be a list of 10 bits")

        # run through initial permutation
        C = self.initial_permutation(P, inverse=False)
        # run through Pi_F with subkey 1
        C = self.permute_substitute(C, self.subkey(K, n=1))
        # run through switch function
        C = self.switch(C)
        # run through Pi_F with subkey 2
        C = self.permute_substitute(C, self.subkey(K, n=2))
        # run through inverse permutation
        C = self.initial_permutation(C, inverse=True)
        # output the ciphertext
        return C

    def initial_permutation(self, B, inverse=False):
        r"""
        Return the initial permutation of ``B``. Denote the initial
        permutation function by `P` and let `(b_0, b_1, b_2, \dots, b_7)`
        be a vector of 8 bits, where each `b_i \in \{ 0, 1 \}`. Then

        .. MATH::


            P(b_0, b_1, b_2, b_3, b_4, b_5, b_6, b_7)
            = (b_1, b_5, b_2, b_0, b_3, b_7, b_4, b_6)

        The inverse permutation is `P^{-1}`:

        .. MATH::

            P^{-1}(b_0, b_1, b_2, b_3, b_4, b_5, b_6, b_7)
            = (b_3, b_0, b_2, b_4, b_6, b_1, b_7, b_5)

        INPUT:

        - ``B`` -- list; a block of 8 bits

        - ``inverse`` -- boolean (default: ``False``); if ``True`` then use the
          inverse permutation `P^{-1}`. If ``False`` then use the initial
          permutation `P`.

        OUTPUT:

        The initial permutation of ``B`` if ``inverse=False``, or the
        inverse permutation of ``B`` if ``inverse=True``.

        EXAMPLES:

        The initial permutation of a list of 8 bits::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: B = [1, 0, 1, 1, 0, 1, 0, 0]
            sage: P = sdes.initial_permutation(B); P
            [0, 1, 1, 1, 1, 0, 0, 0]

        Recovering the original list of 8 bits from the permutation::

            sage: Pinv = sdes.initial_permutation(P, inverse=True)
            sage: Pinv; B
            [1, 0, 1, 1, 0, 1, 0, 0]
            [1, 0, 1, 1, 0, 1, 0, 0]

        We can also work with a string of bits::

            sage: S = "10110100"
            sage: L = sdes.string_to_list(S)
            sage: P = sdes.initial_permutation(L); P
            [0, 1, 1, 1, 1, 0, 0, 0]
            sage: sdes.initial_permutation(sdes.string_to_list("01111000"), inverse=True)
            [1, 0, 1, 1, 0, 1, 0, 0]

        TESTS:

        The input block must be a list::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: sdes.initial_permutation("B")
            Traceback (most recent call last):
            ...
            TypeError: input block must be a list of 8 bits
            sage: sdes.initial_permutation(())
            Traceback (most recent call last):
            ...
            TypeError: input block must be a list of 8 bits

        The input block must be a list of 8 bits::

            sage: sdes.initial_permutation([])
            Traceback (most recent call last):
            ...
            ValueError: input block must be a list of 8 bits
            sage: sdes.initial_permutation([1, 2, 3, 4, 5, 6, 7, 8, 9])
            Traceback (most recent call last):
            ...
            ValueError: input block must be a list of 8 bits

        The value of each element of the list must be either 0 or 1::

            sage: sdes.initial_permutation([1, 2, 3, 4, 5, 6, 7, 8])
            Traceback (most recent call last):
            ...
            TypeError: Argument x (= 2) is not a valid string.
        """
        # sanity check
        if not isinstance(B, list):
            raise TypeError("input block must be a list of 8 bits")
        if len(B) != 8:
            raise ValueError("input block must be a list of 8 bits")

        bin = BinaryStrings()

        # use the initial permutation P
        if not inverse:
            return [ bin(str(B[1])), bin(str(B[5])),
                     bin(str(B[2])), bin(str(B[0])),
                     bin(str(B[3])), bin(str(B[7])),
                     bin(str(B[4])), bin(str(B[6])) ]

        # use the inverse permutation P^-1
        if inverse:
            return [ bin(str(B[3])), bin(str(B[0])),
                     bin(str(B[2])), bin(str(B[4])),
                     bin(str(B[6])), bin(str(B[1])),
                     bin(str(B[7])), bin(str(B[5])) ]

    def left_shift(self, B, n=1):
        r"""
        Return a circular left shift of ``B`` by ``n`` positions. Let
        `B = (b_0, b_1, b_2, b_3, b_4, b_5, b_6, b_7, b_8, b_9)` be a vector
        of 10 bits. Then the left shift operation `L_n` is performed on the
        first 5 bits and the last 5 bits of `B` separately. That is, if the
        number of shift positions is ``n=1``, then `L_1` is defined as

        .. MATH::

            L_1(b_0, b_1, b_2, b_3, b_4, b_5, b_6, b_7, b_8, b_9)
            = (b_1, b_2, b_3, b_4, b_0, b_6, b_7, b_8, b_9, b_5)

        If the number of shift positions is ``n=2``, then `L_2` is given by

        .. MATH::

            L_2(b_0, b_1, b_2, b_3, b_4, b_5, b_6, b_7, b_8, b_9)
            = (b_2, b_3, b_4, b_0, b_1, b_7, b_8, b_9, b_5, b_6)

        INPUT:

        - ``B`` -- list of 10 bits

        - ``n`` -- (default: 1) if ``n=1`` then perform left shift by 1
          position; if ``n=2`` then perform left shift by 2 positions. The
          valid values for ``n`` are 1 and 2, since only up to 2 positions
          are defined for this circular left shift operation.

        OUTPUT: the circular left shift of each half of ``B``

        EXAMPLES:

        Circular left shift by 1 position of a 10-bit string::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: B = [1, 0, 0, 0, 0, 0, 1, 1, 0, 0]
            sage: sdes.left_shift(B)
            [0, 0, 0, 0, 1, 1, 1, 0, 0, 0]
            sage: sdes.left_shift([1, 0, 1, 0, 0, 0, 0, 0, 1, 0])
            [0, 1, 0, 0, 1, 0, 0, 1, 0, 0]

        Circular left shift by 2 positions of a 10-bit string::

            sage: B = [0, 0, 0, 0, 1, 1, 1, 0, 0, 0]
            sage: sdes.left_shift(B, n=2)
            [0, 0, 1, 0, 0, 0, 0, 0, 1, 1]

        Here we work with a string of bits::

            sage: S = "1000001100"
            sage: L = sdes.string_to_list(S)
            sage: sdes.left_shift(L)
            [0, 0, 0, 0, 1, 1, 1, 0, 0, 0]
            sage: sdes.left_shift(sdes.string_to_list("1010000010"), n=2)
            [1, 0, 0, 1, 0, 0, 1, 0, 0, 0]

        TESTS:

        The input block must be a list::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: sdes.left_shift("B")
            Traceback (most recent call last):
            ...
            TypeError: input block must be a list of 10 bits
            sage: sdes.left_shift(())
            Traceback (most recent call last):
            ...
            TypeError: input block must be a list of 10 bits

        The input block must be a list of 10 bits::

            sage: sdes.left_shift([])
            Traceback (most recent call last):
            ...
            ValueError: input block must be a list of 10 bits
            sage: sdes.left_shift([1, 2, 3, 4, 5])
            Traceback (most recent call last):
            ...
            ValueError: input block must be a list of 10 bits

        The value of each element of the list must be either 0 or 1::

            sage: sdes.left_shift([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
            Traceback (most recent call last):
            ...
            TypeError: Argument x (= 2) is not a valid string.

        The number of shift positions must be either 1 or 2::

            sage: B = [0, 0, 0, 0, 1, 1, 1, 0, 0, 0]
            sage: sdes.left_shift(B, n=-1)
            Traceback (most recent call last):
            ...
            ValueError: input n must be either 1 or 2
            sage: sdes.left_shift(B, n=3)
            Traceback (most recent call last):
            ...
            ValueError: input n must be either 1 or 2
        """
        # sanity check
        if not isinstance(B, list):
            raise TypeError("input block must be a list of 10 bits")
        if len(B) != 10:
            raise ValueError("input block must be a list of 10 bits")

        bin = BinaryStrings()
        # circular left shift by 1 position
        if n == 1:
            return [ bin(str(B[1])), bin(str(B[2])),
                     bin(str(B[3])), bin(str(B[4])),
                     bin(str(B[0])), bin(str(B[6])),
                     bin(str(B[7])), bin(str(B[8])),
                     bin(str(B[9])), bin(str(B[5])) ]
        # circular left shift by 2 positions
        elif n == 2:
            return [ bin(str(B[2])), bin(str(B[3])),
                     bin(str(B[4])), bin(str(B[0])),
                     bin(str(B[1])), bin(str(B[7])),
                     bin(str(B[8])), bin(str(B[9])),
                     bin(str(B[5])), bin(str(B[6])) ]
        # an invalid number of shift positions
        else:
            raise ValueError("input n must be either 1 or 2")

    def list_to_string(self, B):
        r"""
        Return a binary string representation of the list ``B``.

        INPUT:

        - ``B`` -- a non-empty list of bits

        OUTPUT: the binary string representation of ``B``

        EXAMPLES:

        A binary string representation of a list of bits::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: L = [0, 0, 0, 0, 1, 1, 0, 1, 0, 0]
            sage: sdes.list_to_string(L)
            0000110100

        TESTS:

        Input ``B`` must be a non-empty list::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: sdes.list_to_string("L")
            Traceback (most recent call last):
            ...
            TypeError: input B must be a non-empty list of bits
            sage: sdes.list_to_string([])
            Traceback (most recent call last):
            ...
            ValueError: input B must be a non-empty list of bits

        Input must be a non-empty list of bits::

            sage: sdes.list_to_string([0, 1, 2])
            <repr(<sage.monoids.string_monoid_element.StringMonoidElement at 0x...>) failed: IndexError: tuple index out of range>
        """
        # sanity check
        if not isinstance(B, list):
            raise TypeError("input B must be a non-empty list of bits")
        if len(B) == 0:
            raise ValueError("input B must be a non-empty list of bits")

        # perform the conversion from list to binary string
        from sage.rings.integer import Integer
        bin = BinaryStrings()
        return bin([Integer(str(b)) for b in B])

    def permutation4(self, B):
        r"""
        Return a permutation of a 4-bit string. This permutation is called
        `P_4` and is specified as follows. Let
        `(b_0, b_1, b_2, b_3)` be a vector of 4 bits where each
        `b_i \in \{ 0, 1 \}`. Then `P_4` is defined by

        .. MATH::

            P_4(b_0, b_1, b_2, b_3) = (b_1, b_3, b_2, b_0)

        INPUT:

        - ``B`` -- a block of 4-bit string

        OUTPUT: a permutation of ``B``

        EXAMPLES:

        Permute a 4-bit string::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: B = [1, 1, 0, 0]
            sage: sdes.permutation4(B)
            [1, 0, 0, 1]
            sage: sdes.permutation4([0, 1, 0, 1])
            [1, 1, 0, 0]

        We can also work with a string of bits::

            sage: S = "1100"
            sage: L = sdes.string_to_list(S)
            sage: sdes.permutation4(L)
            [1, 0, 0, 1]
            sage: sdes.permutation4(sdes.string_to_list("0101"))
            [1, 1, 0, 0]

        TESTS:

        The input block must be a list::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: sdes.permutation4("B")
            Traceback (most recent call last):
            ...
            TypeError: input block must be a list of 4 bits
            sage: sdes.permutation4(())
            Traceback (most recent call last):
            ...
            TypeError: input block must be a list of 4 bits

        The input block must be a list of 4 bits::

            sage: sdes.permutation4([])
            Traceback (most recent call last):
            ...
            ValueError: input block must be a list of 4 bits
            sage: sdes.permutation4([1, 2, 3, 4, 5])
            Traceback (most recent call last):
            ...
            ValueError: input block must be a list of 4 bits

        The value of each element of the list must be either 0 or 1::

            sage: sdes.permutation4([1, 2, 3, 4])
            Traceback (most recent call last):
            ...
            TypeError: Argument x (= 2) is not a valid string.
        """
        # sanity check
        if not isinstance(B, list):
            raise TypeError("input block must be a list of 4 bits")
        if len(B) != 4:
            raise ValueError("input block must be a list of 4 bits")

        # perform the permutation
        bin = BinaryStrings()
        return [ bin(str(B[1])), bin(str(B[3])),
                 bin(str(B[2])), bin(str(B[0])) ]

    def permutation8(self, B):
        r"""
        Return a permutation of an 8-bit string. This permutation is called
        `P_8` and is specified as follows. Let
        `(b_0, b_1, b_2, b_3, b_4, b_5, b_6, b_7, b_8, b_9)` be a vector of
        10 bits where each `b_i \in \{ 0, 1 \}`. Then `P_8` picks out 8 of
        those 10 bits and permutes those 8 bits:

        .. MATH::

            P_8(b_0, b_1, b_2, b_3, b_4, b_5, b_6, b_7, b_8, b_9)
            =
            (b_5, b_2, b_6, b_3, b_7, b_4, b_9, b_8)

        INPUT:

        - ``B`` -- a block of 10-bit string

        OUTPUT: pick out 8 of the 10 bits of ``B`` and permute those 8 bits

        EXAMPLES:

        Permute a 10-bit string::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: B = [1, 1, 0, 0, 1, 0, 0, 1, 0, 1]
            sage: sdes.permutation8(B)
            [0, 0, 0, 0, 1, 1, 1, 0]
            sage: sdes.permutation8([0, 1, 1, 0, 1, 0, 0, 1, 0, 1])
            [0, 1, 0, 0, 1, 1, 1, 0]
            sage: sdes.permutation8([0, 0, 0, 0, 1, 1, 1, 0, 0, 0])
            [1, 0, 1, 0, 0, 1, 0, 0]

        We can also work with a string of bits::

            sage: S = "1100100101"
            sage: L = sdes.string_to_list(S)
            sage: sdes.permutation8(L)
            [0, 0, 0, 0, 1, 1, 1, 0]
            sage: sdes.permutation8(sdes.string_to_list("0110100101"))
            [0, 1, 0, 0, 1, 1, 1, 0]

        TESTS:

        The input block must be a list::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: sdes.permutation8("B")
            Traceback (most recent call last):
            ...
            TypeError: input block must be a list of 10 bits
            sage: sdes.permutation8(())
            Traceback (most recent call last):
            ...
            TypeError: input block must be a list of 10 bits

        The input block must be a list of 10 bits::

            sage: sdes.permutation8([])
            Traceback (most recent call last):
            ...
            ValueError: input block must be a list of 10 bits
            sage: sdes.permutation8([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
            Traceback (most recent call last):
            ...
            ValueError: input block must be a list of 10 bits

        The value of each element of the list must be either 0 or 1::

            sage: sdes.permutation8([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
            Traceback (most recent call last):
            ...
            TypeError: Argument x (= 6) is not a valid string.
        """
        # sanity check
        if not isinstance(B, list):
            raise TypeError("input block must be a list of 10 bits")
        if len(B) != 10:
            raise ValueError("input block must be a list of 10 bits")

        # perform the permutation
        bin = BinaryStrings()
        return [ bin(str(B[5])), bin(str(B[2])),
                 bin(str(B[6])), bin(str(B[3])),
                 bin(str(B[7])), bin(str(B[4])),
                 bin(str(B[9])), bin(str(B[8])) ]

    def permutation10(self, B):
        r"""
        Return a permutation of a 10-bit string. This permutation is called
        `P_{10}` and is specified as follows. Let
        `(b_0, b_1, b_2, b_3, b_4, b_5, b_6, b_7, b_8, b_9)` be a vector of
        10 bits where each `b_i \in \{ 0, 1 \}`. Then `P_{10}` is given by

        .. MATH::

            P_{10}(b_0, b_1, b_2, b_3, b_4, b_5, b_6, b_7, b_8, b_9)
            =
            (b_2, b_4, b_1, b_6, b_3, b_9, b_0, b_8, b_7, b_5)

        INPUT:

        - ``B`` -- a block of 10-bit string

        OUTPUT: a permutation of ``B``

        EXAMPLES:

        Permute a 10-bit string::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: B = [1, 1, 0, 0, 1, 0, 0, 1, 0, 1]
            sage: sdes.permutation10(B)
            [0, 1, 1, 0, 0, 1, 1, 0, 1, 0]
            sage: sdes.permutation10([0, 1, 1, 0, 1, 0, 0, 1, 0, 1])
            [1, 1, 1, 0, 0, 1, 0, 0, 1, 0]
            sage: sdes.permutation10([1, 0, 1, 0, 0, 0, 0, 0, 1, 0])
            [1, 0, 0, 0, 0, 0, 1, 1, 0, 0]

        Here we work with a string of bits::

            sage: S = "1100100101"
            sage: L = sdes.string_to_list(S)
            sage: sdes.permutation10(L)
            [0, 1, 1, 0, 0, 1, 1, 0, 1, 0]
            sage: sdes.permutation10(sdes.string_to_list("0110100101"))
            [1, 1, 1, 0, 0, 1, 0, 0, 1, 0]

        TESTS:

        The input block must be a list::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: sdes.permutation10("B")
            Traceback (most recent call last):
            ...
            TypeError: input block must be a list of 10 bits
            sage: sdes.permutation10(())
            Traceback (most recent call last):
            ...
            TypeError: input block must be a list of 10 bits

        The input block must be a list of 10 bits::

            sage: sdes.permutation10([])
            Traceback (most recent call last):
            ...
            ValueError: input block must be a list of 10 bits
            sage: sdes.permutation10([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
            Traceback (most recent call last):
            ...
            ValueError: input block must be a list of 10 bits

        The value of each element of the list must be either 0 or 1::

            sage: sdes.permutation10([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
            Traceback (most recent call last):
            ...
            TypeError: Argument x (= 3) is not a valid string.
        """
        # sanity check
        if not isinstance(B, list):
            raise TypeError("input block must be a list of 10 bits")
        if len(B) != 10:
            raise ValueError("input block must be a list of 10 bits")

        # perform the permutation
        bin = BinaryStrings()
        return [ bin(str(B[2])), bin(str(B[4])),
                 bin(str(B[1])), bin(str(B[6])),
                 bin(str(B[3])), bin(str(B[9])),
                 bin(str(B[0])), bin(str(B[8])),
                 bin(str(B[7])), bin(str(B[5])) ]

    def permute_substitute(self, B, key):
        r"""
        Apply the function `\Pi_F` on the block ``B`` using subkey ``key``.
        Let `(b_0, b_1, b_2, b_3, b_4, b_5, b_6, b_7)`
        be a vector of 8 bits where each `b_i \in \{ 0, 1 \}`, let `L` and
        `R` be the leftmost 4 bits and rightmost 4 bits of ``B``
        respectively, and let `F` be a function mapping 4-bit strings to
        4-bit strings. Then

        .. MATH::

            \Pi_F(L, R) = (L \oplus F(R, S), R)

        where `S` is a subkey and `\oplus` denotes the bit-wise
        exclusive-OR function.

        The function `F` can be described as follows. Its 4-bit input block
        `(n_0, n_1, n_2, n_3)` is first expanded into an 8-bit block
        to become `(n_3, n_0, n_1, n_2, n_1, n_2, n_3, n_0)`. This is
        usually represented as follows

        .. MATH::

            \begin{tabular}{c|cc|c}
              $n_3$ & $n_0$ & $n_1$ & $n_2$ \\
              $n_1$ & $n_2$ & $n_3$ & $n_0$
            \end{tabular}

        Let `K = (k_0, k_1, k_2, k_3, k_4, k_5, k_6, k_7)` be an 8-bit
        subkey. Then `K` is added to the above expanded input block using
        exclusive-OR to produce

        .. MATH::

            \begin{tabular}{c|cc|c}
              $n_3 + k_0$ & $n_0 + k_1$ & $n_1 + k_2$ & $n_2 + k_3$ \\
              $n_1 + k_4$ & $n_2 + k_5$ & $n_3 + k_6$ & $n_0 + k_7$
            \end{tabular}
            =
            \begin{tabular}{c|cc|c}
              $p_{0,0}$ & $p_{0,1}$ & $p_{0,2}$ & $p_{0,3}$ \\
              $p_{1,0}$ & $p_{1,1}$ & $p_{1,2}$ & $p_{1,3}$
            \end{tabular}

        Now read the first row as the 4-bit string
        `p_{0,0} p_{0,3} p_{0,1} p_{0,2}` and input this 4-bit string through
        S-box `S_0` to get a 2-bit output.

        .. MATH::

            S_0
            =
            \begin{tabular}{cc|cc}            \hline
              Input & Output & Input & Output \\\hline
              0000  & 01     & 1000  & 00     \\
              0001  & 00     & 1001  & 10     \\
              0010  & 11     & 1010  & 01     \\
              0011  & 10     & 1011  & 11     \\
              0100  & 11     & 1100  & 11     \\
              0101  & 10     & 1101  & 01     \\
              0110  & 01     & 1110  & 11     \\
              0111  & 00     & 1111  & 10     \\\hline
            \end{tabular}

        Next read the second row as the 4-bit string
        `p_{1,0} p_{1,3} p_{1,1} p_{1,2}` and input this 4-bit string through
        S-box `S_1` to get another 2-bit output.

        .. MATH::

            S_1
            =
            \begin{tabular}{cc|cc}            \hline
              Input & Output & Input & Output \\\hline
              0000  & 00     & 1000  & 11     \\
              0001  & 01     & 1001  & 00     \\
              0010  & 10     & 1010  & 01     \\
              0011  & 11     & 1011  & 00     \\
              0100  & 10     & 1100  & 10     \\
              0101  & 00     & 1101  & 01     \\
              0110  & 01     & 1110  & 00     \\
              0111  & 11     & 1111  & 11     \\\hline
            \end{tabular}

        Denote the 4 bits produced by `S_0` and `S_1` as `b_0 b_1 b_2 b_3`.
        This 4-bit string undergoes another permutation called `P_4` as
        follows:

        .. MATH::

            P_4(b_0, b_1, b_2, b_3) = (b_1, b_3, b_2, b_0)

        The output of `P_4` is the output of the function `F`.

        INPUT:

        - ``B`` -- list of 8 bits

        - ``key`` -- an 8-bit subkey

        OUTPUT: the result of applying the function `\Pi_F` to ``B``

        EXAMPLES:

        Applying the function `\Pi_F` to an 8-bit block and an 8-bit subkey::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: B = [1, 0, 1, 1, 1, 1, 0, 1]
            sage: K = [1, 1, 0, 1, 0, 1, 0, 1]
            sage: sdes.permute_substitute(B, K)
            [1, 0, 1, 0, 1, 1, 0, 1]

        We can also work with strings of bits::

            sage: B = "10111101"
            sage: K = "11010101"
            sage: B = sdes.string_to_list(B); K = sdes.string_to_list(K)
            sage: sdes.permute_substitute(B, K)
            [1, 0, 1, 0, 1, 1, 0, 1]

        TESTS:

        The input ``B`` must be a block of 8 bits::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: sdes.permute_substitute("B", "K")
            Traceback (most recent call last):
            ...
            TypeError: input B must be an 8-bit string
            sage: sdes.permute_substitute([], "K")
            Traceback (most recent call last):
            ...
            ValueError: input B must be an 8-bit string

        The input ``key`` must be an 8-bit subkey::

            sage: sdes.permute_substitute([0, 1, 0, 0, 1, 1, 1, 0], "K")
            Traceback (most recent call last):
            ...
            TypeError: input key must be an 8-bit subkey
            sage: sdes.permute_substitute([0, 1, 0, 0, 1, 1, 1, 0], [])
            Traceback (most recent call last):
            ...
            ValueError: input key must be an 8-bit subkey

        The value of each element of ``B`` or ``key`` must be either 0 or 1::

            sage: B = [1, 2, 3, 4, 5, 6, 7, 8]
            sage: K = [0, 1, 2, 3, 4, 5, 6, 7]
            sage: sdes.permute_substitute(B, K)
            Traceback (most recent call last):
            ...
            TypeError: Argument x (= 2) is not a valid string.
            sage: B = [0, 1, 0, 0, 1, 1, 1, 0]
            sage: K = [1, 2, 3, 4, 5, 6, 7, 8]
            sage: sdes.permute_substitute(B, K)
            Traceback (most recent call last):
            ...
            TypeError: Argument x (= 2) is not a valid string.
        """
        # sanity check
        if not isinstance(B, list):
            raise TypeError("input B must be an 8-bit string")
        if len(B) != 8:
            raise ValueError("input B must be an 8-bit string")
        if not isinstance(key, list):
            raise TypeError("input key must be an 8-bit subkey")
        if len(key) != 8:
            raise ValueError("input key must be an 8-bit subkey")

        from sage.rings.finite_rings.finite_field_constructor import FiniteField
        GF = FiniteField(2, "x")
        bin = BinaryStrings()
        bin_to_GF2 = {bin("0"): GF(0), bin("1"): GF(1)}

        # the leftmost 4 bits of B
        L = [ bin_to_GF2[bin(str(B[i]))] for i in range(4) ]
        # the rightmost 4 bits of B
        R = [ bin_to_GF2[bin(str(B[i]))] for i in range(4, len(B)) ]
        # get the GF(2) representation of the subkey
        K = [ bin_to_GF2[bin(str(key[i]))] for i in range(len(key)) ]
        # expand the rightmost 4 bits into an 8-bit block
        RX = [ R[3], R[0], R[1], R[2], R[1], R[2], R[3], R[0] ]
        # add the subkey to the expanded 8-bit block using exclusive-OR
        P = [ RX[i] + K[i] for i in range(len(K)) ]
        # run each half of P separately through the S-boxes
        left = self._sbox0([ P[0], P[3], P[1], P[2] ])
        right = self._sbox1([ P[4], P[7], P[5], P[6] ])
        # First concatenate the left and right parts, then get the
        # output of the function F.
        F = self.permutation4(left + right)
        F = [ bin_to_GF2[F[i]] for i in range(len(F)) ]
        # Add L to F using exclusive-OR. Then concatenate the result with
        # the rightmost 4 bits of B. This is the output of the function Pi_F.
        L = [ L[i] + F[i] for i in range(len(F)) ]
        return L + R

    def random_key(self):
        r"""
        Return a random 10-bit key.

        EXAMPLES:

        The size of each key is the same as the block size::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: key = sdes.random_key()
            sage: len(key) == sdes.block_length()
            True
        """
        from sage.misc.prandom import randint
        bin = BinaryStrings()
        return [bin(str(randint(0, 1))) for i in range(self._key_size)]

    def sbox(self):
        r"""
        Return the S-boxes of simplified DES.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: sbox = sdes.sbox()
            sage: sbox[0]; sbox[1]
            (1, 0, 3, 2, 3, 2, 1, 0, 0, 2, 1, 3, 3, 1, 3, 2)
            (0, 1, 2, 3, 2, 0, 1, 3, 3, 0, 1, 0, 2, 1, 0, 3)
        """
        return [self._sbox0, self._sbox1]

    def string_to_list(self, S):
        r"""
        Return a list representation of the binary string ``S``.

        INPUT:

        - ``S`` -- string of bits

        OUTPUT: list representation of the string ``S``

        EXAMPLES:

        A list representation of a string of bits::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: S = "0101010110"
            sage: sdes.string_to_list(S)
            [0, 1, 0, 1, 0, 1, 0, 1, 1, 0]

        TESTS:

        Input must be a non-empty string::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: sdes.string_to_list("")
            Traceback (most recent call last):
            ...
            ValueError: input S must be a non-empty string of bits
            sage: sdes.string_to_list(1)
            Traceback (most recent call last):
            ...
            TypeError: input S must be a non-empty string of bits

        Input must be a non-empty string of bits::

            sage: sdes.string_to_list("0123")
            Traceback (most recent call last):
            ...
            TypeError: Argument x (= 2) is not a valid string.
        """
        # sanity check
        if not isinstance(S, str):
            raise TypeError("input S must be a non-empty string of bits")
        if len(S) == 0:
            raise ValueError("input S must be a non-empty string of bits")

        # perform the conversion from string to list
        bin = BinaryStrings()
        return [bin(s) for s in S]

    def subkey(self, K, n=1):
        r"""
        Return the `n`-th subkey based on the key ``K``.

        INPUT:

        - ``K`` -- a 10-bit secret key of this Simplified DES

        - ``n`` -- (default: 1) if ``n=1`` then return the first subkey based
          on ``K``; if ``n=2`` then return the second subkey. The valid
          values for ``n`` are 1 and 2, since only two subkeys are defined
          for each secret key in Schaefer's S-DES.

        OUTPUT: the `n`-th subkey based on the secret key ``K``

        EXAMPLES:

        Obtain the first subkey from a secret key::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: key = [1, 0, 1, 0, 0, 0, 0, 0, 1, 0]
            sage: sdes.subkey(key, n=1)
            [1, 0, 1, 0, 0, 1, 0, 0]

        Obtain the second subkey from a secret key::

            sage: key = [1, 0, 1, 0, 0, 0, 0, 0, 1, 0]
            sage: sdes.subkey(key, n=2)
            [0, 1, 0, 0, 0, 0, 1, 1]

        We can also work with strings of bits::

            sage: K = "1010010010"
            sage: L = sdes.string_to_list(K)
            sage: sdes.subkey(L, n=1)
            [1, 0, 1, 0, 0, 1, 0, 1]
            sage: sdes.subkey(sdes.string_to_list("0010010011"), n=2)
            [0, 1, 1, 0, 1, 0, 1, 0]

        TESTS:

        Input ``K`` must be a 10-bit key::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: sdes.subkey("K")
            Traceback (most recent call last):
            ...
            TypeError: input K must be a 10-bit key
            sage: sdes.subkey([])
            Traceback (most recent call last):
            ...
            ValueError: input K must be a 10-bit key

        There are only two subkeys::

            sage: key = [1, 0, 1, 0, 0, 0, 0, 0, 1, 0]
            sage: sdes.subkey(key, n=0)
            Traceback (most recent call last):
            ...
            ValueError: input n must be either 1 or 2
            sage: sdes.subkey(key, n=3)
            Traceback (most recent call last):
            ...
            ValueError: input n must be either 1 or 2
        """
        # sanity check
        if not isinstance(K, list):
            raise TypeError("input K must be a 10-bit key")
        if len(K) != self._key_size:
            raise ValueError("input K must be a 10-bit key")

        # get the first subkey
        if n == 1:
            key1 = self.permutation10(K)
            key1 = self.left_shift(key1, n=1)
            return self.permutation8(key1)
        # get the second subkey
        elif n == 2:
            key2 = self.permutation10(K)
            key2 = self.left_shift(key2, n=1)
            key2 = self.left_shift(key2, n=2)
            return self.permutation8(key2)
        # an invalid subkey number
        else:
            raise ValueError("input n must be either 1 or 2")

    def switch(self, B):
        r"""
        Interchange the first 4 bits with the last 4 bits in the list ``B``
        of 8 bits. Let `(b_0, b_1, b_2, b_3, b_4, b_5, b_6, b_7)`
        be a vector of 8 bits, where each `b_i \in \{ 0, 1 \}`. Then the
        switch function `\sigma` is given by

        .. MATH::

            \sigma(b_0, b_1, b_2, b_3, b_4, b_5, b_6, b_7)
            = (b_4, b_5, b_6, b_7, b_0, b_1, b_2, b_3)

        INPUT:

        - ``B`` -- list; a block of 8 bits

        OUTPUT:

        A block of the same dimension, but in which the first 4 bits from
        ``B`` has been switched for the last 4 bits in ``B``.

        EXAMPLES:

        Interchange the first 4 bits with the last 4 bits::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: B = [1, 1, 1, 0, 1, 0, 0, 0]
            sage: sdes.switch(B)
            [1, 0, 0, 0, 1, 1, 1, 0]
            sage: sdes.switch([1, 1, 1, 1, 0, 0, 0, 0])
            [0, 0, 0, 0, 1, 1, 1, 1]

        We can also work with a string of bits::

            sage: S = "11101000"
            sage: L = sdes.string_to_list(S)
            sage: sdes.switch(L)
            [1, 0, 0, 0, 1, 1, 1, 0]
            sage: sdes.switch(sdes.string_to_list("11110000"))
            [0, 0, 0, 0, 1, 1, 1, 1]

        TESTS:

        The input block must be a list::

            sage: from sage.crypto.block_cipher.sdes import SimplifiedDES
            sage: sdes = SimplifiedDES()
            sage: sdes.switch("B")
            Traceback (most recent call last):
            ...
            TypeError: input block must be a list of 8 bits
            sage: sdes.switch(())
            Traceback (most recent call last):
            ...
            TypeError: input block must be a list of 8 bits

        The input block must be a list of 8 bits::

            sage: sdes.switch([])
            Traceback (most recent call last):
            ...
            ValueError: input block must be a list of 8 bits
            sage: sdes.switch([1, 2, 3, 4, 5, 6, 7, 8, 9])
            Traceback (most recent call last):
            ...
            ValueError: input block must be a list of 8 bits

        The value of each element of the list must be either 0 or 1::

            sage: sdes.switch([1, 2, 3, 4, 5, 6, 7, 8])
            Traceback (most recent call last):
            ...
            TypeError: Argument x (= 5) is not a valid string.
        """
        # sanity check
        if not isinstance(B, list):
            raise TypeError("input block must be a list of 8 bits")
        if len(B) != 8:
            raise ValueError("input block must be a list of 8 bits")

        # perform the switch
        bin = BinaryStrings()
        return [ bin(str(B[4])), bin(str(B[5])),
                 bin(str(B[6])), bin(str(B[7])),
                 bin(str(B[0])), bin(str(B[1])),
                 bin(str(B[2])), bin(str(B[3])) ]
