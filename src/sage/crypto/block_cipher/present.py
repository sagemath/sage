r"""
PRESENT

An ultra-lightweight block cipher.

This file implements the PRESENT block cipher and the corresponding key
schedule as described in [BKLPPRSV2007]_. PRESENT is an example of an
SP-network and consists of 31 rounds. The block length is 64 bits and two key
lengths of 80 and 128 bits are supported.

Note, this implementation is ment for experimental and educational usage only,
do not use it in production code!

EXAMPLES:

Encrypt a message::

    sage: from sage.crypto.block_cipher.present import PRESENT
    sage: present = PRESENT()
    sage: present.encrypt(plaintext=0, key=0).hex()
    '2844b365c06992a3'

And decrypt it again::

    sage: present.decrypt(ciphertext=0x2844b365c06992a3, key=0)
    0

Have a look at the used round keys::

    sage: from sage.crypto.block_cipher.present import PRESENT_KS
    sage: ks = PRESENT_KS()
    sage: [k.hex() for k in ks(0)]
    ['0',
     'c000000000000000',
      ...
     '6dab31744f41d700']

Tweak around with the cipher::

    sage: from sage.crypto.sbox import SBox
    sage: cipher = PRESENT(rounds=1, doFinalRound=False)
    sage: cipher.sbox = SBox(range(16))
    sage: cipher.keySchedule = lambda x: [0, 0]  # return the 0 keys as round keys
    sage: cipher.encrypt(plaintext=0x1234, key=0x0).hex()
    '1234'

AUTHORS:

- Lukas Stennes (2019-02-01): initial version
"""

# ****************************************************************************
#       Copyright (C) 2013 Lukas Stennes <lukas.stennes@rub.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.structure.sage_object import SageObject
from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.modules.free_module_element import vector
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.crypto.sboxes import PRESENT as PRESENTSBOX
from sage.modules.vector_mod2_dense import Vector_mod2_dense


def _smallscale_present_linearlayer(nsboxes=16):
    """
    .. TODO::

        switch to sage.crypto.linearlayer
        (:trac:`25735`) as soon as it is included in sage

    EXAMPLES::

        sage: from sage.crypto.block_cipher.present import _smallscale_present_linearlayer
        sage: _smallscale_present_linearlayer(4)
        [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0]
        [0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0]
        [0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0]
        [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]
    """
    from sage.modules.free_module import VectorSpace
    from sage.modules.free_module_element import vector
    from sage.matrix.constructor import Matrix
    from sage.rings.finite_rings.finite_field_constructor import GF

    def present_llayer(n, x):
        dim = 4*n
        y = [0]*dim
        for i in range(dim-1):
            y[i] = x[(n * i) % (dim - 1)]
        y[dim-1] = x[dim-1]
        return vector(GF(2), y)

    m = Matrix(GF(2), [present_llayer(nsboxes, ei)
                       for ei in VectorSpace(GF(2), 4*nsboxes).basis()])
    return m


class PRESENT(SageObject):
    r"""
    This class implements PRESENT described in [BKLPPRSV2007]_.

    EXAMPLES:

    You can invoke PRESENT encryption/decryption either by calling PRESENT with
    an appropriate flag::

        sage: from sage.crypto.block_cipher.present import PRESENT
        sage: present = PRESENT()
        sage: P = 0xFFFFFFFFFFFFFFFF
        sage: K = 0x0
        sage: present(present(P, K, 'encrypt'), K, 'decrypt') == P
        True

    Or by calling encryption/decryption methods directly::

        sage: C = present.encrypt(P, K)
        sage: P == present.decrypt(C, K)
        True

    The number of rounds can be reduced easily::

        sage: present = PRESENT(rounds=15)
        sage: present(present(P, K, 'encrypt'), K, 'decrypt') == P
        True

    You can use integers or a list-like bit representation for the inputs. If
    the input is an integer the output will be too. If it is list-like the
    output will be a bit vector::

        sage: P = ZZ(0).digits(2,padto=64)
        sage: K = ZZ(0).digits(2,padto=80)
        sage: list(present(present(P, K, 'encrypt'), K, 'decrypt')) == P
        True
        sage: P = ZZ(0).digits(2,padto=64)
        sage: K = 0x0
        sage: list(present(present(P, K, 'encrypt'), K, 'decrypt')) == P
        True

    The 80-bit version of PRESENT is used by default but the 128-bit version is
    also implemented::

        sage: present = PRESENT(128)
        sage: P = 0x0123456789abcdef
        sage: K = 0x00112233445566778899aabbccddeeff
        sage: present(present(P, K, 'encrypt'), K, 'decrypt') == P
        True

    .. SEEALSO::

        :class:`PRESENT_KS`
        :mod:`sage.crypto.sboxes`

    .. automethod:: __init__
    .. automethod:: __call__
    """

    def __init__(self, keySchedule=80, rounds=None, doFinalRound=False):
        r"""
        Construct an instance of PRESENT.

        INPUT:

        - ``keySchedule`` -- (default: ``80``); the key schedule that will be
          used for encryption and decryption. Use ``80`` or ``128`` as a
          shortcut for the original key schedules from [BKLPPRSV2007]_.

        - ``rounds``  -- integer (default: ``None``); the number of rounds. If
          ``None`` the number of rounds of the key schedule is used.

        - ``doFinalRound`` -- boolean (default: ``False``); flag to
          control whether the linear layer in the last round should take place
          or not. Since the last linear layer does not add any security, it
          usually does not take place in real world implementations for
          performance reasons.

        EXAMPLES:

        By default a 80-bit version with 31 rounds is created::

            sage: from sage.crypto.block_cipher.present import PRESENT
            sage: PRESENT() # indirect doctest
            PRESENT block cipher with 31 rounds, deactivated linear layer in
            last round and the following key schedule:
            Original PRESENT key schedule with 80-bit keys and 31 rounds

        The 128-bit version is also implemented::

            sage: PRESENT(128) # indirect doctest
            PRESENT block cipher with 31 rounds, deactivated linear layer in
            last round and the following key schedule:
            Original PRESENT key schedule with 128-bit keys and 31 rounds

        Reducing the number of rounds is simple. But increasing it is not
        possible::

            sage: PRESENT(keySchedule=80, rounds=23) # indirect doctest
            PRESENT block cipher with 23 rounds, deactivated linear layer in
            last round and the following key schedule:
            Original PRESENT key schedule with 80-bit keys and 31 rounds
            sage: PRESENT(80, 32) # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: number of rounds must be less or equal to the number
            of rounds of the key schedule


        By default the linear layer operation in the last round is omitted but
        of course you can enable it::

            sage: PRESENT(doFinalRound=True) # indirect doctest
            PRESENT block cipher with 31 rounds, activated linear layer in
            last round and the following key schedule:
            Original PRESENT key schedule with 80-bit keys and 31 rounds

        You can use arbitrary key schedules. Since it is the only one
        implemented here the original key schedule is used for demonstration::

            sage: from sage.crypto.block_cipher.present import PRESENT_KS
            sage: PRESENT(keySchedule=PRESENT_KS(80, 15)) # indirect doctest
            PRESENT block cipher with 15 rounds, deactivated linear layer in
            last round and the following key schedule:
            Original PRESENT key schedule with 80-bit keys and 15 rounds

        .. SEEALSO::

            :class:`PRESENT_KS`
        """
        if keySchedule == 80:
            self.keySchedule = PRESENT_KS()
            self._keysize = 80
        elif keySchedule == 128:
            self.keySchedule = PRESENT_KS(128)
            self._keysize = 128
        else:
            self.keySchedule = keySchedule
        if rounds is None:
            self._rounds = self.keySchedule._rounds
        elif rounds <= self.keySchedule._rounds:
            self._rounds = rounds
        else:
            raise ValueError('number of rounds must be less or equal to the '
                             'number of rounds of the key schedule')
        self._blocksize = 64
        self.sbox = PRESENTSBOX
        self._permutationMatrix = _smallscale_present_linearlayer()
        self._inversePermutationMatrix = self._permutationMatrix.inverse()
        self._doFinalRound = doFinalRound

    def __call__(self, block, key, algorithm='encrypt'):
        r"""
        Apply PRESENT encryption or decryption on ``block`` using ``key``.
        The flag ``algorithm`` controls what action is to be performed on
        ``block``.

        INPUT:

        - ``block`` -- integer or bit list-like; the plaintext or ciphertext

        - ``K`` -- integer or bit list-like; the key

        - ``algorithm`` -- string (default: ``'encrypt'``); a flag to signify
          whether encryption or decryption is to be applied to ``B``. The
          encryption flag is ``'encrypt'`` and the decryption flag is
          ``'decrypt'``

        OUTPUT:

        - The plaintext or ciphertext corresponding to ``block``, obtained using
          the ``key``. If ``block`` is an integer the output will be too. If
          ``block`` is list-like the output will be a bit vector.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.present import PRESENT
            sage: present = PRESENT(doFinalRound=True)
            sage: P = 0xFFFFFFFFFFFFFFFF
            sage: K = 0x0
            sage: present(P, K, 'encrypt').hex()
            'a112ffc72f68417b'
        """
        if algorithm == 'encrypt':
            return self.encrypt(block, key)
        elif algorithm == 'decrypt':
            return self.decrypt(block, key)
        else:
            raise ValueError('Algorithm must be \'encrypt\' or \'decrypt\' and'
                             ' not \'%s\'' % algorithm)

    def __eq__(self, other):
        r"""
        Compare ``self`` with ``other``.

        PRESENT objects are the same if all attributes are the same.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.present import PRESENT
            sage: PRESENT(80) == PRESENT(80) # indirect doctest
            True
            sage: PRESENT(80) == PRESENT(128) # indirect doctest
            False
            sage: PRESENT(80) == 80 # indirect doctest
            False
            sage: present = PRESENT()
            sage: present.sbox = present.sbox.inverse()
            sage: present == PRESENT() # indirect doctest
            False
        """
        if not isinstance(other, PRESENT):
            return False
        else:
            return self.__dict__ == other.__dict__

    def __repr__(self):
        r"""
        A string representation of this PRESENT.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.present import PRESENT
            sage: PRESENT() # indirect doctest
            PRESENT block cipher with 31 rounds, deactivated linear layer in
            last round and the following key schedule:
            Original PRESENT key schedule with 80-bit keys and 31 rounds
        """
        return ('PRESENT block cipher with %s rounds, %s linear layer in last '
                'round and the following key schedule:\n%s'
                % (self._rounds, 'activated' if self._doFinalRound else
                   'deactivated', self.keySchedule.__repr__()))

    def encrypt(self, plaintext, key):
        r"""
        Return the ciphertext corresponding to ``plaintext``, using PRESENT
        encryption with ``key``.

        INPUT:

        - ``plaintext`` -- integer or bit list-like; the plaintext that will be
          encrypted.

        - ``key`` -- integer or bit list-like; the key

        OUTPUT:

        - The ciphertext corresponding to ``plaintext``, obtained using the
          ``key``. If ``plaintext`` is an integer the output will be too. If
          ``plaintext`` is list-like the output will be a bit vector.

        EXAMPLES:

        The test vectors from [BKLPPRSV2007]_ are checked here::

            sage: from sage.crypto.block_cipher.present import PRESENT
            sage: present = PRESENT(doFinalRound=True)
            sage: p1 = 0x0
            sage: k1 = 0x0
            sage: c1 = 0x5579C1387B228445
            sage: present.encrypt(p1, k1) == c1
            True
            sage: p2 = 0x0
            sage: k2 = 0xFFFFFFFFFFFFFFFFFFFF
            sage: c2 = 0xE72C46C0F5945049
            sage: present.encrypt(p2, k2) == c2
            True
            sage: p3 = 0xFFFFFFFFFFFFFFFF
            sage: k3 = 0x0
            sage: c3 = 0xA112FFC72F68417B
            sage: present.encrypt(p3, k3) == c3
            True
            sage: p4 = 0xFFFFFFFFFFFFFFFF
            sage: k4 = 0xFFFFFFFFFFFFFFFFFFFF
            sage: c4 = 0x3333DCD3213210D2
            sage: present.encrypt(p4, k4) == c4
            True

        ALGORITHM:

        Description of the encryption function based on [BKLPPRSV2007]_:

        A top-level algorithmic description of PRESENT encryption::

            generateRoundKeys()
            for i = 1 to 31 do
                addRoundkey(STATE, K_i)
                sBoxLayer(STATE)
                pLayer(STATE)
            end for
            addRoundkey(STATE, K_{32})

        Each of the 31 rounds consists of an XOR operation to introduce a round
        key `K_i` for `1 \leq i \leq 32`, where `K_{32}` is used for
        post-whitening, a linear bitwise permutation and a non-linear
        substitution layer. The non-linear layer uses a single 4-bit S-box
        which is applied 16 times in parallel in each round. Each stage but
        addRoundkey is specified in its corresponding function.

        **addRoundkey:**
        Given round key `K_i = \kappa^i_{63} \dots \kappa^i_0` for `1 \leq i
        \leq 32` and current STATE `b_{63} \dots b_0`, addRoundkey consists of
        the operation for `0 \leq j \leq 63`, `b_j = b_j \oplus \kappa^i_j`.
        """
        if isinstance(plaintext, (list, tuple, Vector_mod2_dense)):
            inputType = 'vector'
        elif isinstance(plaintext, (Integer, int)):
            inputType = 'integer'
        state = convert_to_vector(plaintext, 64)
        key = convert_to_vector(key, self._keysize)
        roundKeys = self.keySchedule(key)
        for r, K in enumerate(roundKeys[:self._rounds]):
            state = self.round(state, r, K)
        state = state + roundKeys[self._rounds]
        return state if inputType == 'vector' else ZZ(list(state), 2)

    def decrypt(self, ciphertext, key):
        r"""
        Return the plaintext corresponding to the ``ciphertext``, using PRESENT decryption with ``key``.

        INPUT:

        - ``ciphertext`` -- integer or bit list-like; the ciphertext that will
          be decrypted

        - ``key`` -- integer or bit list-like; the key

        OUTPUT:

        - The plaintext corresponding to ``ciphertext``, obtained using the
          ``key``. If ``ciphertext`` is an integer the output will be too. If
          ``ciphertext`` is list-like the output will be a bit vector.

        EXAMPLES:

        The test vectors from [BKLPPRSV2007]_ are checked here::

            sage: from sage.crypto.block_cipher.present import PRESENT
            sage: present = PRESENT(doFinalRound=True)
            sage: p1 = 0x0
            sage: k1 = 0x0
            sage: c1 = 0x5579C1387B228445
            sage: present.decrypt(c1, k1) == p1
            True
            sage: p2 = 0x0
            sage: k2 = 0xFFFFFFFFFFFFFFFFFFFF
            sage: c2 = 0xE72C46C0F5945049
            sage: present.decrypt(c2, k2) == p2
            True
            sage: p3 = 0xFFFFFFFFFFFFFFFF
            sage: k3 = 0x0
            sage: c3 = 0xA112FFC72F68417B
            sage: present.decrypt(c3, k3) == p3
            True
            sage: p4 = 0xFFFFFFFFFFFFFFFF
            sage: k4 = 0xFFFFFFFFFFFFFFFFFFFF
            sage: c4 = 0x3333DCD3213210D2
            sage: present.decrypt(c4, k4) == p4
            True
       """
        if isinstance(ciphertext, (list, tuple, Vector_mod2_dense)):
            inputType = 'vector'
        elif isinstance(ciphertext, (Integer, int)):
            inputType = 'integer'
        state = convert_to_vector(ciphertext, 64)
        key = convert_to_vector(key, self._keysize)
        roundKeys = self.keySchedule(key)
        state = state + roundKeys[self._rounds]
        for r, K in enumerate(roundKeys[:self._rounds][::-1]):
            state = self.round(state, r, K, inverse=True)
        return state if inputType == 'vector' else ZZ(list(state), 2)

    def round(self, state, round_counter, round_key, inverse=False):
        """
        Apply one round of PRESENT to ``state`` and return the result.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.present import PRESENT
            sage: from sage.crypto.block_cipher.present import convert_to_vector
            sage: present = PRESENT(128)
            sage: k = convert_to_vector(0x0011223344556677, 64)
            sage: p = convert_to_vector(0x0123456789abcdef, 64)
            sage: ZZ(list(present.round(p, 0, k)), 2).hex()
            'ad0ed4ca386b6559'
        """
        out = state[:]
        if not inverse:
            out = out + round_key
            out = self.sbox_layer(out)
            if self._doFinalRound or round_counter != self._rounds - 1:
                out = self.linear_layer(out)
        else:
            if self._doFinalRound or round_counter != 0:
                out = self.linear_layer(out, inverse=True)
            out = self.sbox_layer(out, inverse=True)
            out = out + round_key
        return out

    def sbox_layer(self, state, inverse=False):
        r"""
        Apply the sBoxLayer of PRESENT to the bit vector ``state`` and return
        the result.

        The S-box used in PRESENT is a 4-bit to 4-bit S-box. The action of this
        box in hexadecimal notation is given by the following table.

        +------+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
        | x    | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | A | B | C | D | E | F |
        +------+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
        | S[x] | C | 5 | 6 | B | 9 | 0 | A | D | 3 | E | F | 8 | 4 | 7 | 1 | 2 |
        +------+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+

        For sBoxLayer the current STATE `b_{63} \dots b_0` is considered as
        sixteen 4-bit words `w_{15} \dots w_0` where `w_i =
        b_{4i+3}||b_{4i+2}||b_{4i+1}||b_{4i}` for `0 \leq i \leq 15` and the
        output nibble S[`w_i`] provides the updated state values in the obvious
        way.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.present import PRESENT
            sage: present = PRESENT()
            sage: state = vector(GF(2), 64, [0]*64)
            sage: present.sbox_layer(state)
            (0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0,
            1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1,
            0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1)
            sage: state = vector(GF(2), 64, [1]+[0]*63)
            sage: present.sbox_layer(state)
            (1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0,
            1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1,
            0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1)

        .. NOTE::

            :mod:`sage.crypto.sbox` uses big endian by default whereas most of
            Sage uses little endian. So to use the big endian PRESENT Sbox from
            :mod:`sage.crypto.sboxes` :func:`sbox_layer` has to do some endian
            conversion (i.e. reverse input and ouput of the Sbox). Keep this in
            mind if you change the Sbox or :func:`sbox_layer`.
        """
        sbox = self.sbox if not inverse else self.sbox.inverse()
        out = vector(GF(2), 64)
        for nibble in [slice(4*j, 4*j+4) for j in range(16)]:
            out[nibble] = sbox(state[nibble][::-1])[::-1]
        return out

    def linear_layer(self, state, inverse=False):
        """
        Apply the pLayer of PRESENT to the bit vector ``state`` and return the
        result.

        The bit permutation used in PRESENT is given by the following
        table. Bit `i` of STATE is moved to bit position `P(i)`.

        +------+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+
        | i    |  0 |  1 |  2 |  3 |  4 |  5 |  6 |  7 |  8 |  9 | 10 | 11 | 12 | 13 | 14 | 15 |
        +------+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+
        | P(i) |  0 | 16 | 32 | 48 |  1 | 17 | 33 | 49 |  2 | 18 | 34 | 50 |  3 | 19 | 35 | 51 |
        +------+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+
        +------+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+
        | i    | 16 | 17 | 18 | 19 | 20 | 21 | 22 | 23 | 24 | 25 | 26 | 27 | 28 | 29 | 30 | 31 |
        +------+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+
        | P(i) |  4 | 20 | 36 | 52 |  5 | 21 | 37 | 53 |  6 | 22 | 38 | 54 |  7 | 23 | 39 | 55 |
        +------+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+
        +------+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+
        | i    | 32 | 33 | 34 | 35 | 36 | 37 | 38 | 39 | 40 | 41 | 42 | 43 | 44 | 45 | 46 | 47 |
        +------+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+
        | P(i) |  8 | 24 | 40 | 56 |  9 | 25 | 41 | 57 | 10 | 26 | 42 | 58 | 11 | 27 | 43 | 59 |
        +------+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+
        +------+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+
        | i    | 48 | 49 | 50 | 51 | 52 | 53 | 54 | 55 | 56 | 57 | 58 | 59 | 60 | 61 | 62 | 63 |
        +------+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+
        | P(i) | 12 | 28 | 44 | 60 | 13 | 29 | 45 | 61 | 14 | 30 | 46 | 62 | 15 | 31 | 47 | 63 |
        +------+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+----+

        EXAMPLES::

            sage: from sage.crypto.block_cipher.present import PRESENT
            sage: present = PRESENT()
            sage: state = vector(GF(2), 64, [0, 1]+[0]*62)
            sage: present.linear_layer(state)
            (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        """
        m = self._permutationMatrix if not inverse else self._inversePermutationMatrix
        return m * state


class PRESENT_KS(SageObject):
    r"""
    This class implements the PRESENT key schedules for both 80-bit and 128-bit
    keys as described in [BKLPPRSV2007]_.

    EXAMPLES:

    Initialise the key schedule with a `master\_key` to use it as an iterable::

        sage: from sage.crypto.block_cipher.present import PRESENT_KS
        sage: ks = PRESENT_KS(master_key=0)
        sage: ks[0] == 0x0
        True
        sage: ks[31] == 0x6dab31744f41d700
        True

    Or omit the `master\_key` and pass a key when calling the key schedule::

        sage: ks = PRESENT_KS(keysize=128)
        sage: K = ks(0x00112233445566778899aabbccddeeff)
        sage: K[0] == 0x0011223344556677
        True
        sage: K[31] == 0x091989a5ae8eab21
        True

    ALGORITHM:

    Description of the key schedule for 64-bit and 128-bit keys from
    [BKLPPRSV2007]_:

    The key schedule for 64-bit keys works as follows:

    At round `i` the 64-bit round key `K_i = \kappa_{63}\kappa_{62} \dots
    \kappa_0` consists of the 64 leftmost bits of the current contents of
    register `K`. Thus at round `i` we have that:

    .. MATH::

        K_i = \kappa_{63}\kappa_{62}\dots \kappa_0 = k_{79}k_{78}\dots k_{16}.

    After extracting the round key `K_i`, the key register `K = k_{79}k_{78}
    \dots k_0` is updated as follows:

    .. MATH::

        \begin{aligned}
        \ [k_{79}k_{78}\dots k_{1}k_{0}] &= [k_{18}k_{17}\dots k_{20}k_{19}] \\
        [k_{79}k_{78}k_{77}k_{76}] &= S[k_{79}k_{78}k_{77}k_{76}] \\
        [k_{19}k_{18}k_{17}k_{16}k_{15}] &= [k_{19}k_{18}k_{17}k_{16}k_{15}]
                                            \oplus round\_counter
        \end{aligned}

    Thus, the key register is rotated by 61 bit positions to the left, the
    left-most four bits are passed through the PRESENT S-box, and the
    round_counter value `i` is exclusive-ored with bits `k_{19}k_{18}k_{17}
    k_{16}k_{15}` of `K` with the least significant bit of round_counter on
    the right.

    The key schedule for 128-bit keys works as follows:

    At round `i` the 64-bit round key `K_i = \kappa_{63}\kappa_{62} \dots
    \kappa_0` consists of the 64 leftmost bits fo the current contents of
    register `K`. Thus at round `i` we have that:

    .. MATH::

        K_i = \kappa_{63}\kappa_{62}\dots \kappa_0 = k_{127}k_{126}\dots k_{64}

    After extracting the round key `K_i`, the key register `K = k_{127}k_{126}
    \dots k_0` is updated as follows:

    .. MATH::

        \begin{aligned}
        \ [k_{127}k_{126}\dots k_{1}k_{0}] &= [k_{66}k_{65}\dots k_{68}k_{67}]
        \\
        [k_{127}k_{126}k_{125}k_{124}] &= S[k_{127}k_{126}k_{125}k_{124}]\\
        [k_{123}k_{122}k_{121}k_{120}] &= S[k_{123}k_{122}k_{121}k_{120}]\\
        [k_{66}k_{65}k_{64}k_{63}k_{62}] &= [k_{66}k_{65}k_{64}k_{63}k_{62}]
                                            \oplus round\_counter
        \end{aligned}

    Thus, the key register is rotated by 61 bit positions to the left, the
    left-most eight bits are passed through two PRESENT S-boxes, and the
    round_counter value `i` is exclusive-ored with bits `k_{66}k_{65}k_{64}
    k_{63}k_{62}` of `K` with the least significant bit of round_counter on
    the right.

    .. SEEALSO::

        :class:`PRESENT`
        :mod:`sage.crypto.sboxes`

    .. NOTE::

        :mod:`sage.crypto.sbox` uses big endian by default whereas most of Sage
        uses little endian. So to use the big endian PRESENT Sbox from
        :mod:`sage.crypto.sboxes` :class:`PRESENT_KS` has to do some endian
        conversion (i.e. reverse input and ouput of the Sbox). Keep this in
        mind if you change the Sbox or :func:`__call__`.

    .. automethod:: __init__
    .. automethod:: __call__
    """

    def __init__(self, keysize=80, rounds=31, master_key=None):
        r"""
        Construct an instance of PRESENT_KS.

        INPUT:

        - ``keysize`` -- integer (default: ``80``); the size of the keys that
          will be used in bits. It must be either 80 or 128.

        - ``rounds`` -- integer (default: ``31``); the number of rounds
          ``self`` can create keys for

        - ``master_key`` -- integer or bit list-like (default: ``None``); the
          key that will be used

        EXAMPLES::

            sage: from sage.crypto.block_cipher.present import PRESENT_KS
            sage: PRESENT_KS()
            Original PRESENT key schedule with 80-bit keys and 31 rounds

        .. NOTE::

            If you want to use a PRESENT_KS object as an iterable you have to
            pass a ``master_key`` value on initialisation. Otherwise you can
            omit ``master_key`` and pass a key when you call the object.
        """
        if keysize != 80 and keysize != 128:
            raise ValueError('keysize must bei either 80 or 128 and not %s'
                             % keysize)
        self._keysize = keysize
        self._rounds = rounds
        self.sbox = PRESENTSBOX
        self._master_key = master_key

    def __call__(self, K):
        r"""
        Return all round keys in a list.

        INPUT:

        - ``K`` -- integer or bit list-like; the key

        OUTPUT:

        - A list containing ``rounds + 1`` round keys. Since addRoundkey takes
          place in every round and after the last round there must be
          ``rounds + 1`` round keys.  If ``K`` is an integer the elements of
          the output list will be too. If ``K`` is list-like the element of the
          output list will be  bit vectors.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.present import PRESENT_KS
            sage: ks = PRESENT_KS()
            sage: ks(0)[31] == 0x6dab31744f41d700
            True

        .. NOTE::

            If you want to use a PRESENT_KS object as an iterable you have to
            pass a ``master_key`` value on initialisation. Otherwise you can
            omit ``master_key`` and pass a key when you call the object.
        """
        if isinstance(K, (list, tuple, Vector_mod2_dense)):
            inputType = 'vector'
        elif isinstance(K, (Integer, int)):
            inputType = 'integer'
        K = convert_to_vector(K, self._keysize)
        roundKeys = []
        if self._keysize == 80:
            for i in range(1, self._rounds+1):
                roundKeys.append(K[16:])
                K[0:] = list(K[19:]) + list(K[:19])
                K[76:] = self.sbox(K[76:][::-1])[::-1]
                rc = vector(GF(2), ZZ(i).digits(2, padto=5))
                K[15:20] = K[15:20] + rc
            roundKeys.append(K[16:])
        elif self._keysize == 128:
            for i in range(1, self._rounds+1):
                roundKeys.append(K[64:])
                K[0:] = list(K[67:]) + list(K[:67])
                K[124:] = self.sbox(K[124:][::-1])[::-1]
                K[120:124] = self.sbox(K[120:124][::-1])[::-1]
                rc = vector(GF(2), ZZ(i).digits(2, padto=5))
                K[62:67] = K[62:67] + rc
            roundKeys.append(K[64:])
        return roundKeys if inputType == 'vector' else [ZZ(list(k), 2) for k in
                                                        roundKeys]

    def __eq__(self, other):
        r"""
        Compare ``self`` with ``other``.

        PRESENT_KS objects are the same if all attributes are the same.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.present import PRESENT_KS
            sage: PRESENT_KS(80) == PRESENT_KS(80) # indirect doctest
            True
            sage: PRESENT_KS(80) == PRESENT_KS(128) # indirect doctest
            False
            sage: PRESENT_KS(80) == 80 # indirect doctest
            False
        """
        if not isinstance(other, PRESENT_KS):
            return False
        else:
            return self.__dict__ == other.__dict__

    def __repr__(self):
        r"""
        A string representation of this PRESENT_KS.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.present import PRESENT_KS
            sage: PRESENT_KS() # indirect doctest
            Original PRESENT key schedule with 80-bit keys and 31 rounds
        """
        return ('Original PRESENT key schedule with %s-bit keys and %s rounds'
                % (self._keysize, self._rounds))

    def __getitem__(self, r):
        r"""
        Computes the sub key for round ``r`` derived from initial master key.

        The key schedule object has to have been initialised with the
        `master_key` argument.

        INPUT:

        - ``r`` integer; the round for which the sub key is computed

        EXAMPLES::

            sage: from sage.crypto.block_cipher.present import PRESENT_KS
            sage: ks = PRESENT_KS(master_key=0x0)
            sage: ks[0] ==  0x0 # indirect doctest
            True
            sage: ks[31] ==  0x6dab31744f41d700 # indirect doctest
            True
        """
        if self._master_key is None:
            raise ValueError('Key not set during initialisation')
        return self(self._master_key)[r]

    def __iter__(self):
        """
        Iterate over the ``self._rounds + 1`` PRESENT round keys, derived from
        `master_key`

        EXAMPLES::

            sage: from sage.crypto.block_cipher.present import PRESENT_KS
            sage: K = [k for k in PRESENT_KS(master_key=0x0)]
            sage: K[0] == 0x0 # indirect doctest
            True
            sage: K[31] == 0x6dab31744f41d700 # indirect doctest
            True
        """
        if self._master_key is None:
            raise ValueError('Key not set during initialisation')
        return iter(self(self._master_key))


def convert_to_vector(I, L):
    r"""
    Convert ``I`` to a bit vector of length ``L``.

    INPUT:

    - ``I`` -- integer or bit list-like

    - ``L`` -- integer; the desired bit length of the ouput

    OUTPUT:

    - the ``L``-bit vector representation of ``I``

    EXAMPLES::

        sage: from sage.crypto.block_cipher.present import convert_to_vector
        sage: convert_to_vector(0x1F, 8)
        (1, 1, 1, 1, 1, 0, 0, 0)
        sage: v = vector(GF(2), 4, [1,0,1,0])
        sage: convert_to_vector(v, 4)
        (1, 0, 1, 0)
    """
    try:
        state = vector(GF(2), L, ZZ(I).digits(2, padto=L))
        return state
    except TypeError:
        # ignore the error and try list-like types
        pass
    state = vector(GF(2), L, ZZ(list(I), 2).digits(2, padto=L))
    return state
