r"""
Enigma

This module implements the enigma and provides predefined rotors and
reflectors for the Enigma I. It also provides the ability to create custom
rotors and reflectors and does not limit the amount of rotors to use in the
enigma.


EXAMPLES:

The easiest way to create an Enigma is by using the predefined Enigma settings::

    sage: enigma = Enigma.create([5, 3, 2], [5, 24, 11], [2, 4, 24], 2)
    sage: one = enigma.en_de_crypt("GOODMORNING")
    sage: one
    'HKVOWVPLLBC'

You can also build a custom enigma. Here we build the enigma we just created by
hand::

    sage: a = Rotor([21, 25, 1, 17, 6, 8, 19, 24, 20, 15, 18, 3, 13, 7, 11, 23,
    ....: 0, 22, 12, 9, 16, 14, 5, 4, 2, 10], 26, 5, 2)
    sage: b = Rotor([1, 3, 5, 7, 9, 11, 2, 15, 17, 19, 23, 21, 25, 13, 24, 4,
    ....: 8, 22, 6, 0, 10, 12, 20, 18, 16, 14], 22, 24, 4)
    sage: c = Rotor([0, 9, 3, 10, 18, 8, 17, 20, 23, 1, 11, 7, 22, 19, 12, 2,
    ....: 16, 6, 25, 13, 15, 24, 5, 21, 14, 4], 5, 11, 24)
    sage: r = Reflector([24, 17, 20, 7, 16, 18, 11, 3, 15, 23, 13, 6, 14, 10,
    ....: 12, 8, 4, 1, 5, 25, 2, 22, 21, 9, 0, 19])
    sage: enigma_custom = Enigma([a, b, c], r)
    sage: two = enigma_custom.en_de_crypt("GOODMORNING")
    sage: two
    'HKVOWVPLLBC'
    sage: one == two
    True

To decrypt an encrypted message we need to know the initial settings of the
enigma used to encrypt the message. This can be done by just recreating that
enigma::

    sage: decrypt_enigma = Enigma.create([5, 3, 2], [5, 24, 11], [2, 4, 24], 2)
    sage: decrypted = decrypt_enigma.en_de_crypt("HKVOWVPLLBC")
    sage: decrypted
    'GOODMORNING'

If we encrypt multiple messages it may be stressful to recreate the exact state
of the enigma to decrypt a specific message. In that case we can duplicate the
enigma before we encrypt a message that we want to decrypt::

    sage: enigma.en_de_crypt("MORE")
    'QKYK'
    sage: enigma.en_de_crypt("WORDS")
    'VZGEV'
    sage: decrypt_enigma = enigma.duplicate()
    sage: encrypted = enigma.en_de_crypt("SECRET")
    sage: encrypted
    'NGZIFC'
    sage: decrypted = decrypt_enigma.en_de_crypt(encrypted)
    sage: decrypted
    'SECRET'

The enigma also supports a plugboard::

    sage: p = Reflector([1, 3, 13, 19, 24, 2, 25, 8, 0, 17, 18, 20, 12, 4, 9,
    ....: 11, 16, 21, 10, 6, 23, 15, 7, 22, 5, 14])
    sage: enigma = Enigma.create([3, 4, 1], [12, 23, 1], [1, 8, 7], 1, p)
    sage: decrypt_enigma = enigma.duplicate()
    sage: encrypted = enigma.en_de_crypt("PLUGBOARD")
    sage: encrypted
    'YHKVADMPH'
    sage: decrypted = decrypt_enigma.en_de_crypt("YHKVADMPH")
    sage: decrypted
    'PLUGBOARD'

It may also be helpful to take a closer look at the representation of the used
permutations. Each letter has a index beginning with 0 representing the letter
A. If we use a plugboard [0, 1, 2, ...] we map every letter to itself (A -> A,
B -> B, and so on)::

    sage: p = Reflector([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
    ....: 16, 17, 18, 19, 20, 21, 22, 23, 24, 25])
    sage: e_wo_p = Enigma.create([3, 4, 1], [12, 23, 1], [1, 8, 7], 1)
    sage: e_w_p = Enigma.create([3, 4, 1], [12, 23, 1], [1, 8, 7], 1, p)
    sage: one = e_wo_p.en_de_crypt("PERMUTATION")
    sage: one
    'KQFCZJUITNF'
    sage: two = e_w_p.en_de_crypt("PERMUTATION")
    sage: two
    'KQFCZJUITNF'
    sage: one == two
    True

AUTHORS:

- Nico Koeckerbauer and Benedek Kovacs (2024-06): initial version
- Tobias Schneider (2024-06): proofreading and testing it via implementation of Turing-Bombe

"""

# ****************************************************************************
#       Copyright (C) 2024 Nico Koeckerbauer <kolando@gmx.de>,
#                           Tobias Schneider <bartv4der@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


class Rotor:
    r"""
    An object representing a rotor of the enigma. The rotor of the enigma
    permutates each letter that passes through it. It rotates after the
    preceding wheel has rotated a fixed amount of times.

    INPUT:

    - ``permutation`` -- list[integer]; the permutation of the letters for this
      rotor, must have a length of 26 with values being a integer between 0 and
      25, also note that the permutation must be a bijection
    - ``notch`` -- integer; sets the position at which the rotor will cause the
      following rotor to rotate, must be a integer between 1 and 26
    - ``initial_position`` -- integer; sets the position at which the rotor
      will start, must be a integer between 1 and 26
    - ``offset`` -- integer; also known as ring position

    EXAMPLES::

        sage: r = Rotor([21, 25, 1, 17, 6, 8, 19, 24, 20, 15, 18, 3, 13, 7, 11,
        ....: 23, 0, 22, 12, 9, 16, 14, 5, 4, 2, 10], 26, 5, 2)
        sage: r
        Rotor with current position 5 (Initially 5), offset of 2, notch at
            26 and using permutation [21, 25, 1, 17, 6, 8, 19, 24, 20, 15, 18,
            3, 13, 7, 11, 23, 0, 22, 12, 9, 16, 14, 5, 4, 2, 10]

    The permutation used in the rotor needs to have a length of exactly 26::

        sage: r = Rotor([21, 25, 1, 17, 6, 8, 19, 24, 20, 15, 18, 3, 13, 7, 11,
        ....: 23, 0, 22, 12, 9, 16, 14], 26, 5, 2)
        Traceback (most recent call last):
        ...
        ValueError: Permutation must have a length of 26
        sage: r = Rotor([21, 25, 1, 17, 6, 8, 19, 24, 20, 15, 18, 3, 13, 7, 11,
        ....: 23, 0, 22, 12, 9, 16, 14, 12, 11, 10, 9, 8, 7, 5, 4], 26, 5, 2)
        Traceback (most recent call last):
        ...
        ValueError: Permutation must have a length of 26

    The notch must be a integer between 1 and 26::

        sage: r = Rotor([21, 25, 1, 17, 6, 8, 19, 24, 20, 15, 18, 3, 13, 7, 11,
        ....: 23, 0, 22, 12, 9, 16, 14, 5, 4, 2, 10], 50, 5, 2)
        Traceback (most recent call last):
        ...
        ValueError: Notch must be a integer between 1 and 26
        sage: r = Rotor([21, 25, 1, 17, 6, 8, 19, 24, 20, 15, 18, 3, 13, 7, 11,
        ....: 23, 0, 22, 12, 9, 16, 14, 5, 4, 2, 10], 0, 5, 2)
        Traceback (most recent call last):
        ...
        ValueError: Notch must be a integer between 1 and 26

    The initial position of the rotor must be a integer between 1 and 26::

        sage: r = Rotor([21, 25, 1, 17, 6, 8, 19, 24, 20, 15, 18, 3, 13, 7, 11,
        ....: 23, 0, 22, 12, 9, 16, 14, 5, 4, 2, 10], 26, 0, 2)
        Traceback (most recent call last):
        ...
        ValueError: Initial position must be a integer between 1 and 26
        sage: r = Rotor([21, 25, 1, 17, 6, 8, 19, 24, 20, 15, 18, 3, 13, 7, 11,
        ....: 23, 0, 22, 12, 9, 16, 14, 5, 4, 2, 10], 26, 30, 2)
        Traceback (most recent call last):
        ...
        ValueError: Initial position must be a integer between 1 and 26

    The values in the permutation must be a integer between 0 and 25::

        sage: r = Rotor([26, 25, 1, 17, 6, 8, 19, 24, 20, 15, 18, 3, 13, 7, 11,
        ....: 23, 0, 22, 12, 9, 16, 14, 5, 4, 2, 10], 26, 1, 2)
        Traceback (most recent call last):
        ...
        ValueError: Values of the permutation must be a integer between 0 and 25
        sage: r = Rotor([-1, 25, 1, 17, 6, 8, 19, 24, 20, 15, 18, 3, 13, 7, 11,
        ....: 23, 0, 22, 12, 9, 16, 14, 5, 4, 2, 10], 26, 1, 2)
        Traceback (most recent call last):
        ...
        ValueError: Values of the permutation must be a integer between 0 and 25

    The permutation also must be a bijection, meaning that no value in the
    permutation can occur more than one time::

        sage: r = Rotor([21, 25, 1, 17, 6, 8, 19, 24, 21, 15, 18, 3, 13, 7, 11,
        ....: 23, 0, 22, 12, 9, 16, 14, 5, 4, 2, 10], 26, 1, 2)
        Traceback (most recent call last):
        ...
        ValueError: Permutation must be a bijection
    """

    def __init__(
        self, permutation: list[int], notch: int, initial_position: int, offset: int
    ):
        r"""
        Initializes the ``Rotor`` class. See the class :class:`Rotor` for full
        documentation on the input of this initialization method.

        EXAMPLES::

        sage: r = Rotor([21, 25, 1, 17, 6, 8, 19, 24, 20, 15, 18, 3, 13, 7, 11,
        ....: 23, 0, 22, 12, 9, 16, 14, 5, 4, 2, 10], 26, 5, 2)
        sage: r
        Rotor with current position 5 (Initially 5), offset of 2, notch at
            26 and using permutation [21, 25, 1, 17, 6, 8, 19, 24, 20, 15, 18,
            3, 13, 7, 11, 23, 0, 22, 12, 9, 16, 14, 5, 4, 2, 10]
        """

        # list used to check if permutation is a bijection
        check_list: list[int] = [int] * 26
        if len(permutation) != 26:
            raise ValueError("Permutation must have a length of 26")
        if notch > 26 or notch < 1:
            raise ValueError("Notch must be a integer between 1 and 26")
        if initial_position > 26 or initial_position < 1:
            raise ValueError("Initial position must be a integer between 1 and " "26")
        for v in permutation:
            if v < 0 or v > 25:
                raise ValueError(
                    "Values of the permutation must be a integer " "between 0 and 25"
                )
            if check_list[v] == 1:
                raise ValueError("Permutation must be a bijection")
            else:
                check_list[v] = 1

        self._position = initial_position - 1
        self._initial_position = initial_position - 1
        self._notch = notch - 1
        self._permutation = permutation
        self._offset = offset - 1

    def _process_letter_forward(self, letter: str, rotate: bool) -> tuple[str, bool]:
        r"""
        Private function which allows a rotor to process a single letter when
        moving towards the reflector of the enigma

        INPUT:

        - ``letter`` -- string; the letter to be processed
        - ``rotate`` -- boolean; indicates whether the rotor should rotate

        OUTPUT: a tuple of

        - the processed ``letter``
        - a boolean indicating whether the next rotor should rotate, this
          happens when the current `position` of the rotor equals the ``notch``

        EXAMPLES::

            sage: r = Rotor([21, 25, 1, 17, 6, 8, 19, 24, 20, 15, 18, 3, 13, 7,
            ....: 11, 23, 0, 22, 12, 9, 16, 14, 5, 4, 2, 10], 26, 5, 2)
            sage: r._process_letter_forward("B", False)
            ('D', False)

        Can only process letters A - Z::

            sage: r._process_letter_forward("1", False)
            Traceback (most recent call last):
            ...
            ValueError: Invalid letter: Input can only contain letters A - Z

        Letters need to be uppercase::

            sage: r._process_letter_forward("a", False)
            Traceback (most recent call last):
            ...
            ValueError: Invalid letter: Input can only contain letters A - Z
        """

        letter_value = ord(letter) - 65
        next_rotate: bool = False

        if letter_value < 0 or letter_value > 25:
            raise ValueError("Invalid letter: Input can only contain letters A - Z")

        if rotate:
            self._position = (self._position + 1) % 26

        if self._position == self._notch:
            next_rotate: bool = True

        current_offset = (self._position - self._offset) % 26
        letter_value = (letter_value + current_offset) % 26
        letter_value = (self._permutation[letter_value] - current_offset) % 26

        return chr(letter_value + 65), next_rotate

    def _process_letter_backward(self, letter: str) -> str:
        r"""
        Private function which allows a rotor to process a single letter when
        moving through the enigma,
        after passing the reflector.

        INPUT:

        - ``letter`` -- string; letter to be processed

        OUTPUT: processed letter

        EXAMPLES::

            sage: r = Rotor([21, 25, 1, 17, 6, 8, 19, 24, 20, 15, 18, 3, 13, 7,
            ....: 11, 23, 0, 22, 12, 9, 16, 14, 5, 4, 2, 10], 26, 5, 2)
            sage: r._process_letter_backward("B")
            'U'

        Can only process letters A - Z::

            sage: r._process_letter_backward("1")
            Traceback (most recent call last):
            ...
            ValueError: Invalid letter: Input can only contain letters A - Z

        Letters need to be uppercase::

            sage: r._process_letter_backward("a")
            Traceback (most recent call last):
            ...
            ValueError: Invalid letter: Input can only contain letters A - Z
        """

        letter_value = ord(letter) - 65

        if letter_value < 0 or letter_value > 25:
            raise ValueError("Invalid letter: Input can only contain letters A - Z")

        current_offset = (self._position - self._offset) % 26
        letter_value = (letter_value + current_offset) % 26
        letter_value = (self._permutation.index(letter_value) - current_offset) % 26

        return chr(letter_value + 65)

    def duplicate(self):
        r"""
        Creates a copy of the current rotor. Note that it will not be a 1:1
        copy of the rotor, since the initial position of the copy will be set
        to the current position of this rotor

        OUTPUT: Copy of this rotor

        EXAMPLES::

            sage: r = Rotor([21, 25, 1, 17, 6, 8, 19, 24, 20, 15, 18, 3, 13, 7,
            ....: 11, 23, 0, 22, 12, 9, 16, 14, 5, 4, 2, 10], 26, 5, 2)
            sage: r._process_letter_forward('A', True)
            ('C', False)
            sage: c = r.duplicate()
            sage: r
            Rotor with current position 6 (Initially 5), offset of 2, notch at
            26 and using permutation [21, 25, 1, 17, 6, 8, 19, 24, 20, 15, 18,
            3, 13, 7, 11, 23, 0, 22, 12, 9, 16, 14, 5, 4, 2, 10]
        """

        return Rotor(
            self._permutation, self._notch + 1, self._position + 1, self._offset + 1
        )

    def __repr__(self) -> str:
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: r = Rotor([21, 25, 1, 17, 6, 8, 19, 24, 20, 15, 18, 3, 13, 7,
            ....: 11, 23, 0, 22, 12, 9, 16, 14, 5, 4, 2, 10], 26, 5, 2)
            sage: r
            Rotor with current position 5 (Initially 5), offset of 2, notch at
            26 and using permutation [21, 25, 1, 17, 6, 8, 19, 24, 20, 15, 18,
            3, 13, 7, 11, 23, 0, 22, 12, 9, 16, 14, 5, 4, 2, 10]
        """

        msg = (
            "Rotor with current position {0} (Initially {1}), "
            "offset of {2}, notch at {3} and using permutation {4}"
        )
        return msg.format(
            self._position + 1,
            self._initial_position + 1,
            self._offset + 1,
            self._notch + 1,
            self._permutation,
        )


class Reflector:
    r"""
    A special rotor that misses the key thing of one: the rotating, which makes
    it a simple permutation.

    This class fulfills the purpose of the needed reflector inside the enigma
    but can also act as the plugboard since the technical requirements are
    almost the same.
    In the case of being the actual reflector it gets used once a letter has
    passed every rotor for the first time, reflecting it back to the rotors
    after which the letter passes every rotor again.
    In the case of being the plugboard it gets used twice. Once before entering
    the first rotor and once after being reflected and passing the first rotor
    again.
    In both cases it permutates the passing letter.

    INPUT:

    - ``permutation`` -- list[integer]; permutation to be used in the reflector
      which must have a length of 26, note that every value in the permutation
      has to be a integer between 0 and 25 and the permutation itself must be
      a bijection

    EXAMPLES::

        sage: r = Reflector([24, 17, 20, 7, 16, 18, 11, 3, 15, 23, 13, 6, 14,
        ....: 10, 12, 8, 4, 1, 5, 25, 2, 22, 21, 9, 0, 19])
        sage: r
        Reflector or Plugboard using permutation [24, 17, 20, 7, 16, 18, 11, 3,
        15, 23, 13, 6, 14, 10, 12, 8, 4, 1, 5, 25, 2, 22, 21, 9, 0, 19]

    The reflector/plugboard needs to have a length of exactly 26::

        sage: r = Reflector([24, 17, 20, 7, 16, 18, 11, 3, 15, 23, 13, 6, 14,
        ....: 10, 12, 8, 4, 1, 5, 25, 2, 22, 21])
        Traceback (most recent call last):
        ...
        ValueError: Permutation must have a length of 26
        sage: p = Reflector([24, 17, 20, 7, 16, 18, 11, 3, 15, 23, 13, 6, 14,
        ....: 10, 12, 8, 4, 1, 5, 25, 2, 22, 21, 9, 0, 19, 2, 1 , 1])
        Traceback (most recent call last):
        ...
        ValueError: Permutation must have a length of 26

    Note that values used in the reflector/plugboard must be a integer
    between 0 and 25::

        sage: r = Reflector([50, 17, 20, 7, 16, 18, 11, 3, 15, 23, 13, 6, 14,
        ....: 10, 12, 8, 4, 1, 5, 25, 2, 22, 21, 12, 2, 19])
        Traceback (most recent call last):
        ...
        ValueError: Values of the permutation must be a integer between 0 and 25
        sage: r = Reflector([-12, 17, 20, 7, 16, 18, 11, 3, 15, 23, 13, 6, 14,
        ....: 10, 12, 8, 4, 1, 5, 25, 2, 22, 21, 12, 2, 19])
        Traceback (most recent call last):
        ...
        ValueError: Values of the permutation must be a integer between 0 and 25

    The permutation for a reflector/plugboard must be a bijection, meaning that
    no value in the permutation can occur more than one time::

        sage: r = Reflector([24, 17, 20, 7, 16, 17, 11, 3, 15, 23, 13, 6, 14,
        ....: 10, 12, 8, 4, 1, 5, 25, 2, 22, 21, 9, 0, 19])
        Traceback (most recent call last):
        ...
        ValueError: Permutation must be a bijection
    """

    def __init__(self, permutation: list[int]):
        r"""
        Initializes the ``Reflector`` class. See the class :class:`Reflector`
        for full documentation on the input of this initialization method.

        EXAMPLES::

        sage: r = Reflector([24, 17, 20, 7, 16, 18, 11, 3, 15, 23, 13, 6, 14,
        ....: 10, 12, 8, 4, 1, 5, 25, 2, 22, 21, 9, 0, 19])
        """

        # list used to check if permutation is a bijection
        check_list: list[int] = [int] * 26
        if len(permutation) != 26:
            raise ValueError("Permutation must have a length of 26")
        for v in permutation:
            if v < 0 or v > 25:
                raise ValueError(
                    "Values of the permutation must be a integer between 0 and 25"
                )
            if check_list[v] == 1:
                raise ValueError("Permutation must be a bijection")
            else:
                check_list[v] = 1

        self._permutation = permutation

    def _process_letter(self, letter: str, forward: bool = True) -> str:
        r"""
        Private function to permutate a given letter

        INPUT:

        - ``letter`` -- string; letter to be processed
        - ``forward`` -- boolean (default: ``True``); indicating the direction of travel, this is
          only used when the reflector is used as a plugboard

        OUTPUT: processed letter

        EXAMPLES::

            sage: r = Reflector([24, 17, 20, 7, 16, 18, 11, 3, 15, 23, 13, 6,
            ....: 14, 10, 12, 8, 4, 1, 5, 25, 2, 22, 21, 9, 0, 19])
            sage: r._process_letter("A")
            'Y'

        When being used as a plugboard we need to change the direction of
        travel when the letter passes the second time::

            sage: r = Reflector([24, 17, 20, 7, 16, 18, 11, 3, 15, 23, 13, 6,
            ....: 14, 10, 12, 8, 4, 1, 5, 25, 2, 22, 21, 9, 0, 19])
            sage: r._process_letter("Y", False)
            'A'

        Can only process letters A - Z::

            sage: r._process_letter("1")
            Traceback (most recent call last):
            ...
            ValueError: Invalid letter: Input can only contain letters A - Z

        Letters need to be uppercase::

            sage: r._process_letter("a")
            Traceback (most recent call last):
            ...
            ValueError: Invalid letter: Input can only contain letters A - Z
        """

        letter_value = ord(letter) - 65

        if letter_value < 0 or letter_value > 25:
            raise ValueError("Invalid letter: Input can only contain letters A - Z")

        if forward:
            letter_value = self._permutation[letter_value]
        else:
            letter_value = self._permutation.index(letter_value)

        return chr(letter_value + 65)

    def __repr__(self) -> str:
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: r = Reflector([24, 17, 20, 7, 16, 18, 11, 3, 15, 23, 13, 6,
            ....: 14, 10, 12, 8, 4, 1, 5, 25, 2, 22, 21, 9, 0, 19])
            sage: r
            Reflector or Plugboard using permutation [24, 17, 20, 7, 16, 18,
            11, 3, 15, 23, 13, 6, 14, 10, 12, 8, 4, 1, 5, 25, 2, 22, 21, 9, 0,
            19]
        """

        msg = "Reflector or Plugboard using permutation {0}"
        return msg.format(self._permutation)


class Enigma:
    r"""
    An object representing the Enigma. The enigma is used to encrypt and
    decrypt messages.
    To decrypt a message a enigma with the same settings as the one used to
    encrypt the message is needed.

    INPUT:

    - ``rotors`` -- list[Rotor]; a list of rotors to be used in this enigma
    - ``reflector`` -- Reflector; the reflector to be used in this enigma
    - ``plugboard`` -- Reflector (default: ``None``); optional plugboard to use
      in this enigma

    EXAMPLES::

        sage: a = Rotor([21, 25, 1, 17, 6, 8, 19, 24, 20, 15, 18, 3, 13,
        ....: 7, 11, 23, 0, 22, 12, 9, 16, 14, 5, 4, 2, 10], 26, 5, 2)
        sage: b = Rotor([1, 3, 5, 7, 9, 11, 2, 15, 17, 19, 23, 21, 25, 13,
        ....: 24, 4, 8, 22, 6, 0, 10, 12, 20, 18, 16, 14], 22, 24, 4)
        sage: c = Rotor([0, 9, 3, 10, 18, 8, 17, 20, 23, 1, 11, 7, 22, 19,
        ....: 12, 2, 16, 6, 25, 13, 15, 24, 5, 21, 14, 4], 5, 11, 24)
        sage: d = Reflector([24, 17, 20, 7, 16, 18, 11, 3, 15, 23, 13, 6,
        ....: 14, 10, 12, 8, 4, 1, 5, 25, 2, 22, 21, 9, 0, 19])
        sage: e = Enigma([a, b, c], d, None)
        sage: e
        Enigma using [Rotor with current position 5 (Initially 5), offset of
        2, notch at 26 and using permutation [21, 25, 1, 17, 6, 8, 19, 24, 20,
        15, 18, 3, 13, 7, 11, 23, 0, 22, 12, 9, 16, 14, 5, 4, 2, 10], Rotor
        with current position 24 (Initially 24), offset of 4, notch at 22 and
        using permutation [1, 3, 5, 7, 9, 11, 2, 15, 17, 19, 23, 21, 25, 13,
        24, 4, 8, 22, 6, 0, 10, 12, 20, 18, 16, 14], Rotor with current
        position 11 (Initially 11), offset of 24, notch at 5 and using
        permutation [0, 9, 3, 10, 18, 8, 17, 20, 23, 1, 11, 7, 22, 19, 12, 2,
        16, 6, 25, 13, 15, 24, 5, 21, 14, 4]], Reflector or Plugboard using
        permutation [24, 17, 20, 7, 16, 18, 11, 3, 15, 23, 13, 6, 14, 10, 12,
        8, 4, 1, 5, 25, 2, 22, 21, 9, 0, 19] as a reflector and None as a
        plugboard

    .. NOTE::

        To use predefined rotors and reflectors see :meth:`create`
    """

    def __init__(
        self,
        rotors: list[Rotor],
        reflector: Reflector,
        plugboard: Reflector = None,
    ):
        r"""
        Initializes the ``Enigma`` class. See the class :class:`Enigma` for
        full documentation on the input of this initialization method.

        EXAMPLES::

        sage: a = Rotor([21, 25, 1, 17, 6, 8, 19, 24, 20, 15, 18, 3, 13,
        ....: 7, 11, 23, 0, 22, 12, 9, 16, 14, 5, 4, 2, 10], 26, 5, 2)
        sage: b = Rotor([1, 3, 5, 7, 9, 11, 2, 15, 17, 19, 23, 21, 25, 13,
        ....: 24, 4, 8, 22, 6, 0, 10, 12, 20, 18, 16, 14], 22, 24, 4)
        sage: c = Rotor([0, 9, 3, 10, 18, 8, 17, 20, 23, 1, 11, 7, 22, 19,
        ....: 12, 2, 16, 6, 25, 13, 15, 24, 5, 21, 14, 4], 5, 11, 24)
        sage: d = Reflector([24, 17, 20, 7, 16, 18, 11, 3, 15, 23, 13, 6,
        ....: 14, 10, 12, 8, 4, 1, 5, 25, 2, 22, 21, 9, 0, 19])
        sage: e = Enigma([a, b, c], d, None)
        sage: e
        Enigma using [Rotor with current position 5 (Initially 5), offset of 2,
        notch at 26 and using permutation [21, 25, 1, 17, 6, 8, 19, 24, 20, 15, 18,
        3, 13, 7, 11, 23, 0, 22, 12, 9, 16, 14, 5, 4, 2, 10], Rotor with current
        position 24 (Initially 24), offset of 4, notch at 22 and using permutation
        [1, 3, 5, 7, 9, 11, 2, 15, 17, 19, 23, 21, 25, 13, 24, 4, 8, 22, 6, 0, 10, 12,
        20, 18, 16, 14], Rotor with current position 11 (Initially 11), offset of 24,
        notch at 5 and using permutation [0, 9, 3, 10, 18, 8, 17, 20, 23, 1, 11, 7, 22,
        19, 12, 2, 16, 6, 25, 13, 15, 24, 5, 21, 14, 4]], Reflector or Plugboard using
        permutation [24, 17, 20, 7, 16, 18, 11, 3, 15, 23, 13, 6, 14, 10, 12, 8, 4,
        1, 5, 25, 2, 22, 21, 9, 0, 19] as a reflector and None as a plugboard
        """

        self._rotors = rotors
        self._reflector = reflector
        self._plugboard = plugboard

    @staticmethod
    def create(
        rotors: list[int],
        initial_positions: list[int],
        offsets: list[int],
        reflector: int,
        plugboard: Reflector = None,
    ):
        r"""
        Function to create an Enigma with predefined rotors and reflectors.
        All the predefined rotors and reflectors are replicated from the
        original Enigma I.

        INPUT:

        - ``rotors`` -- list[integer]; List of integers with each integer
          indicating which rotor to use, must be a value between 1 and 5

        - ``initial_positions`` -- list[integer]; List of integers setting the
          initial position of the ``rotors``. This list must be the same length
          as ``rotors`` since the initial position at index ``i`` relates to
          the rotor at index ``i``. Note that every value in this list has to
          be a integer between 1 and 26

        - ``offsets`` -- list[integer]; List of integers setting the offset of
          ``rotors``. This list must be the same length as ``rotors`` since the
          offset at index ``i`` relates to the rotor at index ``i``

        - ``reflector`` -- integer; Select one of three predefined reflectors
          used in the Enigma I

        - ``plugboard`` -- Reflector (default: None); Optional plugboard to be
          used in this enigma

        OUTPUT: Enigma with defined settings

        EXAMPLES::

            sage: e = Enigma.create([5, 3, 2], [5, 24, 11], [2, 4, 24], 2, None)
            sage: e
            Enigma using [Rotor with current position 5 (Initially 5), offset
            of 2, notch at 26 and using permutation [21, 25, 1, 17, 6, 8, 19,
            24, 20, 15, 18, 3, 13, 7, 11, 23, 0, 22, 12, 9, 16, 14, 5, 4, 2,
            10], Rotor with current position 24 (Initially 24), offset of 4,
            notch at 22 and using permutation [1, 3, 5, 7, 9, 11, 2, 15, 17,
            19, 23, 21, 25, 13, 24, 4, 8, 22, 6, 0, 10, 12, 20, 18, 16, 14],
            Rotor with current position 11 (Initially 11), offset of 24, notch
            at 5 and using permutation [0, 9, 3, 10, 18, 8, 17, 20, 23, 1, 11,
            7, 22, 19, 12, 2, 16, 6, 25, 13, 15, 24, 5, 21, 14, 4]], Reflector
            or Plugboard using permutation [24, 17, 20, 7, 16, 18, 11, 3, 15,
            23, 13, 6, 14, 10, 12, 8, 4, 1, 5, 25, 2, 22, 21, 9, 0, 19] as a
            reflector and None as a plugboard

        Although the original Enigma I did only have 3 rotors, nothing stops
        you from using more rotors than that::

            sage: e = Enigma.create([1, 2, 3, 4, 5, 4, 3, 2, 1],
            ....: [3, 6, 9, 12, 15, 18, 21, 24, 26],
            ....: [2, 4, 6, 8, 10, 12, 14, 16, 18], 3)
            sage: e
            Enigma using [Rotor with current position 3 (Initially 3), offset
            of 2, notch at 17 and using permutation [4, 10, 12, 5, 11, 6, 3,
            16, 21, 25, 13, 19, 14, 22, 24, 7, 23, 20, 18, 15, 0, 8, 1, 17, 2,
            9], Rotor with current position 6 (Initially 6), offset of 4, notch
            at 5 and using permutation [0, 9, 3, 10, 18, 8, 17, 20, 23, 1, 11,
            7, 22, 19, 12, 2, 16, 6, 25, 13, 15, 24, 5, 21, 14, 4], Rotor with
            current position 9 (Initially 9), offset of 6, notch at 22 and
            using permutation [1, 3, 5, 7, 9, 11, 2, 15, 17, 19, 23, 21, 25,
            13, 24, 4, 8, 22, 6, 0, 10, 12, 20, 18, 16, 14], Rotor with current
            position 12 (Initially 12), offset of 8, notch at 10 and using
            permutation [4, 18, 14, 21, 15, 25, 9, 0, 24, 16, 20, 8, 17, 7, 23,
            11, 13, 5, 19, 6, 10, 3, 2, 12, 22, 1], Rotor with current position
            15 (Initially 15), offset of 10, notch at 26 and using permutation
            [21, 25, 1, 17, 6, 8, 19, 24, 20, 15, 18, 3, 13, 7, 11, 23, 0, 22,
            12, 9, 16, 14, 5, 4, 2, 10], Rotor with current position 18
            (Initially 18), offset of 12, notch at 10 and using permutation
            [4, 18, 14, 21, 15, 25, 9, 0, 24, 16, 20, 8, 17, 7, 23, 11, 13, 5,
            19, 6, 10, 3, 2, 12, 22, 1], Rotor with current position 21
            (Initially 21), offset of 14, notch at 22 and using permutation
            [1, 3, 5, 7, 9, 11, 2, 15, 17, 19, 23, 21, 25, 13, 24, 4, 8, 22,
            6, 0, 10, 12, 20, 18, 16, 14], Rotor with current position 24
            (Initially 24), offset of 16, notch at 5 and using permutation
            [0, 9, 3, 10, 18, 8, 17, 20, 23, 1, 11, 7, 22, 19, 12, 2, 16, 6,
            25, 13, 15, 24, 5, 21, 14, 4], Rotor with current position 26
            (Initially 26), offset of 18, notch at 17 and using permutation
            [4, 10, 12, 5, 11, 6, 3, 16, 21, 25, 13, 19, 14, 22, 24, 7, 23,
            20, 18, 15, 0, 8, 1, 17, 2, 9]], Reflector or Plugboard using
            permutation [5, 21, 15, 9, 8, 0, 14, 24, 4, 3, 17, 25, 23, 22, 6,
            2, 19, 10, 20, 16, 18, 1, 13, 12, 7, 11] as a reflector and None
            as a plugboard

        The length of the list for the initials position and the list for the
        offsets must be the same length as the list for the rotors to be used
        since every rotor need a initial position and offset::

            sage: e = Enigma.create([5, 3, 2, 4], [5, 24], [2, 4, 24], 2, None)
            Traceback (most recent call last):
            ...
            ValueError: initial_positions and rotors must be the same length
            sage: e = Enigma.create([5, 3, 2], [5, 24, 1], [2, 4], 2, None)
            Traceback (most recent call last):
            ...
            ValueError: offsets and rotors must be the same length

        There are not unlimited predefined rotors and reflectors::

            sage: e = Enigma.create([2, 6], [2, 5], [12, 13], 2, None)
            Traceback (most recent call last):
            ...
            ValueError: Values in rotors must be a integer between 1 and 5
            sage: e = Enigma.create([2, 5], [2, 5], [12, 13], 4, None)
            Traceback (most recent call last):
            ...
            ValueError: Reflector must be a value between 1 and 3
        """

        wheels = [
            [
                4,
                10,
                12,
                5,
                11,
                6,
                3,
                16,
                21,
                25,
                13,
                19,
                14,
                22,
                24,
                7,
                23,
                20,
                18,
                15,
                0,
                8,
                1,
                17,
                2,
                9,
            ],
            [
                0,
                9,
                3,
                10,
                18,
                8,
                17,
                20,
                23,
                1,
                11,
                7,
                22,
                19,
                12,
                2,
                16,
                6,
                25,
                13,
                15,
                24,
                5,
                21,
                14,
                4,
            ],
            [
                1,
                3,
                5,
                7,
                9,
                11,
                2,
                15,
                17,
                19,
                23,
                21,
                25,
                13,
                24,
                4,
                8,
                22,
                6,
                0,
                10,
                12,
                20,
                18,
                16,
                14,
            ],
            [
                4,
                18,
                14,
                21,
                15,
                25,
                9,
                0,
                24,
                16,
                20,
                8,
                17,
                7,
                23,
                11,
                13,
                5,
                19,
                6,
                10,
                3,
                2,
                12,
                22,
                1,
            ],
            [
                21,
                25,
                1,
                17,
                6,
                8,
                19,
                24,
                20,
                15,
                18,
                3,
                13,
                7,
                11,
                23,
                0,
                22,
                12,
                9,
                16,
                14,
                5,
                4,
                2,
                10,
            ],
        ]
        notches = [17, 5, 22, 10, 26]
        reflectors = [
            [
                4,
                9,
                12,
                25,
                0,
                11,
                24,
                23,
                21,
                1,
                22,
                5,
                2,
                17,
                16,
                20,
                14,
                13,
                19,
                18,
                15,
                8,
                10,
                7,
                6,
                3,
            ],
            [
                24,
                17,
                20,
                7,
                16,
                18,
                11,
                3,
                15,
                23,
                13,
                6,
                14,
                10,
                12,
                8,
                4,
                1,
                5,
                25,
                2,
                22,
                21,
                9,
                0,
                19,
            ],
            [
                5,
                21,
                15,
                9,
                8,
                0,
                14,
                24,
                4,
                3,
                17,
                25,
                23,
                22,
                6,
                2,
                19,
                10,
                20,
                16,
                18,
                1,
                13,
                12,
                7,
                11,
            ],
        ]
        rotor_list: list[Rotor] = []

        if len(rotors) < 0:
            raise ValueError("Enigma must have a rotor")
        if len(initial_positions) != len(rotors):
            raise ValueError("initial_positions and rotors must be the same length")
        if len(offsets) != len(rotors):
            raise ValueError("offsets and rotors must be the same length")
        if reflector < 1 or reflector > 3:
            raise ValueError("Reflector must be a value between 1 and 3")

        for idx, r in enumerate(rotors):
            if r < 0 or r > 5:
                raise ValueError("Values in rotors must be a integer between 1 and 5")
            rotor_list.append(
                Rotor(
                    wheels[r - 1], notches[r - 1], initial_positions[idx], offsets[idx]
                )
            )

        return Enigma(rotor_list, Reflector(reflectors[reflector - 1]), plugboard)

    def en_de_crypt(self, text):
        r"""
        Processes the `text` which results in receiving a encrypted or
        decrypted text.
        To decrypt a text, a enigma is needed with the same starting settings
        as the enigma used to encrypt the text.

        INPUT:

        - ``text`` -- string; The text to be processed by this enigma, can only
          contain letters A - Z

        OUTPUT: the encrypted or decrypted text

        EXAMPLES::

            sage: a = Rotor([21, 25, 1, 17, 6, 8, 19, 24, 20, 15, 18, 3, 13,
            ....: 7, 11, 23, 0, 22, 12, 9, 16, 14, 5, 4, 2, 10], 26, 5, 2)
            sage: b = Rotor([1, 3, 5, 7, 9, 11, 2, 15, 17, 19, 23, 21, 25, 13,
            ....: 24, 4, 8, 22, 6, 0, 10, 12, 20, 18, 16, 14], 22, 24, 4)
            sage: c = Rotor([0, 9, 3, 10, 18, 8, 17, 20, 23, 1, 11, 7, 22, 19,
            ....: 12, 2, 16, 6, 25, 13, 15, 24, 5, 21, 14, 4], 5, 11, 24)
            sage: d = Reflector([24, 17, 20, 7, 16, 18, 11, 3, 15, 23, 13, 6,
            ....: 14, 10, 12, 8, 4, 1, 5, 25, 2, 22, 21, 9, 0, 19])
            sage: e = Enigma([a, b, c], d, None)
            sage: e.en_de_crypt("HELLO")
            'GCIID'

        The text can not contain any other character than the letters A - Z::

            sage: e.en_de_crypt("HELLO HOW ARE YOU")
            Traceback (most recent call last):
            ...
            ValueError: Invalid letter: Input can only contain letters A - Z

        However, lowercase letters are automatically converted to uppercase
        letters::

            sage: e.en_de_crypt("hellohowareyou")
            'XPVIDDZIXUIDVA'
        """

        result = ""
        message = str(text).upper()

        for letter in message:
            predecessor_tuple: tuple[str, bool] = letter, True

            if self._plugboard is not None:
                predecessor_tuple = (
                    self._plugboard._process_letter(predecessor_tuple[0]),
                    True,
                )

            for rotor in self._rotors:
                predecessor_tuple = rotor._process_letter_forward(
                    predecessor_tuple[0], predecessor_tuple[1]
                )

            predecessor_tuple = (
                self._reflector._process_letter(predecessor_tuple[0]),
                False,
            )

            for rotor in reversed(self._rotors):
                predecessor_tuple = (
                    rotor._process_letter_backward(predecessor_tuple[0]),
                    False,
                )

            if self._plugboard is not None:
                predecessor_tuple = (
                    self._plugboard._process_letter(predecessor_tuple[0], False),
                    False,
                )

            result += predecessor_tuple[0]

        return result

    def duplicate(self):
        r"""
        Creates a copy of the current enigma. Note that it will not be a 1:1
        copy of the enigma, since the initial positions of the new rotors will
        be the current positions of the rotors in this enigma.

        OUTPUT: Copy of this enigma

        EXAMPLES::

            sage: e = Enigma.create([5, 3, 2], [5, 24, 11], [2, 4, 24], 2)
            sage: c = e.duplicate()
            sage: text = "Hello".upper()
            sage: encrypted_text = e.en_de_crypt(text)
            sage: decrypted_text = c.en_de_crypt(encrypted_text)
            sage: decrypted_text == text
            True
        """

        # since rotors have a own state (position) we need a deep copy
        new_rotors: list[Rotor] = []
        for rotor in self._rotors:
            new_rotors.append(rotor.duplicate())

        # plugboard and reflector are static permutations anyway
        return Enigma(new_rotors, self._reflector, self._plugboard)

    def __repr__(self) -> str:
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: a = Rotor([21, 25, 1, 17, 6, 8, 19, 24, 20, 15, 18, 3, 13,
            ....: 7, 11, 23, 0, 22, 12, 9, 16, 14, 5, 4, 2, 10], 26, 5, 2)
            sage: b = Rotor([1, 3, 5, 7, 9, 11, 2, 15, 17, 19, 23, 21, 25, 13,
            ....: 24, 4, 8, 22, 6, 0, 10, 12, 20, 18, 16, 14], 22, 24, 4)
            sage: c = Rotor([0, 9, 3, 10, 18, 8, 17, 20, 23, 1, 11, 7, 22, 19,
            ....: 12, 2, 16, 6, 25, 13, 15, 24, 5, 21, 14, 4], 5, 11, 24)
            sage: d = Reflector([24, 17, 20, 7, 16, 18, 11, 3, 15, 23, 13, 6,
            ....: 14, 10, 12, 8, 4, 1, 5, 25, 2, 22, 21, 9, 0, 19])
            sage: e = Enigma([a, b, c], d, None)
            sage: e
            Enigma using [Rotor with current position 5 (Initially 5), offset
            of 2, notch at 26 and using permutation [21, 25, 1, 17, 6, 8, 19,
            24, 20, 15, 18, 3, 13, 7, 11, 23, 0, 22, 12, 9, 16, 14, 5, 4, 2,
            10], Rotor with current position 24 (Initially 24), offset of 4,
            notch at 22 and using permutation [1, 3, 5, 7, 9, 11, 2, 15, 17,
            19, 23, 21, 25, 13, 24, 4, 8, 22, 6, 0, 10, 12, 20, 18, 16, 14],
            Rotor with current position 11 (Initially 11), offset of 24, notch
            at 5 and using permutation [0, 9, 3, 10, 18, 8, 17, 20, 23, 1, 11,
            7, 22, 19, 12, 2, 16, 6, 25, 13, 15, 24, 5, 21, 14, 4]], Reflector
            or Plugboard using permutation [24, 17, 20, 7, 16, 18, 11, 3, 15,
            23, 13, 6, 14, 10, 12, 8, 4, 1, 5, 25, 2, 22, 21, 9, 0, 19] as a
            reflector and None as a plugboard
        """

        msg = "Enigma using {0}, {1} as a reflector and {2} as a plugboard"
        return msg.format(self._rotors, self._reflector, self._plugboard)
