r"""
This file contains a classical implementation of the Playfair cipher

NOTE::
    - I did not inherit the symmetric key cipher class because it required
      us to define a keyspace, and the keyspace contains 25! 5x5 matrices.
    - Helper functions that typically would be written with leading
      underscores were not, because we wanted to allow use of this class
      for pedagogical reasons so one can investigate the inner workings of
      the cipher.

AUTHORS:

- Amy Feaver (2026-06-12): initial version of Playfair Cipher, Sage Days 131
- Brian Heckel (2026-06-12): project group lead, Sage Days 131
- Laura Maddison (2026-06-12): helped in discussing code design
- James Bui (2026-06-12): helped in discussing code design
- Alasdair McAndrew (2010-03-19): we looked at McAndrew's code as published
                                  here https://github.com/sagemath/sage/issues/8559
"""

#*****************************************************************************
#       Copyright (C) 2026 Amy Feaver <mathfeaver@gmail.com>
#                Brian Heckel <heckelbri@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


class PlayfairCryptosystem():
    r"""
    Let `A = \{ a_0, a_1, a_2, \dots, a_{n-1} \}` be a non-empty alphabet
    consisting of `n` unique elements.

    In the classical Playfair cipher, the alphabet is the capital letters
    of the English alphabet with the exception that J and I are identified
    (J is replaced with I before encryption and in the keyword). This
    results in an alphabet of 25 symbols which can be then arranged into a
    `5 \times 5` key square according to a keyword. The keyword is written
    first, omitting repeated letters, and the remaining symbols are then
    inserted in alphabetical order to fill the square.

    A key of the Playfair cipher is therefore a `5 \times 5` arrangement
    of the symbols of the alphabet (without J), with each symbol appearing
    exactly once. Each symbol can be identified with an ordered pair of
    coordinates `(r,c)`, `0 \leq r,c < 5`. Thus the key square defines a
    bijection between the symbols of the alphabet and the set of coordinate
    pairs in the square.

    Unlike ciphers that operate on individual characters, the Playfair
    cipher encrypts pairs of letters, called digraphs. Before encryption,
    the plaintext is first stripped of all non-alphabetic characters, all
    instances of J are replaced with I, and the text is divided into digraphs.
    If a digraph contains two identical letters, a filler symbol is inserted
    between them. If the plaintext contains an odd number of symbols, a filler
    symbol is appended to the final digraph. The filler symbol is usually X,
    unless that would result in a digraph of XX, in which case the filler
    symbol is Z instead.

    Let `(p_1,p_2)` be a plaintext digraph and let `(r_1,c_1)` and `(r_2,c_2)`
    denote the coordinates of `p_1` and `p_2` in the key square. The
    corresponding ciphertext digraph `(c_1,c_2)` is obtained according to the
    following rules:

    - If `r_1 = r_2`, that is, the two symbols lie in the same row, each
        symbol is replaced by the symbol immediately to its right, wrapping
        around to the beginning of the row when necessary.

    - If `c_1 = c_2`, that is, the two symbols lie in the same column,
        each symbol is replaced by the symbol immediately below it,
        wrapping around to the top of the column when necessary.

    - Otherwise, the two symbols form the corners of a rectangle in the
        key square. Each symbol is replaced by the symbol in the same row
        but in the column occupied by the other symbol. Thus, if
        `(r_1,c_1)` and `(r_2,c_2)` are the coordinates of the plaintext
        symbols, the ciphertext symbols occupy positions `(r_1,c_2)` and
        `(r_2,c_1)`.

    To decrypt a ciphertext digraph, the inverse operations are
    applied:

    After decryption, the text will contain the filler symbols and all
    instances of J will still appear in place of I. Traditionally, these
    were manually removed or overlooked by the recipient of the message who
    could easily deduce the meaning of the text.

    EXAMPLES::

        sage: PF = PlayfairCryptosystem()
        sage: PF
        Playfair cryptosystem with alphabet 'ABCDEFGHIKLMNOPQRSTUVWXYZ' (I/J combined)
    """
    def __init__(self):
        r"""
        See ``PlayfairCryptosystem`` for full documentation.
        """
        self.alpha_no_j = 'ABCDEFGHIKLMNOPQRSTUVWXYZ'

    def __repr__(self):
        r"""
        Returns string representation of the cryptosystem.
        Defines the official string representation of an object
        """
        return "Playfair cryptosystem with alphabet 'ABCDEFGHIKLMNOPQRSTUVWXYZ' (I/J combined)"

    def format_word(self, word=''):
        r"""
        Normalize a word for Playfair processing. Will accept any string.
        If no string is provided, the empty string is used and returned.

        INPUT:

        - ``word`` -- input string.

        OUTPUT:

        Processed string where:
        - all letters are made uppercase
        - J replaced by I
        - non-alphabetic characters removed

        EXAMPLES::
        Capitalize all letters::

            sage: PF = PlayfairCryptosystem()
            sage: PF.format_word('hello')
            'HELLO'

            sage: PF = PlayfairCryptosystem()
            sage: PF.format_word('CaMeL')
            'CAMEL'

        'J' is substituted with 'I'::

            sage: PF = PlayfairCryptosystem()
            sage: PF.format_word('jello')
            'IELLO'

        Ignore whitespace::

            sage: PF = PlayfairCryptosystem()
            sage: PF.format_word('  C a M    e    L   ')
            'CAMEL'

        Remove all non-alphabetic characters::

            sage: PF = PlayfairCryptosystem()
            sage: PF.format_word('...')
            ''

            sage: PF = PlayfairCryptosystem()
            sage: PF.format_word('Do you, george, enjoy cheezeburgers? For lunch?')
            'DOYOUGEORGEENIOYCHEEZEBURGERSFORLUNCH'

            sage: PF = PlayfairCryptosystem()
            sage: PF.format_word('999A99')
            'A'
        """
        # Convert all letters to uppercase.
        word = word.upper()
        # Replace J with I for the Playfair alphabet.
        word = word.replace('J', 'I')
        # Build a new string containing only alphabetic characters.
        stripped_word = ''
        for char in word:
            if char.isalpha():
                stripped_word += char
        # Return the cleaned word.
        return stripped_word

    def create_key(self, keyword=''):
        r"""
        Construct a 5x5 Playfair key matrix

        The keyword will be formatted, so no errors will arise regardless
        of the formatting of the keyword string given. Will default to the
        empty string if no keyword is given.

        INPUT:

        - ``keyword`` -- keyword string.

        OUTPUT:

        A 5x5 list-of-lists representing the Playfair key matrix.

        EXAMPLES::
        If there is no keyword passed in, or the processed string is empty
        then the matrix is just filled in alphabetical order (skipping 'J')::

            sage: PF = PlayfairCryptosystem()
            sage: PF.create_key()
            [['A', 'B', 'C', 'D', 'E'],
             ['F', 'G', 'H', 'I', 'K'],
             ['L', 'M', 'N', 'O', 'P'],
             ['Q', 'R', 'S', 'T', 'U'],
             ['V', 'W', 'X', 'Y', 'Z']]

            sage: PF = PlayfairCryptosystem()
            sage: PF.create_key('')
            [['A', 'B', 'C', 'D', 'E'],
             ['F', 'G', 'H', 'I', 'K'],
             ['L', 'M', 'N', 'O', 'P'],
             ['Q', 'R', 'S', 'T', 'U'],
             ['V', 'W', 'X', 'Y', 'Z']]

            sage: PF = PlayfairCryptosystem()
            sage: PF.create_key('999888777666')
            [['A', 'B', 'C', 'D', 'E'],
             ['F', 'G', 'H', 'I', 'K'],
             ['L', 'M', 'N', 'O', 'P'],
             ['Q', 'R', 'S', 'T', 'U'],
             ['V', 'W', 'X', 'Y', 'Z']]

        The matrix is filled with the processed string first with no repeat
        letters, then the remainder of the alphabet put in::

             sage: PF = PlayfairCryptosystem()
             sage: PF.create_key('zyxwvut')
             [['Z', 'Y', 'X', 'W', 'V'],
              ['U', 'T', 'A', 'B', 'C'],
              ['D', 'E', 'F', 'G', 'H'],
              ['I', 'K', 'L', 'M', 'N'],
              ['O', 'P', 'Q', 'R', 'S']]

             sage: PF = PlayfairCryptosystem()
             sage: PF.create_key('jello wiggles')
             [['I', 'E', 'L', 'O', 'W'],
              ['G', 'S', 'A', 'B', 'C'],
              ['D', 'F', 'H', 'K', 'M'],
              ['N', 'P', 'Q', 'R', 'T'],
              ['U', 'V', 'X', 'Y', 'Z']]

            sage: PF = PlayfairCryptosystem()
            sage: PF.create_key('i AM the G.O.A.T.')
            [['I', 'A', 'M', 'T', 'H'],
             ['E', 'G', 'O', 'B', 'C'],
             ['D', 'F', 'K', 'L', 'N'],
             ['P', 'Q', 'R', 'S', 'U'],
             ['V', 'W', 'X', 'Y', 'Z']]
        """
        # Format the keyword by removing invalid characters
        # and replacing J with I.
        word = self.format_word(keyword)
        # Build a string of unique characters, starting with
        # the keyword and then the remaining alphabet letters.
        matrix_entries = ''
        for char in word + self.alpha_no_j:
            if char not in matrix_entries:
                matrix_entries += char
        # Arrange the characters into a 5x5 key matrix.
        key = [[matrix_entries[5*i + j] for j in range(5)] for i in range(5)]
        # Return the completed key matrix.
        return key

    def make_digraphs(self, text=''):
        r"""
        Convert text into Playfair digraphs.

        This method prepares the input text for encryption by transforming
        it into valid Playfair digraphs (pairs of letters). The process
        includes normalizing the text, handling repeated letters within
        pairs, and ensuring the final output has even length.

        If both letters in a digraph are the same, insert a filler character
        (by default, 'X') between them. Inserts 'Z' when the repeated letter
        is 'X' to avoid creating a double 'X pair

        INPUT:

        - ``text`` -- input string.

        OUTPUT:

        List of digraph strings (length 2).

        EXAMPLES::

            sage: PF = PlayfairCryptosystem()
            sage: PF.make_digraphs('i AM the G.O.A.T.')
            ['IA', 'MT', 'HE', 'GO', 'AT']

            sage: PF = PlayfairCryptosystem()
            sage: PF.make_digraphs('......')
            []

        A filler letter is inserted to prevent double letters::

            sage: PF = PlayfairCryptosystem()
            sage: PF.make_digraphs('ddubbel')
            ['DX', 'DU', 'BX', 'BE', 'LX']

            sage: PF = PlayfairCryptosystem()
            sage: PF.make_digraphs('XXXXX')
            ['XZ', 'XZ', 'XZ', 'XZ', 'XZ']

        A filler letter also inserted if there's an odd number of letters::

            sage: PF = PlayfairCryptosystem()
            sage: PF.make_digraphs('odd number of letters, why?')
            ['OD', 'DN', 'UM', 'BE', 'RO', 'FL', 'ET', 'TE', 'RS', 'WH', 'YX']
        """
        # Format the plaintext by converting to uppercase, removing
        # invalid characters, and replacing J with I as needed.
        tmp = self.format_word(text)
        # Process the string two characters at a time to create
        # valid Playfair digraphs.
        i = 0
        while 2*i + 1 < len(tmp):
            if tmp[2*i] == tmp[2*i + 1]:
                insert = 'Z' if tmp[2*i] == 'X' else 'X'
                tmp = tmp[:2*i + 1] + insert + tmp[2*i + 1:]
            else:
                i += 1
        # If the final string has odd length, append a filler
        # character so that every digraph contains exactly two letters.
        if len(tmp) % 2 == 1:
            tmp += 'Z' if tmp[-1] == 'X' else 'X'
        # Split the processed string into a list of digraphs.
        return [tmp[i:i+2] for i in range(0, len(tmp), 2)]

    def _single_digraph_encrypt_(self, di, key_matrix):
        r"""
        Encrypt a single digraph using Playfair rules.

        This method performs encryption on a single pair of letters
        (a digraph) using a 5x5 Playfair key matrix. It locates each
        character of the digraph within the matrix and applies the
        standard Playfair transformation rules:

        - If both letters are in the same row, each is replaced by the
          letter immediately to its right (wrapping around the row).
        - If both letters are in the same column, each is replaced by
          the letter immediately below it (wrapping around the column).
        - Otherwise, the letters form opposite corners of a rectangle,
          and each is replaced by the letter in the same row but the
          other letter’s column.

        INPUT:

        - ``di`` -- digraph string (two different capital letters as a string,
                    the letter 'J' will not appear)
        - ``key_matrix`` -- 5x5 key matrix

        OUTPUT:

        Encrypted digraph string.

        EXAMPLES::

            sage: PF = PlayfairCryptosystem()
            sage: key_mat = PF.create_key('zyxwvut')
            sage: PF._single_digraph_encrypt_('ZY', key_mat)
            'YX'

            sage: PF = PlayfairCryptosystem()
            sage: key_mat = PF.create_key('zyxwvut')
            sage: PF._single_digraph_encrypt_('XX', key_mat)
            'WW'

            sage: PF = PlayfairCryptosystem()
            sage: key_mat = PF.create_key('zyxwvut')
            sage: PF._single_digraph_encrypt_('AS', key_mat)
            'CQ'

            sage: PF = PlayfairCryptosystem()
            sage: key_mat = PF.create_key('zyxwvut')
            sage: PF._single_digraph_encrypt_('ZU', key_mat)
            'UD'
        """
        i0 = j0 = i1 = j1 = None
        for i in range(5):
            for j in range(5):
                if key_matrix[i][j] == di[0]:
                    i0, j0 = i, j
                if key_matrix[i][j] == di[1]:
                    i1, j1 = i, j
        # same row
        if i0 == i1:
            return key_matrix[i0][(j0+1)%5] + key_matrix[i1][(j1+1)%5]
        # same column
        if j0 == j1:
            return key_matrix[(i0+1)%5][j0] + key_matrix[(i1+1)%5][j1]
        # rectangle
        return key_matrix[i0][j1] + key_matrix[i1][j0]

    def enciphering(self, plaintext='', keyword=''):
        r"""
        Encrypt plaintext using Playfair cipher.

        INPUT:

        - ``plaintext`` -- plaintext string
        - ``keyword`` -- keyword

        OUTPUT:

        Ciphertext string.

        EXAMPLES::

            sage: PF = PlayfairCryptosystem()
            sage: PF.enciphering('encode this please using playfair', 'keyword')
            'DULCGDVFLQSHYDNOZNGQHNHCKHBHBT'

        If no keyword provided, will default to the empty string as keyword::

            sage: PF = PlayfairCryptosystem()
            sage: PF.enciphering('Do you... like secrets, jane?')
            'ITDTQPKFCUADUBUTFDPC'

        If no parameters provided, just returns the empty string::
        sage: PF = PlayfairCryptosystem()
        sage: PF.enciphering()
        ''
        """
        key_matrix = self.create_key(keyword)
        digraphs = self.make_digraphs(plaintext)
        return ''.join(self._single_digraph_encrypt_(d, key_matrix) for d in digraphs)

    def _single_digraph_decrypt_(self, di, key_matrix):
        r"""
        Decrypt a single digraph using Playfair cipher rules.

        This method reverses the Playfair encryption process for a
        single pair of letters (digraph) using a 5x5 key matrix.

        Each character of the digraph is located in the key matrix,
        and the inverse transformation rules are applied:

        - If both letters are in the same row, each is replaced by the
          letter immediately to its left (wrapping around the row).
        - If both letters are in the same column, each is replaced by
          the letter immediately above it (wrapping around the column).
        - Otherwise, the letters form opposite corners of a rectangle,
          and each is replaced by the letter in the same row but the
          other letter’s column.

        INPUT:

        - ``di`` -- two-character string (digraph to decrypt)
        - ``key_matrix`` -- 5x5 Playfair key matrix

        OUTPUT:

        Decrypted two-character digraph string.

        EXAMPLES::

            sage: PF = PlayfairCryptosystem()
            sage: key_mat = PF.create_key('zyxwvut')
            sage: PF._single_digraph_decrypt_('YX', key_mat)
            'ZY'

            sage: PF = PlayfairCryptosystem()
            sage: key_mat = PF.create_key('zyxwvut')
            sage: PF._single_digraph_decrypt_('WW', key_mat)
            'XX'

            sage: PF = PlayfairCryptosystem()
            sage: key_mat = PF.create_key('zyxwvut')
            sage: PF._single_digraph_decrypt_('CQ', key_mat)
            'AS'

            sage: PF = PlayfairCryptosystem()
            sage: key_mat = PF.create_key('zyxwvut')
            sage: PF._single_digraph_decrypt_('UD', key_mat)
            'ZU'
        """
        i0 = j0 = i1 = j1 = None
        # Locate positions of both characters in the key matrix
        for i in range(5):
            for j in range(5):
                if key_matrix[i][j] == di[0]:
                    i0, j0 = i, j
                if key_matrix[i][j] == di[1]:
                    i1, j1 = i, j

        # same row: shift left (inverse of encryption shift right)
        if i0 == i1:
            return key_matrix[i0][(j0-1)%5] + key_matrix[i1][(j1-1)%5]
        # same column: shift up (inverse of encryption shift down)
        if j0 == j1:
            return key_matrix[(i0-1)%5][j0] + key_matrix[(i1-1)%5][j1]
        # rectangle: swap columns (same rule as encryption)
        return key_matrix[i0][j1] + key_matrix[i1][j0]

    def deciphering(self, ciphertext='', keyword=''):
        r"""
        Decrypt ciphertext using Playfair cipher.

        INPUT:

        - ``ciphertext`` -- ciphertext string
        - ``keyword`` -- keyword used to construct key matrix

        OUTPUT:

        Decrypted plaintext string.

        EXAMPLES::

            sage: PF = PlayfairCryptosystem()
            sage: PF.deciphering('DULCGDVFLQSHYDNOZNGQHNHCKHBHBT', 'keyword')
            'ENCODETHISPLEASEUSINGPLAYFAIRX'

        If no keyword provided, will default to the empty string as keyword::

            sage: PF = PlayfairCryptosystem()
            sage: PF.deciphering('ITDTQPKFCUADUBUTFDPC')
            'DOYOULIKESECRETSIANE'

        If no parameters provided, just returns the empty string::

            sage: PF = PlayfairCryptosystem()
            sage: PF.deciphering()
            ''
        """
        # Build the same key matrix used for encryption
        key_matrix = self.create_key(keyword)
        # Split ciphertext into digraphs (pairs of letters)
        digraphs = self.make_digraphs(ciphertext)
        # Decrypt each digraph using Playfair rules
        return ''.join(self._single_digraph_decrypt_(d, key_matrix) for d in digraphs)
