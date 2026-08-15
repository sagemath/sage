# sage.doctest: needs sage.combinat
"""
Ciphers
"""

#*****************************************************************************
#       Copyright (C) 2007 David Kohel <kohel@maths.usyd.edu.au>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# Ciphers should inherit from morphisms (of sets).
# Specific cipher types will implement their functions in terms of the key

from sage.structure.element import Element


class Cipher(Element):
    """
    Cipher class
    """
    def __init__(self, parent, key):
        r"""
        Create a cipher.

        EXAMPLES::

            sage: S = ShiftCryptosystem(AlphabeticStrings())
            sage: E = S(3)
            sage: E.parent()
            Shift cryptosystem on Free alphabetic string monoid on A-Z
            sage: E.key()
            3
        """
        Element.__init__(self, parent)
        self._key = key

    def __eq__(self, right):
        r"""
        Return whether ``self`` and ``right`` are equal.

        Two ciphers are equal when they have the same type, the same
        parent, and the same key.

        EXAMPLES::

            sage: S = ShiftCryptosystem(AlphabeticStrings())
            sage: S(13) == S(13)
            True
            sage: S(13) == S(2)
            False
        """
        return type(self) is type(right) and self.parent() == right.parent() and self._key == right._key

    def _repr_(self):
        r"""
        Return the string representation of this cipher.

        EXAMPLES::

            sage: S = ShiftCryptosystem(AlphabeticStrings())
            sage: S(13)
            Shift cipher on Free alphabetic string monoid on A-Z
        """
        # return str(self._key)
        return "Cipher on %s" % self.parent().cipher_domain()

    def key(self):
        r"""
        Return the key of this cipher.

        EXAMPLES::

            sage: S = ShiftCryptosystem(HexadecimalStrings())
            sage: S(5).key()
            5
        """
        return self._key # was str(self._key)

    def domain(self):
        r"""
        Return the domain (plaintext space) of this cipher.

        EXAMPLES:

        For the ciphers shipped with Sage the plaintext and ciphertext
        alphabets coincide::

            sage: S = ShiftCryptosystem(AlphabeticStrings())
            sage: S(13).domain()
            Free alphabetic string monoid on A-Z

        They are distinct attributes of the parent, though, and need not
        agree if the cryptosystem was built with different alphabets::

            sage: from sage.crypto.cipher import Cipher
            sage: from sage.crypto.cryptosystem import SymmetricKeyCryptosystem
            sage: cs = SymmetricKeyCryptosystem(AlphabeticStrings(), HexadecimalStrings(), None)
            sage: Cipher(cs, None).domain()
            Free alphabetic string monoid on A-Z
        """
        return self.parent().cipher_domain()

    def codomain(self):
        r"""
        Return the codomain (ciphertext space) of this cipher.

        EXAMPLES:

        For the ciphers shipped with Sage the plaintext and ciphertext
        alphabets coincide::

            sage: S = ShiftCryptosystem(AlphabeticStrings())
            sage: S(13).codomain()
            Free alphabetic string monoid on A-Z

        The codomain reads the *ciphertext* alphabet, which may differ from
        the domain::

            sage: from sage.crypto.cipher import Cipher
            sage: from sage.crypto.cryptosystem import SymmetricKeyCryptosystem
            sage: cs = SymmetricKeyCryptosystem(AlphabeticStrings(), HexadecimalStrings(), None)
            sage: Cipher(cs, None).codomain()
            Free hexadecimal string monoid
        """
        return self.parent().cipher_codomain()


class SymmetricKeyCipher(Cipher):
    """
    Symmetric key cipher class
    """
    def __init__(self, parent, key):
        r"""
        Create a symmetric key cipher.

        EXAMPLES::

            sage: from sage.crypto.cipher import SymmetricKeyCipher
            sage: A = AffineCryptosystem(AlphabeticStrings())
            sage: E = A(3, 5)
            sage: E
            Affine cipher on Free alphabetic string monoid on A-Z
            sage: isinstance(E, SymmetricKeyCipher)
            True
        """
        Cipher.__init__(self, parent, key)


class PublicKeyCipher(Cipher):
    """
    Public key cipher class
    """
    def __init__(self, parent, key, public=True):
        """
        Create a public key cipher.
        """
        Cipher.__init__(self, parent, key)
        self._public = public
