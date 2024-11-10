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
        """
        Create a cipher.
        """
        Element.__init__(self, parent)
        self._key = key

    def __eq__(self, right):
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
        return self._key # was str(self._key)

    def domain(self):
        return self.parent().cipher_domain()

    def codomain(self):
        return self.parent().cipher_codomain()


class SymmetricKeyCipher(Cipher):
    """
    Symmetric key cipher class
    """
    def __init__(self, parent, key):
        """
        Create a symmetric cipher.
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
