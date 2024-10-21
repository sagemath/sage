# sage_setup: distribution = sagemath-modules
# sage.doctest: needs sage.combinat
r"""
Cryptosystems

This module contains base classes for various cryptosystems, including
symmetric key and public-key cryptosystems. The classes defined in this
module should not be called directly. It is the responsibility of child
classes to implement specific cryptosystems.

Take for example the
Hill or matrix cryptosystem as implemented in
:class:`HillCryptosystem <sage.crypto.classical.HillCryptosystem>`. It is a
symmetric key cipher so
:class:`HillCryptosystem <sage.crypto.classical.HillCryptosystem>` is a child
class of
:class:`SymmetricKeyCryptosystem <sage.crypto.cryptosystem.SymmetricKeyCryptosystem>`,
which in turn is a child class of
:class:`Cryptosystem <sage.crypto.cryptosystem.Cryptosystem>`. The following
diagram shows the inheritance relationship of particular cryptosystems::

    Cryptosystem
    + SymmetricKeyCryptosystem
    | + HillCryptosystem
    | + LFSRCryptosystem
    | + ShiftCryptosystem
    | + ShrinkingGeneratorCryptosystem
    | + SubstitutionCryptosystem
    | + TranspositionCryptosystem
    | + VigenereCryptosystem
    + PublicKeyCryptosystem
"""
# ****************************************************************************
#       Copyright (C) 2007 David Kohel <kohel@maths.usyd.edu.au>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.sets.set import Set_generic


class Cryptosystem(Set_generic):
    r"""
    A base cryptosystem class. This is meant to be extended by other
    specialized child classes that implement specific cryptosystems.

    A cryptosystem is a pair of maps

    .. MATH::

        E : {\mathcal K} \rightarrow {\rm Hom}({\mathcal M},{\mathcal C})

    .. MATH::

        D : {\mathcal K} \rightarrow {\rm Hom}({\mathcal C},{\mathcal M})

    where `{\mathcal K}` is the key space,
    `{\mathcal M}` is the plaintext or message space, and
    `{\mathcal C}` is the ciphertext space. In many instances
    `{\mathcal M} = {\mathcal C}` and the images will lie in
    `{\rm Aut}({\mathcal M})`. An element of the image of
    `E` is called a cipher.

    We may assume that `E` and `D` are injective, hence
    identify a key `K` in `{\mathcal K}` with its image
    `E_K := E(K)` in
    `\mathrm{Hom}({\mathcal M},{\mathcal C})`.

    The cryptosystem has the property that for every encryption key
    `K_1` there is a decryption key `K_2` such that
    `D_{K_2} \circ E_{K_1}`. A cryptosystem with the
    property that `K := K_2 = K_1`, is called a symmetric
    cryptosystem. Otherwise, if the key `K_2 \ne K_1`, nor is
    `K_2` easily derived from `K_1`, we call the
    cryptosystem asymmetric or public key. In that case, `K_1`
    is called the public key and `K_2` is called the private
    key.

    INPUT:

    - ``plaintext_space`` -- the plaintext alphabet

    - ``ciphertext_space`` -- the ciphertext alphabet

    - ``key_space`` -- the key alphabet

    - ``block_length`` -- (default: 1) the block length

    - ``period`` -- (default: ``None``) the period

    EXAMPLES:

    Various classical cryptosystems::

        sage: ShiftCryptosystem(AlphabeticStrings())
        Shift cryptosystem on Free alphabetic string monoid on A-Z
        sage: SubstitutionCryptosystem(HexadecimalStrings())
        Substitution cryptosystem on Free hexadecimal string monoid
        sage: HillCryptosystem(BinaryStrings(), 3)
        Hill cryptosystem on Free binary string monoid of block length 3
        sage: TranspositionCryptosystem(OctalStrings(), 5)                              # needs sage.groups
        Transposition cryptosystem on Free octal string monoid of block length 5
        sage: VigenereCryptosystem(Radix64Strings(), 7)
        Vigenere cryptosystem on Free radix 64 string monoid of period 7
    """
    def __init__(self, plaintext_space, ciphertext_space, key_space,
                 block_length=1, period=None):
        r"""
        Create a ``Cryptosystem`` object. See the class ``Cryptosystem``
        for detailed documentation.

        INPUT:

        - ``plaintext_space`` -- the plaintext alphabet

        - ``ciphertext_space`` -- the ciphertext alphabet

        - ``key_space`` -- the key alphabet

        - ``block_length`` -- (default: 1) the block length

        - ``period`` -- (default: ``None``) the period

        EXAMPLES:

        Various classical cryptosystems::

            sage: ShiftCryptosystem(AlphabeticStrings())
            Shift cryptosystem on Free alphabetic string monoid on A-Z
            sage: SubstitutionCryptosystem(HexadecimalStrings())
            Substitution cryptosystem on Free hexadecimal string monoid
            sage: HillCryptosystem(BinaryStrings(), 3)
            Hill cryptosystem on Free binary string monoid of block length 3
            sage: TranspositionCryptosystem(OctalStrings(), 5)                          # needs sage.groups
            Transposition cryptosystem on Free octal string monoid of block length 5
            sage: VigenereCryptosystem(Radix64Strings(), 7)
            Vigenere cryptosystem on Free radix 64 string monoid of period 7
        """
        self._cipher_domain = plaintext_space
        self._cipher_codomain = ciphertext_space
        self._key_space = key_space
        self._block_length = block_length
        self._period = period

    def __eq__(self, right):
        r"""
        Comparing ``self`` with ``right``.

        Two ``Cryptosystem`` objects
        are the same if they satisfy all of these conditions:

        - share the same type
        - have the same cipher domain
        - have the same cipher codomain
        - share the same key space
        - share the same block length
        - have the same period

        INPUT:

        - ``right`` -- a ``Cryptosystem`` object

        EXAMPLES:

        Pairs of equivalent classical cryptosystems::

            sage: sub1 = SubstitutionCryptosystem(AlphabeticStrings())
            sage: sub2 = SubstitutionCryptosystem(AlphabeticStrings())
            sage: sub1 == sub2
            True
            sage: shift1 = ShiftCryptosystem(HexadecimalStrings())
            sage: shift2 = ShiftCryptosystem(HexadecimalStrings())
            sage: shift1 == shift2
            True
            sage: hill1 = HillCryptosystem(AlphabeticStrings(), 4)
            sage: hill2 = HillCryptosystem(AlphabeticStrings(), 4)
            sage: hill1 == hill2
            True
            sage: tran1 = TranspositionCryptosystem(HexadecimalStrings(), 5)            # needs sage.groups
            sage: tran2 = TranspositionCryptosystem(HexadecimalStrings(), 5)            # needs sage.groups
            sage: tran1 == tran2                                                        # needs sage.groups
            True
            sage: vig1 = VigenereCryptosystem(AlphabeticStrings(), 7)
            sage: vig2 = VigenereCryptosystem(AlphabeticStrings(), 7)
            sage: vig1 == vig2
            True

        Pairs of different classical cryptosystems::

            sage: sub1 = SubstitutionCryptosystem(AlphabeticStrings())
            sage: sub2 = SubstitutionCryptosystem(OctalStrings())
            sage: sub1 == sub2
            False
            sage: shift1 = ShiftCryptosystem(HexadecimalStrings())
            sage: shift2 = ShiftCryptosystem(BinaryStrings())
            sage: shift1 == shift2
            False
            sage: hill1 = HillCryptosystem(Radix64Strings(), 4)
            sage: hill2 = HillCryptosystem(Radix64Strings(), 5)
            sage: hill1 == hill2
            False
            sage: tran1 = TranspositionCryptosystem(Radix64Strings(), 3)                # needs sage.groups
            sage: tran2 = TranspositionCryptosystem(HexadecimalStrings(), 3)            # needs sage.groups
            sage: tran1 == tran2                                                        # needs sage.groups
            False
            sage: vig1 = VigenereCryptosystem(AlphabeticStrings(), 7)
            sage: vig2 = VigenereCryptosystem(Radix64Strings(), 7)
            sage: vig1 == vig2
            False
        """
        return (type(self) is type(right) and
            self._cipher_domain == right._cipher_domain and
            self._cipher_codomain == right._cipher_codomain and
            self._key_space == right._key_space and
            self._block_length == right._block_length and
            self._period == right._period)

    def plaintext_space(self):
        r"""
        Return the plaintext alphabet of this cryptosystem.

        EXAMPLES:

        The plaintext spaces of various classical cryptosystems::

            sage: ShiftCryptosystem(AlphabeticStrings()).plaintext_space()
            Free alphabetic string monoid on A-Z
            sage: SubstitutionCryptosystem(HexadecimalStrings()).plaintext_space()
            Free hexadecimal string monoid
            sage: HillCryptosystem(BinaryStrings(), 3).plaintext_space()
            Free binary string monoid
            sage: TranspositionCryptosystem(OctalStrings(), 5).plaintext_space()        # needs sage.groups
            Free octal string monoid
            sage: VigenereCryptosystem(Radix64Strings(), 7).plaintext_space()
            Free radix 64 string monoid
        """
        return self._cipher_domain

    def cipher_domain(self):
        r"""
        Return the alphabet used by this cryptosystem for encoding plaintexts.
        This is the same as the plaintext space.

        EXAMPLES:

        The cipher domains, or plaintext spaces, of various classical
        cryptosystems::

            sage: ShiftCryptosystem(AlphabeticStrings()).cipher_domain()
            Free alphabetic string monoid on A-Z
            sage: SubstitutionCryptosystem(HexadecimalStrings()).cipher_domain()
            Free hexadecimal string monoid
            sage: HillCryptosystem(BinaryStrings(), 3).cipher_domain()
            Free binary string monoid
            sage: TranspositionCryptosystem(OctalStrings(), 5).cipher_domain()          # needs sage.groups
            Free octal string monoid
            sage: VigenereCryptosystem(Radix64Strings(), 7).cipher_domain()
            Free radix 64 string monoid
        """
        return self._cipher_domain

    def ciphertext_space(self):
        r"""
        Return the ciphertext alphabet of this cryptosystem.

        EXAMPLES:

        The ciphertext spaces of various classical cryptosystems::

            sage: ShiftCryptosystem(AlphabeticStrings()).ciphertext_space()
            Free alphabetic string monoid on A-Z
            sage: SubstitutionCryptosystem(HexadecimalStrings()).ciphertext_space()
            Free hexadecimal string monoid
            sage: HillCryptosystem(BinaryStrings(), 3).ciphertext_space()
            Free binary string monoid
            sage: TranspositionCryptosystem(OctalStrings(), 5).ciphertext_space()       # needs sage.groups
            Free octal string monoid
            sage: VigenereCryptosystem(Radix64Strings(), 7).ciphertext_space()
            Free radix 64 string monoid
        """
        return self._cipher_codomain

    def cipher_codomain(self):
        r"""
        Return the alphabet used by this cryptosystem for encoding ciphertexts.
        This is the same as the ciphertext space.

        EXAMPLES:

        The cipher codomains, or ciphertext spaces, of various classical
        cryptosystems::

            sage: ShiftCryptosystem(AlphabeticStrings()).cipher_codomain()
            Free alphabetic string monoid on A-Z
            sage: SubstitutionCryptosystem(HexadecimalStrings()).cipher_codomain()
            Free hexadecimal string monoid
            sage: HillCryptosystem(BinaryStrings(), 3).cipher_codomain()
            Free binary string monoid
            sage: TranspositionCryptosystem(OctalStrings(), 5).cipher_codomain()        # needs sage.groups
            Free octal string monoid
            sage: VigenereCryptosystem(Radix64Strings(), 7).cipher_codomain()
            Free radix 64 string monoid
        """
        return self._cipher_codomain

    def key_space(self):
        r"""
        Return the alphabet used by this cryptosystem for encoding keys.

        EXAMPLES:

        The key spaces of various classical cryptosystems::

            sage: ShiftCryptosystem(AlphabeticStrings()).key_space()
            Ring of integers modulo 26
            sage: SubstitutionCryptosystem(HexadecimalStrings()).key_space()
            Free hexadecimal string monoid
            sage: HillCryptosystem(BinaryStrings(), 3).key_space()
            Full MatrixSpace of 3 by 3 dense matrices over Ring of integers modulo 2
            sage: TranspositionCryptosystem(OctalStrings(), 5).key_space()              # needs sage.groups
            Symmetric group of order 5! as a permutation group
            sage: VigenereCryptosystem(Radix64Strings(), 7).key_space()
            Free radix 64 string monoid
        """
        return self._key_space

    def block_length(self):
        r"""
        Return the block length of this cryptosystem.

        For some cryptosystems
        this is not relevant, in which case the block length defaults to 1.

        EXAMPLES:

        The block lengths of various classical cryptosystems::

            sage: ShiftCryptosystem(AlphabeticStrings()).block_length()
            1
            sage: SubstitutionCryptosystem(HexadecimalStrings()).block_length()
            1
            sage: HillCryptosystem(BinaryStrings(), 3).block_length()
            3
            sage: TranspositionCryptosystem(OctalStrings(), 5).block_length()           # needs sage.groups
            5
            sage: VigenereCryptosystem(Radix64Strings(), 7).block_length()
            1
        """
        return self._block_length

    def period(self):
        if self._period is None:
            raise TypeError("Argument has no associated period.")
        return self._period


class SymmetricKeyCryptosystem(Cryptosystem):
    r"""
    The base class for symmetric key, or secret key, cryptosystems.
    """
    def alphabet_size(self):
        r"""
        Return the number of elements in the alphabet of this
        cryptosystem. This only applies to any cryptosystem whose plaintext
        and ciphertext spaces are the same alphabet.

        EXAMPLES::

            sage: ShiftCryptosystem(AlphabeticStrings()).alphabet_size()
            26
            sage: ShiftCryptosystem(BinaryStrings()).alphabet_size()
            2
            sage: ShiftCryptosystem(HexadecimalStrings()).alphabet_size()
            16
            sage: SubstitutionCryptosystem(OctalStrings()).alphabet_size()
            8
            sage: SubstitutionCryptosystem(Radix64Strings()).alphabet_size()
            64
        """
        return self._cipher_domain.ngens()


class PublicKeyCryptosystem(Cryptosystem):
    r"""
    The base class for asymmetric or public-key cryptosystems.
    """
