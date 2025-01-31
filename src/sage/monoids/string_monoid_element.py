"""
String Monoid Elements

AUTHORS:

- David Kohel <kohel@maths.usyd.edu.au>, 2007-01

Elements of free string monoids, internal representation subject to change.

These are special classes of free monoid elements with distinct printing.

The internal representation of elements does not use the exponential
compression of FreeMonoid elements (a feature), and could be packed into words.
"""
# ****************************************************************************
#       Copyright (C) 2007 David Kohel <kohel@maths.usyd.edu.au>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.misc.lazy_import import lazy_import
from sage.rings.integer import Integer
from sage.structure.richcmp import richcmp

lazy_import('sage.rings.real_mpfr', 'RealField')

from .free_monoid_element import FreeMonoidElement


def is_StringMonoidElement(x):
    from sage.misc.superseded import deprecation
    deprecation(38280,
                "The function is_StringMonoidElement is deprecated; "
                "use 'isinstance(..., StringMonoidElement)' instead.")
    return isinstance(x, StringMonoidElement)


def is_AlphabeticStringMonoidElement(x):
    from sage.misc.superseded import deprecation
    deprecation(38280,
                "The function is_AlphabeticStringMonoidElement is deprecated; "
                "use 'isinstance(..., StringMonoidElement) and isinstance(x.parent(), AlphabeticStringMonoid)' instead.")
    from .string_monoid import AlphabeticStringMonoid
    return isinstance(x, StringMonoidElement) and \
        isinstance(x.parent(), AlphabeticStringMonoid)


def is_BinaryStringMonoidElement(x):
    from sage.misc.superseded import deprecation
    deprecation(38280,
                "The function is_BinaryStringMonoidElement is deprecated; "
                "use 'isinstance(..., StringMonoidElement) and isinstance(x.parent(), BinaryStringMonoid)' instead.")
    from .string_monoid import BinaryStringMonoid
    return isinstance(x, StringMonoidElement) and \
        isinstance(x.parent(), BinaryStringMonoid)


def is_OctalStringMonoidElement(x):
    from sage.misc.superseded import deprecation
    deprecation(38280,
                "The function is_OctalStringMonoidElement is deprecated; "
                "use 'isinstance(..., StringMonoidElement) and isinstance(x.parent(), OctalStringMonoid)' instead.")
    from .string_monoid import OctalStringMonoid
    return isinstance(x, StringMonoidElement) and \
        isinstance(x.parent(), OctalStringMonoid)


def is_HexadecimalStringMonoidElement(x):
    from sage.misc.superseded import deprecation
    deprecation(38280,
                "The function is_HexadecimalStringMonoidElement is deprecated; "
                "use 'isinstance(..., StringMonoidElement) and isinstance(x.parent(), HexadecimalStringMonoid)' instead.")
    from .string_monoid import HexadecimalStringMonoid
    return isinstance(x, StringMonoidElement) and \
        isinstance(x.parent(), HexadecimalStringMonoid)


def is_Radix64StringMonoidElement(x):
    from sage.misc.superseded import deprecation
    deprecation(38280,
                "The function is_Radix64StringMonoidElement is deprecated; "
                "use 'isinstance(..., StringMonoidElement) and isinstance(x.parent(), Radix64StringMonoid)' instead.")
    from .string_monoid import Radix64StringMonoid
    return isinstance(x, StringMonoidElement) and \
        isinstance(x.parent(), Radix64StringMonoid)


class StringMonoidElement(FreeMonoidElement):
    """
    Element of a free string monoid.
    """

    def __init__(self, S, x, check=True):
        """
        Create the element ``x`` of the StringMonoid ``S``.

        This should typically be called by a StringMonoid.
        """
        FreeMonoidElement.__init__(self, S, [])
        if isinstance(x, list):
            if check:
                for b in x:
                    if not isinstance(b, (int, Integer)):
                        raise TypeError(
                            "x (= %s) must be a list of integers." % x)
            self._element_list = list(x)  # make copy
        elif isinstance(x, str):
            alphabet = list(self.parent().alphabet())
            self._element_list = []
            for i in range(len(x)):
                try:
                    b = alphabet.index(x[i])
                except ValueError:
                    raise TypeError(
                        "Argument x (= %s) is not a valid string." % x)
                self._element_list += [b]
        else:
            raise TypeError("Argument x (= %s) is of the wrong type." % x)

    def _richcmp_(self, other, op):
        """
        Compare two free monoid elements with the same parents.

        The ordering is the one on the underlying sorted list of
        (monomial, coefficients) pairs.

        EXAMPLES::

            sage: S = BinaryStrings()
            sage: (x,y) = S.gens()
            sage: x * y < y * x
            True
            sage: S("01") < S("10")
            True
        """
        return richcmp(self._element_list, other._element_list, op)

    def _repr_(self):
        """
        The self-representation of a string monoid element. Unlike Python
        strings we print without the enclosing quotes.
        """
        S = self.parent()
        alphabet = S.alphabet()
        s = ''
        for b in self._element_list:
            c = alphabet[b]
            s += c
        return s

    def _latex_(self):
        """
        Return latex representation of ``self``.

        EXAMPLES::

            sage: S = BinaryStrings()
            sage: s = S('101111000')
            sage: latex(s)
            101111000
        """
        return self._repr_()

    def __mul__(self, y):
        """
        Multiply 2 free string monoid elements.

        EXAMPLES::

            sage: S = BinaryStrings()
            sage: (x,y) = S.gens()
            sage: x*y
            01
        """
        if not isinstance(y, StringMonoidElement):
            raise TypeError("Argument y (= %s) is of wrong type." % y)
        S = self.parent()
        x_elt = self._element_list
        y_elt = y._element_list
        z = S('')
        z._element_list = x_elt + y_elt
        return z

    def __pow__(self, n):
        """
        Return the `n`-th power of the string element.

        EXAMPLES::

            sage: (x,y) = BinaryStrings().gens()
            sage: x**3 * y**5 * x**7
            000111110000000
            sage: x**0


        Note that raising to a negative power is *not* a constructor
        for an element of the corresponding free group (yet).

        ::

            sage: x**(-1)
            Traceback (most recent call last):
            ...
            IndexError: Argument n (= -1) must be nonnegative.
        """
        if not isinstance(n, (int, Integer)):
            raise TypeError("Argument n (= %s) must be an integer." % n)
        if n < 0:
            raise IndexError("Argument n (= %s) must be nonnegative." % n)
        elif n == 0:
            return self.parent()('')
        elif n == 1:
            return self
        z = self.parent()('')
        z._element_list = self._element_list * n
        return z

    def __len__(self):
        """
        Return the number of products that occur in this monoid element.
        For example, the length of the identity is 0, and the length
        of the monoid `x_0^2x_1` is three.

        EXAMPLES::

            sage: S = BinaryStrings()
            sage: z = S('')
            sage: len(z)
            0
            sage: (x,y) = S.gens()
            sage: len(x**2 * y**3)
            5
        """
        return len(self._element_list)

    def __iter__(self):
        """
        Return an iterator over this element as a string.

        EXAMPLES::

            sage: t = AlphabeticStrings()('SHRUBBERY')
            sage: next(t.__iter__())
            S
            sage: list(t)
            [S, H, R, U, B, B, E, R, Y]
        """
        l = len(self._element_list)
        i = 0
        while i < l:
            yield self[i]
            i += 1

    def __getitem__(self, n):
        """
        Return the n-th string character.

        EXAMPLES::

            sage: t = AlphabeticStrings()('SHRUBBERY')
            sage: t[0]
            S
            sage: t[3]
            U
            sage: t[-1]
            Y
        """
        try:
            c = self._element_list[n]
        except Exception:
            raise IndexError("Argument n (= %s) is not a valid index." % n)
        if not isinstance(c, list):
            c = [c]
        return self.parent()(c)

    def decoding(self, padic=False):
        r"""
        The byte string associated to a binary or hexadecimal string
        monoid element.

        EXAMPLES::

            sage: S = HexadecimalStrings()
            sage: s = S.encoding("A..Za..z"); s
            412e2e5a612e2e7a
            sage: s.decoding()
            'A..Za..z'
            sage: s = S.encoding("A..Za..z",padic=True); s
            14e2e2a516e2e2a7
            sage: s.decoding()
            '\x14\xe2\xe2\xa5\x16\xe2\xe2\xa7'
            sage: s.decoding(padic=True)
            'A..Za..z'
            sage: S = BinaryStrings()
            sage: s = S.encoding("A..Za..z"); s
            0100000100101110001011100101101001100001001011100010111001111010
            sage: s.decoding()
            'A..Za..z'
            sage: s = S.encoding("A..Za..z",padic=True); s
            1000001001110100011101000101101010000110011101000111010001011110
            sage: s.decoding()
            '\x82ttZ\x86tt^'
            sage: s.decoding(padic=True)
            'A..Za..z'
        """
        S = self.parent()
        from .string_monoid import (AlphabeticStringMonoid,
                                    BinaryStringMonoid,
                                    HexadecimalStringMonoid)
        if isinstance(S, AlphabeticStringMonoid):
            return ''.join(chr(65 + i) for i in self._element_list)
        n = len(self)
        if isinstance(S, HexadecimalStringMonoid):
            if n % 2:
                "String %s must have even length to determine a byte character string." % str(self)
            s = []
            x = self._element_list
            for k in range(n//2):
                m = 2*k
                if padic:
                    c = chr(x[m]+16*x[m+1])
                else:
                    c = chr(16*x[m]+x[m+1])
                s.append(c)
            return ''.join(s)
        if isinstance(S, BinaryStringMonoid):
            if not n % 8 == 0:
                "String %s must have even length 0 mod 8 to determine a byte character string." % str(self)
            pows = [2**i for i in range(8)]
            s = []
            x = self._element_list
            for k in range(n//8):
                m = 8*k
                if padic:
                    c = chr(sum([x[m+i] * pows[i] for i in range(8)]))
                else:
                    c = chr(sum([x[m+7-i] * pows[i] for i in range(8)]))
                s.append(c)
            return ''.join(s)
        raise TypeError(
            "Argument %s must be an alphabetic, binary, or hexadecimal string." % str(self))

    def coincidence_index(self, prec=0):
        """
        Return the probability of two randomly chosen characters being equal.
        """
        if prec == 0:
            RR = RealField()
        else:
            RR = RealField(prec)
        char_dict = {}
        for i in self._element_list:
            if i in char_dict:
                char_dict[i] += 1
            else:
                char_dict[i] = 1
        nn = 0
        ci_num = 0
        for ni in char_dict.values():
            nn += ni
            ci_num += ni * (ni - 1)
        ci_den = nn * (nn - 1)
        return RR(ci_num) / ci_den

    def character_count(self):
        r"""
        Return the count of each unique character.

        EXAMPLES:

        Count the character frequency in an object comprised of capital
        letters of the English alphabet::

            sage: M = AlphabeticStrings().encoding("abcabf")
            sage: sorted(M.character_count().items())
            [(A, 2), (B, 2), (C, 1), (F, 1)]

        In an object comprised of binary numbers::

            sage: M = BinaryStrings().encoding("abcabf")
            sage: sorted(M.character_count().items())
            [(0, 28), (1, 20)]

        In an object comprised of octal numbers::

            sage: A = OctalStrings()
            sage: M = A([1, 2, 3, 2, 5, 3])
            sage: sorted(M.character_count().items())
            [(1, 1), (2, 2), (3, 2), (5, 1)]

        In an object comprised of hexadecimal numbers::

            sage: A = HexadecimalStrings()
            sage: M = A([1, 2, 4, 6, 2, 4, 15])
            sage: sorted(M.character_count().items())
            [(1, 1), (2, 2), (4, 2), (6, 1), (f, 1)]

        In an object comprised of radix-64 characters::

            sage: A = Radix64Strings()
            sage: M = A([1, 2, 63, 45, 45, 10]); M
            BC/ttK
            sage: sorted(M.character_count().items())
            [(B, 1), (C, 1), (K, 1), (t, 2), (/, 1)]

        TESTS:

        Empty strings return no counts of character frequency::

            sage: M = AlphabeticStrings().encoding("")
            sage: M.character_count()
            {}
            sage: M = BinaryStrings().encoding("")
            sage: M.character_count()
            {}
            sage: A = OctalStrings()
            sage: M = A([])
            sage: M.character_count()
            {}
            sage: A = HexadecimalStrings()
            sage: M = A([])
            sage: M.character_count()
            {}
            sage: A = Radix64Strings()
            sage: M = A([])
            sage: M.character_count()
            {}
        """
        # the character frequency, i.e. the character count
        CF = {}
        for e in self:
            if e in CF:
                CF[e] += 1
            else:
                CF.setdefault(e, 1)
        return CF

    def frequency_distribution(self, length=1, prec=0):
        """
        Return the probability space of character frequencies.

        The output of this method is different from that of the method
        :func:`characteristic_frequency()
        <sage.monoids.string_monoid.AlphabeticStringMonoid.characteristic_frequency>`.

        One can think of the characteristic frequency probability of an
        element in an alphabet `A` as the expected probability of that element
        occurring. Let `S` be a string encoded using elements of `A`. The
        frequency probability distribution corresponding to `S` provides us
        with the frequency probability of each element of `A` as observed
        occurring in `S`. Thus one distribution provides expected
        probabilities, while the other provides observed probabilities.

        INPUT:

        - ``length`` -- (default: ``1``) if ``length=1`` then consider the
          probability space of monogram frequency, i.e. probability
          distribution of single characters. If ``length=2`` then consider
          the probability space of digram frequency, i.e. probability
          distribution of pairs of characters. This method currently
          supports the generation of probability spaces for monogram
          frequency (``length=1``) and digram frequency (``length=2``).

        - ``prec`` -- (default: ``0``) a nonnegative integer representing
          the precision (in number of bits) of a floating-point number. The
          default value ``prec=0`` means that we use 53 bits to represent
          the mantissa of a floating-point number. For more information on
          the precision of floating-point numbers, see the function
          :func:`RealField() <sage.rings.real_mpfr.RealField>` or refer to the module
          :mod:`real_mpfr <sage.rings.real_mpfr>`.

        EXAMPLES:

        Capital letters of the English alphabet::

            sage: M = AlphabeticStrings().encoding("abcd")
            sage: L = M.frequency_distribution().function()
            sage: sorted(L.items())
            <BLANKLINE>
            [(A, 0.250000000000000),
            (B, 0.250000000000000),
            (C, 0.250000000000000),
            (D, 0.250000000000000)]

        The binary number system::

            sage: M = BinaryStrings().encoding("abcd")
            sage: L = M.frequency_distribution().function()
            sage: sorted(L.items())
            [(0, 0.593750000000000), (1, 0.406250000000000)]

        The hexadecimal number system::

            sage: M = HexadecimalStrings().encoding("abcd")
            sage: L = M.frequency_distribution().function()
            sage: sorted(L.items())
            <BLANKLINE>
            [(1, 0.125000000000000),
            (2, 0.125000000000000),
            (3, 0.125000000000000),
            (4, 0.125000000000000),
            (6, 0.500000000000000)]

        Get the observed frequency probability distribution of digrams in the
        string "ABCD". This string consists of the following digrams: "AB",
        "BC", and "CD". Now find out the frequency probability of each of
        these digrams as they occur in the string "ABCD"::

            sage: M = AlphabeticStrings().encoding("abcd")
            sage: D = M.frequency_distribution(length=2).function()
            sage: sorted(D.items())
            [(AB, 0.333333333333333), (BC, 0.333333333333333), (CD, 0.333333333333333)]
        """
        from sage.probability.random_variable import DiscreteProbabilitySpace
        if length not in (1, 2):
            raise NotImplementedError("Not implemented")
        if prec == 0:
            RR = RealField()
        else:
            RR = RealField(prec)
        S = self.parent()
        if length == 1:
            Alph = S.gens()
        else:
            Alph = tuple(x * y for x in S.gens() for y in S.gens())
        X = {}
        N = len(self) - length + 1
        eps = RR(Integer(1) / N)
        for i in range(N):
            c = self[i:i+length]
            if c in X:
                X[c] += eps
            else:
                X[c] = eps
        # Return a dictionary of probability distribution. This should
        # allow for easier parsing of the dictionary.
        return DiscreteProbabilitySpace(Alph, X, RR)
