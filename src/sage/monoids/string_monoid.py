r"""
Free String Monoids

AUTHORS:

- David Kohel <kohel@maths.usyd.edu.au>, 2007-01

Sage supports a wide range of specific free string monoids.
"""
#*****************************************************************************
#       Copyright (C) 2007 David Kohel <kohel@maths.usyd.edu.au>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# from sage.rings.integer import Integer
# from sage.structure.parent_gens import ParentWithGens, normalize_names
from free_monoid import FreeMonoid_class
from string_monoid_element import StringMonoidElement
from string_ops import strip_encoding

import weakref

_cache = {}

def BinaryStrings():
    r"""
    Returns the free binary string monoid on generators `\{ 0, 1 \}`.

    OUTPUT:

    - Free binary string monoid.

    EXAMPLES::

        sage: S = BinaryStrings(); S
        Free binary string monoid
        sage: u = S('')
        sage: u

        sage: x = S('0')
        sage: x
        0
        sage: y = S('1')
        sage: y
        1
        sage: z = S('01110')
        sage: z
        01110
        sage: x*y^3*x == z
        True
        sage: u*x == x*u
        True
    """
    # Here we cache the binary strings to make them unique
    if 2 in _cache:
        S = _cache[2]()
        if not S is None:
            return S
    S = BinaryStringMonoid()
    _cache[2] = weakref.ref(S)
    return S

def OctalStrings():
    r"""
    Returns the free octal string monoid on generators `\{ 0, 1, \dots, 7 \}`.

    OUTPUT:

    - Free octal string monoid.

    EXAMPLES::

        sage: S = OctalStrings(); S
        Free octal string monoid
        sage: x = S.gens()
        sage: x[0]
        0
        sage: x[7]
        7
        sage: x[0] * x[3]^3 * x[5]^4 * x[6]
        033355556
    """
    # Here we cache the octal strings to make them unique
    if 8 in _cache:
        S = _cache[8]()
        if not S is None:
            return S
    S = OctalStringMonoid()
    _cache[8] = weakref.ref(S)
    return S

def HexadecimalStrings():
    r"""
    Returns the free hexadecimal string monoid on generators
    `\{ 0, 1, \dots , 9, a, b, c, d, e, f \}`.

    OUTPUT:

    - Free hexadecimal string monoid.

    EXAMPLES::

        sage: S = HexadecimalStrings(); S
        Free hexadecimal string monoid
        sage: x = S.gen(0)
        sage: y = S.gen(10)
        sage: z = S.gen(15)
        sage: z
        f
        sage: x*y^3*z
        0aaaf
    """
    # Here we cache the hexadecimal strings to make them unique
    if 16 in _cache:
        S = _cache[16]()
        if not S is None:
            return S
    S = HexadecimalStringMonoid()
    _cache[16] = weakref.ref(S)
    return S

def Radix64Strings():
    r"""
    Returns the free radix 64 string monoid on 64 generators

    ::

        A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z,
        a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,
        0,1,2,3,4,5,6,7,8,9,+,/

    OUTPUT:

    - Free radix 64 string monoid.

    EXAMPLES::

        sage: S = Radix64Strings(); S
        Free radix 64 string monoid
        sage: x = S.gens()
        sage: x[0]
        A
        sage: x[62]
        +
        sage: x[63]
        /
    """
    # Here we cache the radix-64 strings to make them unique
    if 64 in _cache:
        S = _cache[64]()
        if not S is None:
            return S
    S = Radix64StringMonoid()
    _cache[64] = weakref.ref(S)
    return S

def AlphabeticStrings():
    r"""
    Returns the string monoid on generators A-Z:
    `\{ A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z \}`.

    OUTPUT:

    - Free alphabetic string monoid on A-Z.

    EXAMPLES::

        sage: S = AlphabeticStrings(); S
        Free alphabetic string monoid on A-Z
        sage: x = S.gens()
        sage: x[0]
        A
        sage: x[25]
        Z
    """
    # Here we cache the alphabetic strings to make them unique
    if 26 in _cache:
        S = _cache[26]()
        if not S is None:
            return S
    S = AlphabeticStringMonoid()
    _cache[26] = weakref.ref(S)
    return S

#*****************************************************************************

class StringMonoid_class(FreeMonoid_class):
    r"""
    A free string monoid on `n` generators.
    """

    def __init__(self, n, alphabet=()):
        r"""
        Create free binary string monoid on `n` generators.

        INPUT:

        - ``n`` -- Integer

        - ``alphabet`` -- String or tuple whose characters or elements denote
          the generators.

        EXAMPLES::

            sage: S = BinaryStrings(); S
            Free binary string monoid
            sage: x = S.gens()
            sage: x[0]*x[1]**5 * (x[0]*x[1])
            01111101
        """
        # Names must be alphabetical -- omitted since printing is
        # defined locally.
        # FreeMonoid_class.__init__(self, n, names = alphabet)
        FreeMonoid_class.__init__(self, n)
        self._alphabet = alphabet

    def __contains__(self, x):
        return isinstance(x, StringMonoidElement) and x.parent() == self

    def alphabet(self):
        return tuple(self._alphabet)

    def gen(self, i=0):
        r"""
        The `i`-th generator of the monoid.

        INPUT:

        - ``i`` -- integer (default: 0)

        EXAMPLES::

            sage: S = BinaryStrings()
            sage: S.gen(0)
            0
            sage: S.gen(1)
            1
            sage: S.gen(2)
            Traceback (most recent call last):
            ...
            IndexError: Argument i (= 2) must be between 0 and 1.
            sage: S = HexadecimalStrings()
            sage: S.gen(0)
            0
            sage: S.gen(12)
            c
            sage: S.gen(16)
            Traceback (most recent call last):
            ...
            IndexError: Argument i (= 16) must be between 0 and 15.
        """
        n = self.ngens()
        if i < 0 or not i < n:
            raise IndexError(
                "Argument i (= %s) must be between 0 and %s." % (i, n-1))
        return StringMonoidElement(self, [int(i)])

#*****************************************************************************
# Specific global string monoids
#*****************************************************************************

class BinaryStringMonoid(StringMonoid_class):
    r"""
    The free binary string monoid on generators `\{ 0, 1 \}`.
    """

    def __init__(self):
        r"""
        Create free binary string monoid on generators `\{ 0, 1 \}`.

        EXAMPLES::

            sage: S = BinaryStrings(); S
            Free binary string monoid
            sage: x = S.gens()
            sage: x[0]*x[1]**5 * (x[0]*x[1])
            01111101
        """
        StringMonoid_class.__init__(self, 2, ['0', '1'])

    def __cmp__(self, other):
        if not isinstance(other, BinaryStringMonoid):
            return -1
        return 0

    def __repr__(self):
        return "Free binary string monoid"

    def __call__(self, x, check=True):
        r"""
        Return ``x`` coerced into this free monoid.

        One can create a free binary string monoid element from a
        Python string of 0's and 1's or list of integers.

        NOTE: Due to the ambiguity of the second generator '1' with
        the identity element '' of the monoid, the syntax S(1) is not
        permissible.

        EXAMPLES::

            sage: S = BinaryStrings()
            sage: S('101')
            101
            sage: S.gen(0)
            0
            sage: S.gen(1)
            1
        """
        ## There should really some careful type checking here...
        if isinstance(x, StringMonoidElement) and x.parent() == self:
            return x
        elif isinstance(x, list):
            return StringMonoidElement(self, x, check)
        elif isinstance(x, str):
            return StringMonoidElement(self, x, check)
        else:
            raise TypeError("Argument x (= %s) is of the wrong type." % x)

    def encoding(self, S, padic=False):
        r"""
        The binary encoding of the string ``S``, as a binary string element.

        The default is to keep the standard ASCII byte encoding, e.g.

        ::

            A = 65 -> 01000001
            B = 66 -> 01000010
            .
            .
            .
            Z = 90 -> 01001110

        rather than a 2-adic representation 65 -> 10000010.

        Set ``padic=True`` to reverse the bit string.

        EXAMPLES::

            sage: S = BinaryStrings()
            sage: S.encoding('A')
            01000001
            sage: S.encoding('A',padic=True)
            10000010
            sage: S.encoding(' ',padic=True)
            00000100
        """
        bit_string = []
        for i in range(len(S)):
            n = ord(S[i])
            bits = []
            for i in range(8):
                bits.append(n%2)
                n = n >> 1
            if not padic:
                bits.reverse()
            bit_string.extend(bits)
        return self(bit_string)

    # def ngens(self):
    #     r"""
    #     Return the number of generators of this free binary string monoid.
    #     There are only 2 elements in the binary number system. Hence, this
    #     is the number of generators.

    #     EXAMPLES::

    #         sage: S = BinaryStrings()
    #         sage: S.ngens()
    #         2
    #     """
    #     return 2

class OctalStringMonoid(StringMonoid_class):
    r"""
    The free octal string monoid on generators `\{ 0, 1, \dots, 7 \}`.
    """

    def __init__(self):
        r"""
        Create free octal string monoid on generators `\{ 0, 1, \dots, 7 \}`.

        EXAMPLES::

            sage: S = OctalStrings(); S
            Free octal string monoid
            sage: x = S.gens()
            sage: (x[0]*x[7])**3 * (x[0]*x[1]*x[6]*x[5])**2
            07070701650165
            sage: S([ i for i in range(8) ])
            01234567
        """
        StringMonoid_class.__init__(self, 8, [ str(i) for i in range(8) ])

    def __cmp__(self, other):
        if not isinstance(other, OctalStringMonoid):
            return -1
        return 0

    def __repr__(self):
        return "Free octal string monoid"

    def __call__(self, x, check=True):
        r"""
        Return ``x`` coerced into this free monoid.

        One can create a free octal string monoid element from a
        Python string of 0's to 7's or list of integers.

        EXAMPLES::

            sage: S = OctalStrings()
            sage: S('07070701650165')
            07070701650165
            sage: S.gen(0)
            0
            sage: S.gen(1)
            1
            sage: S([ i for i in range(8) ])
            01234567
        """
        ## There should really some careful type checking here...
        if isinstance(x, StringMonoidElement) and x.parent() == self:
            return x
        elif isinstance(x, list):
            return StringMonoidElement(self, x, check)
        elif isinstance(x, str):
            return StringMonoidElement(self, x, check)
        else:
            raise TypeError("Argument x (= %s) is of the wrong type." % x)

class HexadecimalStringMonoid(StringMonoid_class):
    r"""
    The free hexadecimal string monoid on generators
    `\{ 0, 1, \dots, 9, a, b, c, d, e, f \}`.
    """

    def __init__(self):
        r"""
        Create free hexadecimal string monoid on generators
        `\{ 0, 1, \dots, 9, a, b, c, d, e, f \}`.

        EXAMPLES::

            sage: S = HexadecimalStrings(); S
            Free hexadecimal string monoid
            sage: x = S.gens()
            sage: (x[0]*x[10])**3 * (x[0]*x[1]*x[9]*x[15])**2
            0a0a0a019f019f
            sage: S([ i for i in range(16) ])
            0123456789abcdef
        """
        alph = '0123456789abcdef'
        StringMonoid_class.__init__(self, 16, [ alph[i] for i in range(16) ])

    def __cmp__(self, other):
        if not isinstance(other, HexadecimalStringMonoid):
            return -1
        return 0

    def __repr__(self):
        return "Free hexadecimal string monoid"

    def __call__(self, x, check=True):
        r"""
        Return ``x`` coerced into this free monoid.

        One can create a free hexadecimal string monoid element from a
        Python string of a list of integers in `\{ 0, \dots, 15 \}`.

        EXAMPLES::

            sage: S = HexadecimalStrings()
            sage: S('0a0a0a019f019f')
            0a0a0a019f019f
            sage: S.gen(0)
            0
            sage: S.gen(1)
            1
            sage: S([ i for i in range(16) ])
            0123456789abcdef
        """
        ## There should really some careful type checking here...
        if isinstance(x, StringMonoidElement) and x.parent() == self:
            return x
        elif isinstance(x, list):
            return StringMonoidElement(self, x, check)
        elif isinstance(x, str):
            return StringMonoidElement(self, x, check)
        else:
            raise TypeError("Argument x (= %s) is of the wrong type." % x)

    def encoding(self, S, padic=False):
        r"""
        The encoding of the string ``S`` as a hexadecimal string element.

        The default is to keep the standard right-to-left byte encoding, e.g.

        ::

            A = '\x41' -> 41
            B = '\x42' -> 42
            .
            .
            .
            Z = '\x5a' -> 5a

        rather than a left-to-right representation A = 65 -> 14.
        Although standard (e.g., in the Python constructor '\xhh'),
        this can be confusing when the string reads left-to-right.

        Set ``padic=True`` to reverse the character encoding.

        EXAMPLES::

            sage: S = HexadecimalStrings()
            sage: S.encoding('A')
            41
            sage: S.encoding('A',padic=True)
            14
            sage: S.encoding(' ',padic=False)
            20
            sage: S.encoding(' ',padic=True)
            02
        """
        hex_string = []
        for i in range(len(S)):
            n = ord(S[i])
            n0 = n % 16
            n1 = n // 16
            if not padic:
                hex_chars = [n1, n0]
            else:
                hex_chars = [n0, n1]
            hex_string.extend(hex_chars)
        return self(hex_string)

class Radix64StringMonoid(StringMonoid_class):
    r"""
    The free radix 64 string monoid on 64 generators.
    """

    def __init__(self):
        r"""
        Create free radix 64 string monoid on 64 generators.

        EXAMPLES::

            sage: S = Radix64Strings(); S
            Free radix 64 string monoid
            sage: x = S.gens()
            sage: (x[50]*x[10])**3 * (x[60]*x[1]*x[19]*x[35])**2
            yKyKyK8BTj8BTj
            sage: S([ i for i in range(64) ])
            ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/
        """
        alph = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/'
        StringMonoid_class.__init__(self, 64, [ alph[i] for i in range(64) ])

    def __cmp__(self, other):
        if not isinstance(other, Radix64StringMonoid):
            return -1
        return 0

    def __repr__(self):
        return "Free radix 64 string monoid"

    def __call__(self, x, check=True):
        r"""
        Return ``x`` coerced into this free monoid.

        One can create a free radix 64 string monoid element from a
        Python string or a list of integers in `0, \dots, 63`, as for
        generic ``FreeMonoids``.

        EXAMPLES::

            sage: S = Radix64Strings()
            sage: S.gen(0)
            A
            sage: S.gen(1)
            B
            sage: S.gen(62)
            +
            sage: S.gen(63)
            /
            sage: S([ i for i in range(64) ])
            ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/
        """
        ## There should really some careful type checking here...
        if isinstance(x, StringMonoidElement) and x.parent() == self:
            return x
        elif isinstance(x, list):
            return StringMonoidElement(self, x, check)
        elif isinstance(x, str):
            return StringMonoidElement(self, x, check)
        else:
            raise TypeError("Argument x (= %s) is of the wrong type." % x)

class AlphabeticStringMonoid(StringMonoid_class):
    """
    The free alphabetic string monoid on generators A-Z.

    EXAMPLES::

        sage: S = AlphabeticStrings(); S
        Free alphabetic string monoid on A-Z
        sage: S.gen(0)
        A
        sage: S.gen(25)
        Z
        sage: S([ i for i in range(26) ])
        ABCDEFGHIJKLMNOPQRSTUVWXYZ
    """

    def __init__(self):
        r"""
        Create free alphabetic string monoid on generators A-Z.

        EXAMPLES::

            sage: S = AlphabeticStrings(); S
            Free alphabetic string monoid on A-Z
            sage: S.gen(0)
            A
            sage: S.gen(25)
            Z
            sage: S([ i for i in range(26) ])
            ABCDEFGHIJKLMNOPQRSTUVWXYZ
        """
        from sage.rings.all import RealField
        RR = RealField()
        # The characteristic frequency probability distribution of
        # Robert Edward Lewand.
        self._characteristic_frequency_lewand = {
            "A": RR(0.08167), "B": RR(0.01492),
            "C": RR(0.02782), "D": RR(0.04253),
            "E": RR(0.12702), "F": RR(0.02228),
            "G": RR(0.02015), "H": RR(0.06094),
            "I": RR(0.06966), "J": RR(0.00153),
            "K": RR(0.00772), "L": RR(0.04025),
            "M": RR(0.02406), "N": RR(0.06749),
            "O": RR(0.07507), "P": RR(0.01929),
            "Q": RR(0.00095), "R": RR(0.05987),
            "S": RR(0.06327), "T": RR(0.09056),
            "U": RR(0.02758), "V": RR(0.00978),
            "W": RR(0.02360), "X": RR(0.00150),
            "Y": RR(0.01974), "Z": RR(0.00074)}
        # The characteristic frequency probability distribution of
        # H. Beker and F. Piper.
        self._characteristic_frequency_beker_piper = {
            "A": RR(0.082), "B": RR(0.015),
            "C": RR(0.028), "D": RR(0.043),
            "E": RR(0.127), "F": RR(0.022),
            "G": RR(0.020), "H": RR(0.061),
            "I": RR(0.070), "J": RR(0.002),
            "K": RR(0.008), "L": RR(0.040),
            "M": RR(0.024), "N": RR(0.067),
            "O": RR(0.075), "P": RR(0.019),
            "Q": RR(0.001), "R": RR(0.060),
            "S": RR(0.063), "T": RR(0.091),
            "U": RR(0.028), "V": RR(0.010),
            "W": RR(0.023), "X": RR(0.001),
            "Y": RR(0.020), "Z": RR(0.001)}
        alph = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        StringMonoid_class.__init__(self, 26, [ alph[i] for i in range(26) ])

    def __cmp__(self, other):
        if not isinstance(other, AlphabeticStringMonoid):
            return -1
        return 0

    def __repr__(self):
        return "Free alphabetic string monoid on A-Z"

    def __call__(self, x, check=True):
        r"""
        Return ``x`` coerced into this free monoid.

        One can create a free alphabetic string monoid element from a
        Python string, or a list of integers in `0, \dots,25`.

        EXAMPLES::

            sage: S = AlphabeticStrings()
            sage: S.gen(0)
            A
            sage: S.gen(1)
            B
            sage: S.gen(25)
            Z
            sage: S([ i for i in range(26) ])
            ABCDEFGHIJKLMNOPQRSTUVWXYZ
        """
        ## There should really some careful type checking here...
        if isinstance(x, StringMonoidElement) and x.parent() == self:
            return x
        elif isinstance(x, list):
            return StringMonoidElement(self, x, check)
        elif isinstance(x, str):
            return StringMonoidElement(self, x, check)
        else:
            raise TypeError("Argument x (= %s) is of the wrong type." % x)

    def characteristic_frequency(self, table_name="beker_piper"):
        r"""
        Return a table of the characteristic frequency probability
        distribution of the English alphabet. In written English, various
        letters of the English alphabet occur more frequently than others.
        For example, the letter "E" appears more often than other
        vowels such as "A", "I", "O", and "U". In long works of written
        English such as books, the probability of a letter occurring tends
        to stabilize around a value. We call this value the characteristic
        frequency probability of the letter under consideration. When this
        probability is considered for each letter of the English alphabet,
        the resulting probabilities for all letters of this alphabet is
        referred to as the characteristic frequency probability distribution.
        Various studies report slightly different values for the
        characteristic frequency probability of an English letter. For
        instance, [Lew00]_ reports that "E" has a characteristic
        frequency probability of 0.12702, while [BekPip82]_ reports this
        value as 0.127. The concepts of characteristic frequency probability
        and characteristic frequency probability distribution can also be
        applied to non-empty alphabets other than the English alphabet.

        The output of this method is different from that of the method
        :func:`frequency_distribution()
        <sage.monoids.string_monoid_element.StringMonoidElement.frequency_distribution>`.
        One can think of the characteristic frequency probability of an
        element in an alphabet `A` as the expected probability of that element
        occurring. Let `S` be a string encoded using elements of `A`. The
        frequency probability distribution corresponding to `S` provides us
        with the frequency probability of each element of `A` as observed
        occurring in `S`. Thus one distribution provides expected
        probabilities, while the other provides observed probabilities.

        INPUT:

        - ``table_name`` -- (default ``"beker_piper"``) the table of
          characteristic frequency probability distribution to use. The
          following tables are supported:

          - ``"beker_piper"`` -- the table of characteristic frequency
            probability distribution by Beker and Piper [BekPip82]_. This is
            the default table to use.

          - ``"lewand"`` -- the table of characteristic frequency
            probability distribution by Lewand as described on page 36
            of [Lew00]_.

        OUTPUT:

        - A table of the characteristic frequency probability distribution
          of the English alphabet. This is a dictionary of letter/probability
          pairs.

        EXAMPLES:

        The characteristic frequency probability distribution table of
        Beker and Piper [BekPip82]_::

            sage: A = AlphabeticStrings()
            sage: table = A.characteristic_frequency(table_name="beker_piper")
            sage: sorted(table.items())
            <BLANKLINE>
            [('A', 0.0820000000000000),
            ('B', 0.0150000000000000),
            ('C', 0.0280000000000000),
            ('D', 0.0430000000000000),
            ('E', 0.127000000000000),
            ('F', 0.0220000000000000),
            ('G', 0.0200000000000000),
            ('H', 0.0610000000000000),
            ('I', 0.0700000000000000),
            ('J', 0.00200000000000000),
            ('K', 0.00800000000000000),
            ('L', 0.0400000000000000),
            ('M', 0.0240000000000000),
            ('N', 0.0670000000000000),
            ('O', 0.0750000000000000),
            ('P', 0.0190000000000000),
            ('Q', 0.00100000000000000),
            ('R', 0.0600000000000000),
            ('S', 0.0630000000000000),
            ('T', 0.0910000000000000),
            ('U', 0.0280000000000000),
            ('V', 0.0100000000000000),
            ('W', 0.0230000000000000),
            ('X', 0.00100000000000000),
            ('Y', 0.0200000000000000),
            ('Z', 0.00100000000000000)]

        The characteristic frequency probability distribution table
        of Lewand [Lew00]_::

            sage: table = A.characteristic_frequency(table_name="lewand")
            sage: sorted(table.items())
            <BLANKLINE>
            [('A', 0.0816700000000000),
            ('B', 0.0149200000000000),
            ('C', 0.0278200000000000),
            ('D', 0.0425300000000000),
            ('E', 0.127020000000000),
            ('F', 0.0222800000000000),
            ('G', 0.0201500000000000),
            ('H', 0.0609400000000000),
            ('I', 0.0696600000000000),
            ('J', 0.00153000000000000),
            ('K', 0.00772000000000000),
            ('L', 0.0402500000000000),
            ('M', 0.0240600000000000),
            ('N', 0.0674900000000000),
            ('O', 0.0750700000000000),
            ('P', 0.0192900000000000),
            ('Q', 0.000950000000000000),
            ('R', 0.0598700000000000),
            ('S', 0.0632700000000000),
            ('T', 0.0905600000000000),
            ('U', 0.0275800000000000),
            ('V', 0.00978000000000000),
            ('W', 0.0236000000000000),
            ('X', 0.00150000000000000),
            ('Y', 0.0197400000000000),
            ('Z', 0.000740000000000000)]

        Illustrating the difference between :func:`characteristic_frequency`
        and :func:`frequency_distribution() <sage.monoids.string_monoid_element.StringMonoidElement.frequency_distribution>`::

            sage: A = AlphabeticStrings()
            sage: M = A.encoding("abcd")
            sage: FD = M.frequency_distribution().function()
            sage: sorted(FD.items())
            <BLANKLINE>
            [(A, 0.250000000000000),
            (B, 0.250000000000000),
            (C, 0.250000000000000),
            (D, 0.250000000000000)]
            sage: CF = A.characteristic_frequency()
            sage: sorted(CF.items())
            <BLANKLINE>
            [('A', 0.0820000000000000),
            ('B', 0.0150000000000000),
            ('C', 0.0280000000000000),
            ('D', 0.0430000000000000),
            ('E', 0.127000000000000),
            ('F', 0.0220000000000000),
            ('G', 0.0200000000000000),
            ('H', 0.0610000000000000),
            ('I', 0.0700000000000000),
            ('J', 0.00200000000000000),
            ('K', 0.00800000000000000),
            ('L', 0.0400000000000000),
            ('M', 0.0240000000000000),
            ('N', 0.0670000000000000),
            ('O', 0.0750000000000000),
            ('P', 0.0190000000000000),
            ('Q', 0.00100000000000000),
            ('R', 0.0600000000000000),
            ('S', 0.0630000000000000),
            ('T', 0.0910000000000000),
            ('U', 0.0280000000000000),
            ('V', 0.0100000000000000),
            ('W', 0.0230000000000000),
            ('X', 0.00100000000000000),
            ('Y', 0.0200000000000000),
            ('Z', 0.00100000000000000)]

        TESTS:

        The table name must be either "beker_piper" or "lewand"::

            sage: table = A.characteristic_frequency(table_name="")
            Traceback (most recent call last):
            ...
            ValueError: Table name must be either 'beker_piper' or 'lewand'.
            sage: table = A.characteristic_frequency(table_name="none")
            Traceback (most recent call last):
            ...
            ValueError: Table name must be either 'beker_piper' or 'lewand'.

        REFERENCES:

        .. [BekPip82] H. Beker and F. Piper. *Cipher Systems: The
          Protection of Communications*. John Wiley and Sons, 1982.

        .. [Lew00] Robert Edward Lewand. *Cryptological Mathematics*.
          The Mathematical Association of America, 2000.
        """
        supported_tables = ["beker_piper", "lewand"]
        if table_name not in supported_tables:
            raise ValueError(
                "Table name must be either 'beker_piper' or 'lewand'.")
        from copy import copy
        if table_name == "beker_piper":
            return copy(self._characteristic_frequency_beker_piper)
        if table_name == "lewand":
            return copy(self._characteristic_frequency_lewand)

    def encoding(self, S):
        r"""
        The encoding of the string ``S`` in the alphabetic string monoid,
        obtained by the monoid homomorphism

        ::

            A -> A, ..., Z -> Z, a -> A, ..., z -> Z

        and stripping away all other characters. It should be noted that
        this is a non-injective monoid homomorphism.

        EXAMPLES::

            sage: S = AlphabeticStrings()
            sage: s = S.encoding("The cat in the hat."); s
            THECATINTHEHAT
            sage: s.decoding()
            'THECATINTHEHAT'
        """
        return self(strip_encoding(S))
