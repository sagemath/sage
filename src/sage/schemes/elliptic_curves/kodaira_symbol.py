r"""
Kodaira symbols

Kodaira symbols encode the type of reduction of an elliptic curve at a
(finite) place.

The standard notation for Kodaira Symbols is as a string which is one
of `\rm{I}_m`, `\rm{II}`, `\rm{III}`, `\rm{IV}`, `\rm{I}^*_m`,
`\rm{II}^*`, `\rm{III}^*`, `\rm{IV}^*`, where `m` denotes a
non-negative integer.  These have been encoded by single integers by
different people.  For convenience we give here the conversion table
between strings, the eclib coding and the PARI encoding.

+----------------------+----------------+--------------------+
| Kodaira Symbol       |  Eclib coding  |  PARI Coding       |
+======================+================+====================+
| `\rm{I}_0`           |      `0`       |   `1`              |
+----------------------+----------------+--------------------+
| `\rm{I}^*_0`         |      `1`       |   `-1`             |
+----------------------+----------------+--------------------+
| `\rm{I}_m`  `(m>0)`  |      `10m`     |   `m+4`            |
+----------------------+----------------+--------------------+
| `\rm{I}^*_m` `(m>0)` |      `10m+1`   |   `-(m+4)`         |
+----------------------+----------------+--------------------+
| `\rm{II}`            | `2`            |   `2`              |
+----------------------+----------------+--------------------+
| `\rm{III}`           | `3`            |   `3`              |
+----------------------+----------------+--------------------+
| `\rm{IV}`            | `4`            |   `4`              |
+----------------------+----------------+--------------------+
| `\rm{II}^*`          | `7`            |   `-2`             |
+----------------------+----------------+--------------------+
| `\rm{III}^*`         | `6`            |   `-3`             |
+----------------------+----------------+--------------------+
| `\rm{IV}^*`          | `5`            |   `-4`             |
+----------------------+----------------+--------------------+


AUTHORS:

- David Roe       <roed@math.harvard.edu>

- John Cremona

- Katie Ahrens (unicode art method)
"""

# ****************************************************************************
#       Copyright (C) 2007 David Roe       <roed@math.harvard.edu>
#                          William Stein   <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.structure.sage_object import SageObject
from sage.structure.richcmp import richcmp_method, richcmp
from sage.rings.integer import Integer
from sage.env import SAGE_ROOT
import weakref


@richcmp_method
class KodairaSymbol_class(SageObject):
    r"""
    Class to hold a Kodaira symbol of an elliptic curve over a
    `p`-adic local field.

    Users should use the ``KodairaSymbol()`` function to construct
    Kodaira Symbols rather than use the class constructor directly.
    """

    def __init__(self, symbol):
        r"""
        Constructor for Kodaira Symbol class.

        INPUT:

        - ``symbol`` (string or integer) -- The string should be a
          standard string representation (e.g. III*) of a Kodaira
          symbol, which will be parsed.  Alternatively, use the PARI
          encoding of Kodaira symbols as integers.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.kodaira_symbol import KodairaSymbol_class
            sage: KodairaSymbol_class(14)
            I10
            sage: KodairaSymbol_class('III*')
            III*
            sage: latex(KodairaSymbol_class('In'))
            I_n
            sage: KodairaSymbol_class('In')
            In

        Check that :issue:`31147` is fixed::

            sage: latex(KodairaSymbol_class(-14))
            I_{10}^{*}
        """
        if not isinstance(symbol, str):
            n = Integer(symbol)
            self._n = None
            if n == 0:
                raise ValueError("Kodaira Symbol code number must be nonzero.")
            if n == 1:
                self._n = 0
                self._roman = 1
                self._str = 'I0'
                self._latex = 'I_0'
            elif n == 2:
                self._roman = 2
                self._str = 'II'
                self._latex = 'II'
            elif n == 3:
                self._roman = 3
                self._str = 'III'
                self._latex = 'III'
            elif n == 4:
                self._roman = 4
                self._str = 'IV'
                self._latex = 'IV'
            elif n > 4:
                nu = n - 4
                self._n = nu
                self._roman = 1
                self._str = 'I%s' % nu
                self._latex = 'I_{%s}' % nu
            elif n == -1:
                self._roman = 1
                self._n = 0
                self._str = 'I0*'
                self._latex = 'I_0^{*}'
            elif n == -2:
                self._roman = 2
                self._str = 'II*'
                self._latex = 'II^{*}'
            elif n == -3:
                self._roman = 3
                self._str = 'III*'
                self._latex = 'III^{*}'
            elif n == -4:
                self._roman = 4
                self._str = 'IV*'
                self._latex = 'IV^{*}'
            elif n < -4:
                nu = -n - 4
                self._roman = 1
                self._n = nu
                self._str = 'I%s*' % nu
                self._latex = 'I_{%s}^{*}' % nu
            self._starred = (n < 0)
            self._pari = n
            return
        elif not symbol:
            raise TypeError("symbol must be a nonempty string")
        if symbol[0] == "I":
            symbol = symbol[1:]
        starred = False
        if symbol[-1] == "*":
            starred = True
            symbol = symbol[:-1]
        self._starred = starred
        if symbol in ["I", "II", "V"]:    # NB we have already stripped off the leading 'I'
            self._roman = ["I", "II", "V"].index(symbol) + 2   # =2, 3 or 4
            self._n = None
            if starred:
                sign = -1
                self._str = "I" + symbol + "*"
                self._latex = "I" + symbol + "^*"
            else:
                sign = 1
                self._str = "I" + symbol
                self._latex = "" + self._str + ""
            if symbol == "I":
                self._pari = 2 * sign
            elif symbol == "II":
                self._pari = 3 * sign
            elif symbol == "V":
                self._pari = 4 * sign
        elif symbol == "n":
            self._roman = 1
            self._pari = None
            self._n = "generic"
            if starred:
                self._str = "In*"
                self._latex = "I_n^*"
            else:
                self._str = "In"
                self._latex = "I_n"
        elif symbol.isdigit():
            self._roman = 1
            self._n = Integer(symbol)
            if starred:
                if self._n == 0:
                    self._pari = -1
                else:
                    self._pari = -self._n - 4
                self._str = "I" + symbol + "*"
                self._latex = "I_{%s}^*" % (symbol)
            else:
                if self._n == 0:
                    self._pari = 1
                else:
                    self._pari = self._n + 4
                self._str = "I" + symbol
                self._latex = "I_{%s}" % (symbol)
        else:
            raise ValueError("input is not a Kodaira symbol")

    def __repr__(self):
        r"""
        Return the string representation of this Kodaira Symbol.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.kodaira_symbol import KodairaSymbol_class
            sage: KS = KodairaSymbol_class(15)
            sage: str(KS) # indirect doctest
            'I11'
        """
        return self._str

    def _latex_(self):
        r"""
        Return the string representation of this Kodaira Symbol.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.kodaira_symbol import KodairaSymbol_class
            sage: KS = KodairaSymbol_class(15)
            sage: latex(KS)
            I_{11}
        """
        return self._latex

    def __richcmp__(self, other, op):
        r"""
        Standard comparison function for Kodaira Symbols.

        EXAMPLES::

            sage: from sage.schemes.elliptic_curves.kodaira_symbol import KodairaSymbol_class
            sage: KS1 = KodairaSymbol_class(15); KS1
            I11
            sage: KS2 = KodairaSymbol_class(-34); KS2
            I30*
            sage: KS1 < KS2
            True
            sage: KS2 < KS1
            False

        ::

            sage: Klist = [KodairaSymbol_class(i) for i in [-10..10] if i!=0]
            sage: Klist.sort()
            sage: Klist
            [I0,
            I0*,
            I1,
            I1*,
            I2,
            I2*,
            I3,
            I3*,
            I4,
            I4*,
            I5,
            I5*,
            I6,
            I6*,
            II,
            II*,
            III,
            III*,
            IV,
            IV*]
        """
        if isinstance(other, KodairaSymbol_class):
            if (self._n == "generic" and other._n is not None) or (other._n == "generic" and self._n is not None):
                return richcmp(self._starred, other._starred, op)
            return richcmp(self._str, other._str, op)
        else:
            return NotImplemented

    def _pari_code(self):
        """
        Return the PARI encoding of this Kodaira Symbol.

        EXAMPLES::

            sage: KodairaSymbol('I0')._pari_code()
            1
            sage: KodairaSymbol('I10')._pari_code()
            14
            sage: KodairaSymbol('I10*')._pari_code()
            -14
            sage: [KodairaSymbol(s)._pari_code() for s in ['II','III','IV']]
            [2, 3, 4]
            sage: [KodairaSymbol(s)._pari_code() for s in ['II*','III*','IV*']]
            [-2, -3, -4]
        """
        return self._pari

    def _unicode_art_(self):
        r"""
        Return a unicode art representation of the Kodaira symbol.

        EXAMPLES::

            sage: unicode_art(KodairaSymbol('In'))
            _/___\_
            /  1  \ 1
            :     :
            :     :
            _\____/_
            1 \  / 1
               1
            sage: unicode_art(KodairaSymbol('III*'))
            \   \ / 3
            1 \ 2/ \
              / \   \
             /   \  |\
                    |
             _______|
                2   |
                   4|
                    |/
               \/   /
               /\  /|
             1/ 2\/3|
                /   /\
        """
        from sage.typeset.unicode_art import unicode_art
        string = str(self)
        if string.startswith('I') and list(string)[1] == '0' and list(string)[-1] == '*':
            #I_0*
            return unicode_art("1\n\\ |\n \\|\n\\ |\\\n1\\|\n  \\\n  |\n  |\n  |/\n1/| /\n/ |/\n /|\n/ |\n1")
        elif string.startswith('I') and list(string)[1].isdigit() and list(string)[-1] == '*':
            #I_n*
            return unicode_art("   \\ 1\n    \\\n  \\  \\\n 1 \\  \\\n    \\  \\/\n     \\/ \\\n    / \\\n \\ / 2\n  /\\\n /  \\\n    2 \\\n       :\n \\  2 /\n  \\  /\n   \\/\n   /\\\n  /  \\2\n      \\/\n     / \\\n   1/   \\/\n   /   / \\\n      /   \\\n     /\n    /")
        elif string.startswith('I') and list(string)[1].isdigit():
            #I_n
            return unicode_art(" _/___\\_\n /  1  \\ 1\n :     :\n :     :\n_\\____/_\n1 \\  /  1\n    1")
        elif string == "In":
            #also I_n
            return unicode_art(" _/___\\_\n /  1  \\ 1\n :     :\n :     :\n_\\____/_\n1 \\  /  1\n    1")
            #string = "I_n"
        elif string == "II":
            return unicode_art("\\     /\n \\   /\n  \\ /\n   |\n1  |")
        elif string == "III":
            return unicode_art("  \\     /\n   \\   /\n    \\ /\n     |\n     |\n    /\\\n   /  \\\n1 /    \\ 1")
        elif string == "IV":
            return unicode_art(" \\    /\n  \\  /\n1__\\/__\n   /\\\n1 /  \\ 1\n /    \\")
        elif string == "IV*":
            return unicode_art("     \\ /\n    /  \\\n   /     \\|\n  /       |\\\n          |\n 1|       |\n _|_______|\n  |   2   |\n          |\n          |\n        3 |/\n  \\       /\n1  \\    / |\n     \\/ 2 |\n     /\\")
        elif string == "III*":
            return unicode_art("  \\   \\ / 3\n1  \\ 2/ \\\n   / \\   \\\n  /   \\  |\\\n         |\n  _______|\n     2   |\n        4|\n         |/\n    \\/   /\n    /\\  /|\n  1/ 2\\/3|\n  /   /\\")
        elif string == "II*":
            return unicode_art(" \\   2/\\   \\/ 5\n1  \\/  3\\4/ \\\n  / \\  / \\   \\\n /   /\\   \\  |\\\n             |\n             |\n             |\n             |\n       \\     |/\n        \\   /|\n        2\\ /4\n         / \\")
        else:
            return


_ks_cache = {}


def KodairaSymbol(symbol):
    r"""
    Return the specified Kodaira symbol.

    INPUT:

    - ``symbol`` (string or integer) -- Either a string of the form
      "I0", "I1", ..., "In", "II", "III", "IV", "I0*", "I1*", ..., "In*", "II*", "III*", or "IV*",
      or an integer encoding a Kodaira symbol using PARI's conventions.

    OUTPUT:

    (KodairaSymbol)  The corresponding Kodaira symbol.

    EXAMPLES::

        sage: KS = KodairaSymbol
        sage: [KS(n) for n in range(1,10)]
        [I0, II, III, IV, I1, I2, I3, I4, I5]
        sage: [KS(-n) for n in range(1,10)]
        [I0*, II*, III*, IV*, I1*, I2*, I3*, I4*, I5*]
        sage: all(KS(str(KS(n))) == KS(n) for n in range(-10,10) if n != 0)
        True
        sage: unicode_art(KS(2))
         \     /
          \   /
           \ /
            |
         1  |

    """

    if symbol in _ks_cache:
        ks = _ks_cache[symbol]()
        if ks is not None:
            return ks

    ks = KodairaSymbol_class(symbol)
    _ks_cache[symbol] = weakref.ref(ks)
    return ks
