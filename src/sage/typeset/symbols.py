# sage_setup: distribution = sagemath-categories
"""
Symbols for Character Art

This module defines single- and multi-line symbols.

EXAMPLES::

    sage: from sage.typeset.symbols import *

    sage: symbols = ascii_art('')
    sage: for i in range(1, 5):
    ....:     symbols += ascii_left_parenthesis.character_art(i)
    ....:     symbols += ascii_art(' ')
    ....:     symbols += ascii_right_parenthesis.character_art(i)
    ....:     symbols += ascii_art(' ')
    sage: for i in range(1, 5):
    ....:     symbols += ascii_left_square_bracket.character_art(i)
    ....:     symbols += ascii_art(' ')
    ....:     symbols += ascii_right_square_bracket.character_art(i)
    ....:     symbols += ascii_art(' ')
    sage: for i in range(1, 5):
    ....:     symbols += ascii_left_curly_brace.character_art(i)
    ....:     symbols += ascii_art(' ')
    ....:     symbols += ascii_right_curly_brace.character_art(i)
    ....:     symbols += ascii_art(' ')
    sage: symbols
                ( )             [ ]             { }
            ( ) ( )         [ ] [ ]         { } { }
        ( ) ( ) ( )     [ ] [ ] [ ]     { } { } { }
    ( ) ( ) ( ) ( ) [ ] [ ] [ ] [ ] { } { } { } { }

    sage: symbols = unicode_art('')
    sage: for i in range(1, 5):
    ....:     symbols += unicode_left_parenthesis.character_art(i)
    ....:     symbols += unicode_art(' ')
    ....:     symbols += unicode_right_parenthesis.character_art(i)
    ....:     symbols += unicode_art(' ')
    sage: for i in range(1, 5):
    ....:     symbols += unicode_left_square_bracket.character_art(i)
    ....:     symbols += unicode_art(' ')
    ....:     symbols += unicode_right_square_bracket.character_art(i)
    ....:     symbols += unicode_art(' ')
    sage: for i in range(1, 5):
    ....:     symbols += unicode_left_curly_brace.character_art(i)
    ....:     symbols += unicode_art(' ')
    ....:     symbols += unicode_right_curly_brace.character_art(i)
    ....:     symbols += unicode_art(' ')
    sage: symbols
                ⎛ ⎞             ⎡ ⎤             ⎧ ⎫
            ⎛ ⎞ ⎜ ⎟         ⎡ ⎤ ⎢ ⎥         ⎧ ⎫ ⎭ ⎩
        ⎛ ⎞ ⎜ ⎟ ⎜ ⎟     ⎡ ⎤ ⎢ ⎥ ⎢ ⎥     ⎰ ⎱ ⎨ ⎬ ⎫ ⎧
    ( ) ⎝ ⎠ ⎝ ⎠ ⎝ ⎠ [ ] ⎣ ⎦ ⎣ ⎦ ⎣ ⎦ { } ⎱ ⎰ ⎩ ⎭ ⎩ ⎭
"""
import unicodedata
from sage.structure.sage_object import SageObject


class CompoundSymbol(SageObject):

    def __init__(self, character, top, extension, bottom,
                 middle=None,
                 middle_top=None, middle_bottom=None,
                 top_2=None, bottom_2=None):
        """
        A multi-character (ascii/unicode art) symbol.

        INPUT:

        Instead of a string, each of these can be unicode in Python 2:

        - ``character`` -- string; the single-line version of the symbol

        - ``top`` -- string; the top line of a multi-line symbol

        - ``extension`` -- string; the extension line of a multi-line symbol (will
          be repeated)

        - ``bottom`` -- string; the bottom line of a multi-line symbol

        - ``middle`` -- (optional) string; the middle part, for example
          in curly braces. Will be used only once for the symbol, and
          only if its height is odd.

        - ``middle_top`` -- (optional) string; the upper half of the
          2-line middle part if the height of the symbol is even.
          Will be used only once for the symbol.

        - ``middle_bottom`` -- (optional) string; the lower half of the
          2-line middle part if the height of the symbol is even.
          Will be used only once for the symbol.

        - ``top_2`` -- (optional) string; the upper half of a 2-line symbol

        - ``bottom_2`` -- (optional) string; the lower half of a 2-line symbol

        EXAMPLES::

            sage: from sage.typeset.symbols import CompoundSymbol
            sage: i = CompoundSymbol('I', '+', '|', '+', '|')
            sage: i.print_to_stdout(1)
            I
            sage: i.print_to_stdout(3)
            +
            |
            +
        """
        self.character = character
        self.top = top
        self.extension = extension
        self.bottom = bottom
        self.middle = middle or extension
        self.middle_top = middle_top or extension
        self.middle_bottom = middle_bottom or extension
        self.top_2 = top_2 or top
        self.bottom_2 = bottom_2 or bottom

    def _repr_(self):
        """
        Return string representation.

        EXAMPLES::

            sage: from sage.typeset.symbols import unicode_left_parenthesis
            sage: unicode_left_parenthesis
            multi_line version of "("
        """
        return 'multi_line version of "{0}"'.format(self.character)

    def __call__(self, num_lines):
        r"""
        Return the lines for a multi-line symbol.

        INPUT:

        - ``num_lines`` -- integer; the total number of lines

        OUTPUT: list of strings / unicode strings

        EXAMPLES::

            sage: from sage.typeset.symbols import unicode_left_parenthesis
            sage: unicode_left_parenthesis(4)
            ['\u239b', '\u239c', '\u239c', '\u239d']
        """
        if num_lines <= 0:
            raise ValueError('number of lines must be positive')
        elif num_lines == 1:
            return [self.character]
        elif num_lines == 2:
            return [self.top_2, self.bottom_2]
        elif num_lines == 3:
            return [self.top, self.middle, self.bottom]
        elif num_lines % 2 == 0:
            ext = [self.extension] * ((num_lines - 4) // 2)
            return [self.top] + ext + [self.middle_top, self.middle_bottom] + ext + [self.bottom]
        else:  # num_lines %2 == 1
            ext = [self.extension] * ((num_lines - 3) // 2)
            return [self.top] + ext + [self.middle] + ext + [self.bottom]

    def print_to_stdout(self, num_lines):
        """
        Print the multi-line symbol.

        This method is for testing purposes.

        INPUT:

        - ``num_lines`` -- integer; the total number of lines

        EXAMPLES::

            sage: from sage.typeset.symbols import *
            sage: unicode_integral.print_to_stdout(1)
            ∫
            sage: unicode_integral.print_to_stdout(2)
            ⌠
            ⌡
            sage: unicode_integral.print_to_stdout(3)
            ⌠
            ⎮
            ⌡
            sage: unicode_integral.print_to_stdout(4)
            ⌠
            ⎮
            ⎮
            ⌡
        """
        print('\n'.join(self(num_lines)))


class CompoundAsciiSymbol(CompoundSymbol):

    def character_art(self, num_lines):
        """
        Return the ASCII art of the symbol.

        EXAMPLES::

            sage: from sage.typeset.symbols import *
            sage: ascii_left_curly_brace.character_art(3)
            {
            {
            {
        """
        from sage.typeset.ascii_art import AsciiArt
        return AsciiArt(self(num_lines))


class CompoundUnicodeSymbol(CompoundSymbol):

    def character_art(self, num_lines):
        """
        Return the unicode art of the symbol.

        EXAMPLES::

            sage: from sage.typeset.symbols import *
            sage: unicode_left_curly_brace.character_art(3)
            ⎧
            ⎨
            ⎩
        """
        from sage.typeset.unicode_art import UnicodeArt
        return UnicodeArt(self(num_lines))


ascii_integral = CompoundAsciiSymbol(
    'int',
    r' /\\',
    r' | ',
    r'\\/ ',
)

unicode_integral = CompoundUnicodeSymbol(
    unicodedata.lookup('INTEGRAL'),
    unicodedata.lookup('TOP HALF INTEGRAL'),
    unicodedata.lookup('INTEGRAL EXTENSION'),
    unicodedata.lookup('BOTTOM HALF INTEGRAL'),
)

ascii_left_parenthesis = CompoundAsciiSymbol(
    '(',
    '(',
    '(',
    '(',
)

ascii_right_parenthesis = CompoundAsciiSymbol(
    ')',
    ')',
    ')',
    ')',
)

unicode_left_parenthesis = CompoundUnicodeSymbol(
    unicodedata.lookup('LEFT PARENTHESIS'),
    unicodedata.lookup('LEFT PARENTHESIS UPPER HOOK'),
    unicodedata.lookup('LEFT PARENTHESIS EXTENSION'),
    unicodedata.lookup('LEFT PARENTHESIS LOWER HOOK'),
)

unicode_right_parenthesis = CompoundUnicodeSymbol(
    unicodedata.lookup('RIGHT PARENTHESIS'),
    unicodedata.lookup('RIGHT PARENTHESIS UPPER HOOK'),
    unicodedata.lookup('RIGHT PARENTHESIS EXTENSION'),
    unicodedata.lookup('RIGHT PARENTHESIS LOWER HOOK'),
)

ascii_left_square_bracket = CompoundAsciiSymbol(
    '[',
    '[',
    '[',
    '[',
)

ascii_right_square_bracket = CompoundAsciiSymbol(
    ']',
    ']',
    ']',
    ']',
)

unicode_left_square_bracket = CompoundUnicodeSymbol(
    unicodedata.lookup('LEFT SQUARE BRACKET'),
    unicodedata.lookup('LEFT SQUARE BRACKET UPPER CORNER'),
    unicodedata.lookup('LEFT SQUARE BRACKET EXTENSION'),
    unicodedata.lookup('LEFT SQUARE BRACKET LOWER CORNER'),
)

unicode_right_square_bracket = CompoundUnicodeSymbol(
    unicodedata.lookup('RIGHT SQUARE BRACKET'),
    unicodedata.lookup('RIGHT SQUARE BRACKET UPPER CORNER'),
    unicodedata.lookup('RIGHT SQUARE BRACKET EXTENSION'),
    unicodedata.lookup('RIGHT SQUARE BRACKET LOWER CORNER'),
)

ascii_left_curly_brace = CompoundAsciiSymbol(
    '{',
    '{',
    '{',
    '{',
)

ascii_right_curly_brace = CompoundAsciiSymbol(
    '}',
    '}',
    '}',
    '}',
)

unicode_left_curly_brace = CompoundUnicodeSymbol(
    unicodedata.lookup('LEFT CURLY BRACKET'),
    unicodedata.lookup('LEFT CURLY BRACKET UPPER HOOK'),
    unicodedata.lookup('CURLY BRACKET EXTENSION'),
    unicodedata.lookup('LEFT CURLY BRACKET LOWER HOOK'),
    unicodedata.lookup('LEFT CURLY BRACKET MIDDLE PIECE'),
    unicodedata.lookup('RIGHT CURLY BRACKET LOWER HOOK'),
    unicodedata.lookup('RIGHT CURLY BRACKET UPPER HOOK'),
    unicodedata.lookup('UPPER LEFT OR LOWER RIGHT CURLY BRACKET SECTION'),
    unicodedata.lookup('UPPER RIGHT OR LOWER LEFT CURLY BRACKET SECTION'),
)

unicode_right_curly_brace = CompoundUnicodeSymbol(
    unicodedata.lookup('RIGHT CURLY BRACKET'),
    unicodedata.lookup('RIGHT CURLY BRACKET UPPER HOOK'),
    unicodedata.lookup('CURLY BRACKET EXTENSION'),
    unicodedata.lookup('RIGHT CURLY BRACKET LOWER HOOK'),
    unicodedata.lookup('RIGHT CURLY BRACKET MIDDLE PIECE'),
    unicodedata.lookup('LEFT CURLY BRACKET LOWER HOOK'),
    unicodedata.lookup('LEFT CURLY BRACKET UPPER HOOK'),
    unicodedata.lookup('UPPER RIGHT OR LOWER LEFT CURLY BRACKET SECTION'),
    unicodedata.lookup('UPPER LEFT OR LOWER RIGHT CURLY BRACKET SECTION'),
)
