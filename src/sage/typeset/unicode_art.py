# -*- coding: utf-8 -*-
r"""
Unicode Art

This module implements ascii art using unicode characters. It is a
strict superset of :mod:`~sage.typeset.ascii_art`.
"""

# ******************************************************************************
#       Copyright (C) 2013 Jean-Baptiste Priez <jbp@kerios.fr>,
#                     2015 Volker Braun <vbraun.name@gmail.com>
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
# ******************************************************************************

from sage.typeset.character_art import CharacterArt
from sage.typeset.character_art_factory import CharacterArtFactory
import sage.typeset.symbols as symbol


class UnicodeArt(CharacterArt):
    r"""
    An Ascii art object is an object with some specific representation for
    *printing*.

    INPUT:

    - ``lines`` -- the list of lines of the representation of the ascii art
      object

    - ``breakpoints`` -- the list of points where the representation can be
      split

    - ``baseline`` -- the reference line (from the bottom)

    EXAMPLES::

        sage: i = var('i')
        sage: unicode_art(sum(pi^i/factorial(i)*x^i, i, 0, oo))
         π⋅x
        ℯ
    """
    _string_type = str


_unicode_art_factory = CharacterArtFactory(
    UnicodeArt, str, '_unicode_art_',
    (symbol.unicode_left_parenthesis, symbol.unicode_right_parenthesis),
    (symbol.unicode_left_square_bracket, symbol.unicode_right_square_bracket),
    (symbol.unicode_left_curly_brace, symbol.unicode_right_curly_brace),
)


empty_unicode_art = _unicode_art_factory.build_empty()


def unicode_art(*obj, **kwds):
    r"""
    Return an unicode art representation

    INPUT:


    - ``*obj`` -- any number of positional arguments, of arbitrary
      type. The objects whose ascii art representation we want.

    - ``sep`` -- optional ``'sep=...'`` keyword argument (or ``'separator'``).
      Anything that can be converted to unicode art (default: empty unicode
      art). The separator in-between a list of objects. Only used if
      more than one object given.

    - ``baseline`` -- (default: 0) the baseline for the object

    - ``sep_baseline`` -- (default: 0) the baseline for the separator

    OUTPUT:

    :class:`UnicodeArt` instance.

    EXAMPLES::

        sage: result = unicode_art(integral(exp(sqrt(x))/(x+pi), x))
        ...
        sage: result
            ⌠
            ⎮   √x
            ⎮  ℯ
            ⎮ ───── dx
            ⎮ x + π
            ⌡
        sage: ident = lambda n: identity_matrix(ZZ, n)
        sage: unicode_art(ident(1), ident(2), ident(3), sep=' : ')
                      ⎛1 0 0⎞
              ⎛1 0⎞   ⎜0 1 0⎟
        (1) : ⎝0 1⎠ : ⎝0 0 1⎠

    If specified, the ``sep_baseline`` overrides the baseline of
    an unicode art separator::

        sage: sep_line = unicode_art('\n'.join(' ⎟ ' for _ in range(5)), baseline=5)
        sage: unicode_art(*AlternatingSignMatrices(3),
        ....:             separator=sep_line, sep_baseline=1)
                ⎟         ⎟         ⎟            ⎟         ⎟         ⎟
        ⎛1 0 0⎞ ⎟ ⎛0 1 0⎞ ⎟ ⎛1 0 0⎞ ⎟ ⎛ 0  1  0⎞ ⎟ ⎛0 0 1⎞ ⎟ ⎛0 1 0⎞ ⎟ ⎛0 0 1⎞
        ⎜0 1 0⎟ ⎟ ⎜1 0 0⎟ ⎟ ⎜0 0 1⎟ ⎟ ⎜ 1 -1  1⎟ ⎟ ⎜1 0 0⎟ ⎟ ⎜0 0 1⎟ ⎟ ⎜0 1 0⎟
        ⎝0 0 1⎠ ⎟ ⎝0 0 1⎠ ⎟ ⎝0 1 0⎠ ⎟ ⎝ 0  1  0⎠ ⎟ ⎝0 1 0⎠ ⎟ ⎝1 0 0⎠ ⎟ ⎝1 0 0⎠
                ⎟         ⎟         ⎟            ⎟         ⎟         ⎟

    TESTS::

        sage: n = var('n')
        sage: unicode_art(sum(binomial(2 * n, n + 1) * x^n, n, 0, oo))
         ⎛        _________    ⎞
        -⎝2⋅x + ╲╱ 1 - 4⋅x  - 1⎠
        ─────────────────────────
                   _________
             2⋅x⋅╲╱ 1 - 4⋅x
        sage: unicode_art(list(DyckWords(3)))
        ⎡                                   ╱╲   ⎤
        ⎢            ╱╲    ╱╲      ╱╲╱╲    ╱  ╲  ⎥
        ⎣ ╱╲╱╲╱╲, ╱╲╱  ╲, ╱  ╲╱╲, ╱    ╲, ╱    ╲ ⎦
        sage: unicode_art(1)
        1
    """
    separator, baseline, sep_baseline = _unicode_art_factory.parse_keywords(kwds)
    if kwds:
        raise ValueError('unknown keyword arguments: {0}'.format(list(kwds)))
    if len(obj) == 1:
        return _unicode_art_factory.build(obj[0], baseline=baseline)
    if not isinstance(separator, UnicodeArt):
        separator = _unicode_art_factory.build(separator, baseline=sep_baseline)
    elif sep_baseline is not None:
        from copy import copy
        separator = copy(separator)
        separator._baseline = sep_baseline
    return _unicode_art_factory.concatenate(obj, separator, empty_unicode_art,
                                            baseline=baseline)


_subscript_dict = {'0': '₀', '1': '₁', '2': '₂', '3': '₃', '4': '₄',
                   '5': '₅', '6': '₆', '7': '₇', '8': '₈', '9': '₉',
                   '-': '₋', '+': '₊'}


_superscript_dict = {'0': '⁰', '1': '¹', '2': '²', '3': '³', '4': '⁴',
                     '5': '⁵', '6': '⁶', '7': '⁷', '8': '⁸', '9': '⁹',
                     '-': '⁻', '+': '⁺', '/': 'ᐟ'}


def unicode_superscript(x):
    r"""
    Return the rational number ``x`` as a superscript.

    EXAMPLES::

        sage: from sage.typeset.unicode_art import unicode_superscript
        sage: unicode_superscript(15123902)
        '¹⁵¹²³⁹⁰²'
        sage: unicode_superscript(-712/5)
        '⁻⁷¹²ᐟ⁵'
    """
    return ''.join(_superscript_dict[i] for i in str(x))


def unicode_subscript(x):
    r"""
    Return the integer ``x`` as a superscript.

    EXAMPLES::

        sage: from sage.typeset.unicode_art import unicode_subscript
        sage: unicode_subscript(15123902)
        '₁₅₁₂₃₉₀₂'
        sage: unicode_subscript(-712)
        '₋₇₁₂'
    """
    return ''.join(_subscript_dict[i] for i in str(x))
