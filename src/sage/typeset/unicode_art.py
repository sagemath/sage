# -*- coding: utf-8 -*-
r"""
Unicode Art

This module implements ascii art using unicode characters. It is a
strict superset of :mod:`~sage.typeset.ascii_art`.
"""

#*******************************************************************************
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
#                  http://www.gnu.org/licenses/
#*******************************************************************************

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
    _string_type = unicode


_unicode_art_factory = CharacterArtFactory(
    UnicodeArt, unicode, '_unicode_art_',
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

    - ``sep`` -- optional ``'sep=...'`` keyword argument. Anything
      that can be converted to unicode art (default: empty unicode
      art). The separator in-between a list of objects. Only used if
      more than one object given.

    OUTPUT:

    :class:`UnicodeArt` instance.

    EXAMPLES::

        sage: unicode_art(integral(exp(sqrt(x))/(x+pi), x))
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

    TESTS::

        sage: n = var('n')
        sage: unicode_art(sum(binomial(2 * n, n + 1) * x^n, n, 0, oo))
             ⎛        __________    ⎞
            -⎝2⋅x + ╲╱ -4⋅x + 1  - 1⎠
            ──────────────────────────
                       __________
                 2⋅x⋅╲╱ -4⋅x + 1
        sage: unicode_art(list(DyckWords(3)))
        ⎡                                   ╱╲   ⎤
        ⎢            ╱╲    ╱╲      ╱╲╱╲    ╱  ╲  ⎥
        ⎣ ╱╲╱╲╱╲, ╱╲╱  ╲, ╱  ╲╱╲, ╱    ╲, ╱    ╲ ⎦
        sage: unicode_art(1)
        1
    """
    separator = kwds.pop('sep', empty_unicode_art)
    if kwds:
        raise ValueError('unknown keyword arguments: {0}'.format(list(kwds)))
    if len(obj) == 1:
        return _unicode_art_factory.build(obj[0])
    if not isinstance(separator, UnicodeArt):
        separator = _unicode_art_factory.build(separator)
    obj = map(_unicode_art_factory.build, obj)
    return _unicode_art_factory.concatenate(obj, separator, empty_unicode_art)
