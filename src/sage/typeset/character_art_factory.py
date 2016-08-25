# -*- coding: utf-8 -*-
r"""
Factory for Character-Based Art
"""
#*******************************************************************************
#       Copyright (C) 2013 Jean-Baptiste Priez <jbp@kerios.fr>,
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

from sage.structure.sage_object import SageObject

class CharacterArtFactory(SageObject):

    def __init__(self,
                 art_type, string_type, magic_method_name,
                 parenthesis, square_bracet, curly_brace):
        r"""
        Abstract base class for character art factory

        This class is the common implementation behind
        :func:`~sage.typeset.ascii_art.ascii_art` and
        :func:`~sage.typeset.unicode_art.unicode_art` .

        INPUT:

        - ``art_type`` -- type of the character art (i.e. a subclass of
          :class:`~sage.typeset.character_art.CharacterArt`)

        - ``string_type`` -- type of strings (the lines in the
          character art, e.g. ``str`` or ``unicode``).

        - ``magic_method_name`` -- name of the Sage magic method (e.g.
          ``'_ascii_art_'`` or ``'_unicode_art_'``).

        - ``parenthesis`` -- left/right pair of two multi-line
          symbols. The parenthesis, a.k.a. round brackets (used for printing
          tuples).

        - ``square_bracket`` -- left/right pair of two multi-line
          symbols. The square_brackets (used for printing lists).

        - ``curly_brace`` -- left/right pair of two multi-line
          symbols. The curly braces (used for printing sets).

        EXAMPLES::

            sage: from sage.typeset.ascii_art import _ascii_art_factory as factory
            sage: type(factory)
            <class 'sage.typeset.character_art_factory.CharacterArtFactory'>
        """
        self.art_type = art_type
        assert string_type in [str, unicode]
        self.string_type = string_type
        assert magic_method_name in ['_ascii_art_', '_unicode_art_']
        self.magic_method_name = magic_method_name
        self.left_parenthesis, self.right_parenthesis = parenthesis
        self.left_square_bracket, self.right_square_bracket = square_bracet
        self.left_curly_brace, self.right_curly_brace = curly_brace

    def build(self, obj):
        r"""
        Construct a character art representation

        INPUT:

        - ``obj`` -- anything. The object whose ascii art representation
          we want.

        OUTPUT:

        Character art object.

        EXAMPLES::

            sage: ascii_art(integral(exp(x+x^2)/(x+1), x))
                /
               |
               |   2
               |  x  + x
               | e
               | ------- dx
               |  x + 1
               |
              /

        TESTS::

            sage: n = var('n')
            sage: ascii_art(sum(binomial(2 * n, n + 1) * x^n, n, 0, oo))
             /        __________    \
            -\2*x + \/ -4*x + 1  - 1/
            --------------------------
                       __________
                 2*x*\/ -4*x + 1
            sage: ascii_art(list(DyckWords(3)))
            [                                   /\   ]
            [            /\    /\      /\/\    /  \  ]
            [ /\/\/\, /\/  \, /  \/\, /    \, /    \ ]
            sage: ascii_art(1)
            1
        """
        if isinstance(obj, self.art_type):
            return obj
        elif isinstance(obj, (tuple, list, dict, set)):
            if obj.__class__ is tuple:
                return self.build_tuple(obj)
            elif obj.__class__ is dict:
                return self.build_dict(obj)
            elif obj.__class__ is list:
                return self.build_list(obj)
            else:
                return self.build_set(obj)
        elif isinstance(obj, SageObject):
            return self.build_from_magic_method(obj)
        else:
            return self.build_from_string(obj)

    def build_empty(self):
        """
        Return the empty character art object

        OUTPUT:

        Character art instance.
        """
        return self.art_type.empty()

    def build_from_magic_method(self, obj):
        """
        Return the character art object created by the object's magic method

        OUTPUT:

        Character art instance.

        EXAMPLES::

            sage: from sage.typeset.ascii_art import _ascii_art_factory as factory
            sage: out = factory.build_from_magic_method(identity_matrix(2));  out
            [1 0]
            [0 1]
            sage: type(out)
            <class 'sage.typeset.ascii_art.AsciiArt'>
        """
        magic_method = getattr(obj, self.magic_method_name)
        return magic_method()

    def build_from_string(self, obj):
        r"""
        Return the character art object created from splitting the object's string representation

        INPUT:

        - ``obj`` -- utf-8 encoded byte string or unicode.

        OUTPUT:

        Character art instance.

        EXAMPLES::

            sage: from sage.typeset.ascii_art import _ascii_art_factory as factory
            sage: out = factory.build_from_string('a\nbb\nccc')
            sage: out + out + out
            a  a  a
            bb bb bb
            ccccccccc
            sage: type(out)
            <class 'sage.typeset.ascii_art.AsciiArt'>

        TESTS::

            sage: factory.build_from_string(u'a\nbb\nccc')  # same with unicode
            a
            bb
            ccc
        """
        if self.string_type is unicode and not isinstance(obj, unicode):
            obj = str(obj).decode('utf-8')
        if self.string_type is str and not isinstance(obj, str):
            obj = unicode(obj).encode('utf-8')
        return self.art_type(obj.splitlines())

    def build_container(self, content, left_border, right_border):
        r"""
        Return character art for a container

        INPUT:

        - ``content`` --
          :class:`~sage.typeset.character_art.CharacterArt`. The
          content of the container, usually comma-separated entries.

        - ``left_border`` --
          :class:`~sage.typeset.symbols.CompoundSymbol`. The left
          border of the container.

        - ``right_border`` --
          :class:`~sage.typeset.symbols.CompoundSymbol`. The right
          border of the container.

        TESTS::

            sage: l = ascii_art(list(DyckWords(3)));  l
            [                                   /\   ]
            [            /\    /\      /\/\    /  \  ]
            [ /\/\/\, /\/  \, /  \/\, /    \, /    \ ]
            sage: l.get_breakpoints()
            [9, 17, 25, 33]
        """
        w = content.width()
        h = content.height()
        left_border = left_border.character_art(h)
        right_border = right_border.character_art(h)
        lines = []
        pad = self.string_type(' ')
        for left, line, right in zip(left_border, content, right_border):
            lines.append(left + pad + line.ljust(w) + pad + right)
        shift = len(left_border) + len(pad)
        basepoints = [bp + shift for bp in content.get_breakpoints()]
        return self.art_type(lines, basepoints, baseline=0)

    def build_set(self, s):
        r"""
        Return an character art output of a set.

        TESTS::

            sage: ascii_art(set(DyckWords(3)))
            {                                   /\   }
            {  /\      /\/\              /\    /  \  }
            { /  \/\, /    \, /\/\/\, /\/  \, /    \ }
        """
        comma = self.art_type([self.string_type(', ')], baseline=0)
        repr_elems = self.concatenate(s, comma)
        return self.build_container(
            repr_elems, self.left_curly_brace, self.right_curly_brace)

    def build_dict(self, d):
        r"""
        Return an character art output of a dictionary.

        TESTS::

            sage: d = ascii_art({i:dw for i,dw in enumerate(DyckWords(3))})
            sage: d
            {                                             /\   }
            {                /\      /\        /\/\      /  \  }
            { 0:/\/\/\, 1:/\/  \, 2:/  \/\, 3:/    \, 4:/    \ }
            sage: d.get_breakpoints()
            [11, 21, 31, 41]
        """
        comma = self.art_type([self.string_type(', ')],
                              baseline=0,
                              breakpoints=[1])
        colon = self.art_type([self.string_type(':')], baseline=0)
        def concat_no_breakpoint(k,v):
            k = self.build(k)
            v = self.build(v)
            elt = k + colon + v
            elt._breakpoints.remove(k._l)
            elt._breakpoints.remove(k._l + 1)
            return elt
        repr_elems = self.concatenate(
                (concat_no_breakpoint(k,v) for k,v in d.iteritems()),
                comma)
        return self.build_container(repr_elems,
                self.left_curly_brace, self.right_curly_brace)

    def build_list(self, l):
        r"""
        Return an character art output of a list.

        TESTS::

            sage: l = ascii_art(list(DyckWords(3)));  l
            [                                   /\   ]
            [            /\    /\      /\/\    /  \  ]
            [ /\/\/\, /\/  \, /  \/\, /    \, /    \ ]
            sage: l.get_breakpoints()
            [9, 17, 25, 33]

        The breakpoints of the object are used as breakpoints::

            sage: l = ascii_art([DyckWords(2).list(), DyckWords(2).list()])
            sage: l.get_breakpoints()
            [9, 17, 25]
        """
        comma = self.art_type([self.string_type(', ')],
                              baseline=0,
                              breakpoints=[1])
        repr_elems = self.concatenate(l, comma)
        return self.build_container(
            repr_elems, self.left_square_bracket, self.right_square_bracket)

    def build_tuple(self, t):
        r"""
        Return an character art output of a tuple.

        TESTS::

            sage: ascii_art(tuple(DyckWords(3)))
            (                                   /\   )
            (            /\    /\      /\/\    /  \  )
            ( /\/\/\, /\/  \, /  \/\, /    \, /    \ )
        """
        comma = self.art_type([self.string_type(', ')],
                              baseline=0,
                              breakpoints=[1])
        repr_elems = self.concatenate(t, comma)
        return self.build_container(
            repr_elems, self.left_parenthesis, self.right_parenthesis)

    def concatenate(self, iterable, separator, empty=None):
        """
        Concatenate multiple character art instances

        The breakpoints are set as the breakpoints of the ``separator`` together
        with the breakpoints of the objects in ``iterable``. If there is
        ``None``, the end of the separator is used.

        INPUT:

        - ``iterable`` -- iterable of character art.

        - ``separable`` -- character art. The separator in-between the
          iterable.

        - ``empty`` -- an optional character art which is returned if
          ``iterable`` is empty

        EXAMPLES::

            sage: i2 = identity_matrix(2)
            sage: ascii_art(i2, i2, i2, sep=ascii_art(1/x))
                 1     1
            [1 0]-[1 0]-[1 0]
            [0 1]x[0 1]x[0 1]
        """
        if empty is None:
            empty = self.build_empty()
        result = empty
        breakpoints = []
        separator = self.build(separator)
        bk = separator.get_breakpoints()
        if not bk:
            bk = [separator._l]
        for obj in iterable:
            if result is not empty:
                l = result._l
                result += separator
                breakpoints.extend([l+x for x in bk])
            l = result._l
            obj = self.build(obj)
            result += self.build(obj)
            breakpoints.extend([l+x for x in obj.get_breakpoints()])
        result._breakpoints = breakpoints
        return result
