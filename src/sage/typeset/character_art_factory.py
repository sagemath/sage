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
from six import iteritems, string_types, text_type, binary_type

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
        assert isinstance(string_type('a'), string_types)
        self.string_type = string_type
        assert magic_method_name in ['_ascii_art_', '_unicode_art_']
        self.magic_method_name = magic_method_name
        self.left_parenthesis, self.right_parenthesis = parenthesis
        self.left_square_bracket, self.right_square_bracket = square_bracet
        self.left_curly_brace, self.right_curly_brace = curly_brace

    def build(self, obj, baseline=None):
        r"""
        Construct a character art representation.

        INPUT:

        - ``obj`` -- anything; the object whose ascii art representation
          we want
        - ``baseline`` -- (optional) the baseline of the object

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
             /        _________    \
            -\2*x + \/ 1 - 4*x  - 1/
            -------------------------
                       _________
                 2*x*\/ 1 - 4*x
            sage: ascii_art(list(DyckWords(3)))
            [                                   /\   ]
            [            /\    /\      /\/\    /  \  ]
            [ /\/\/\, /\/  \, /  \/\, /    \, /    \ ]
            sage: ascii_art(1)
            1
        """
        if isinstance(obj, self.art_type):
            if baseline is not None:
                from copy import copy
                obj = copy(obj)
                obj._baseline = baseline
            return obj
        if isinstance(obj, SageObject):
            return self.build_from_magic_method(obj, baseline)
        if baseline is None:
            baseline = 0
        if isinstance(obj, tuple):
            return self.build_tuple(obj, baseline)
        elif isinstance(obj, dict):
            return self.build_dict(obj, baseline)
        elif isinstance(obj, list):
            return self.build_list(obj, baseline)
        elif isinstance(obj, set):
            return self.build_set(obj, baseline)
        else:
            return self.build_from_string(obj, baseline)

    def build_empty(self):
        """
        Return the empty character art object

        OUTPUT:

        Character art instance.

        EXAMPLES::

            sage: from sage.typeset.ascii_art import _ascii_art_factory as factory
            sage: str(factory.build_empty())
            ''
        """
        return self.art_type.empty()

    def build_from_magic_method(self, obj, baseline=None):
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
        ret = magic_method()
        if baseline is not None:
            ret._baseline = baseline
        return ret

    def build_from_string(self, obj, baseline=0):
        r"""
        Return the character art object created from splitting
        the object's string representation.

        INPUT:

        - ``obj`` -- utf-8 encoded byte string or unicode
        - ``baseline`` -- (default: 0) the baseline of the object

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

            sage: from sage.typeset.ascii_art import _ascii_art_factory as factory
            sage: factory.build_from_string(u'a\nbb\nccc')  # same with unicode
            a
            bb
            ccc

        ::

            sage: a = factory.build_from_string('a\nbb\nccc', baseline=2)
            sage: a + ascii_art('<-')
            a  <-
            bb
            ccc
        """
        if self.string_type is text_type and not isinstance(obj, text_type):
            if isinstance(obj, binary_type):
                obj = obj.decode('utf-8')
            else:
                obj = text_type(obj)
        if self.string_type is binary_type and not isinstance(obj, binary_type):
            obj = text_type(obj).encode('utf-8')
        return self.art_type(obj.splitlines(), baseline=baseline)

    def build_container(self, content, left_border, right_border, baseline=0):
        r"""
        Return character art for a container.

        INPUT:

        - ``content`` --
          :class:`~sage.typeset.character_art.CharacterArt`; the
          content of the container, usually comma-separated entries

        - ``left_border`` --
          :class:`~sage.typeset.symbols.CompoundSymbol`; the left
          border of the container

        - ``right_border`` --
          :class:`~sage.typeset.symbols.CompoundSymbol`; the right
          border of the container

        - ``baseline`` -- (default: 0) the baseline of the object

        TESTS::

            sage: l = ascii_art(list(DyckWords(3)))  # indirect doctest
            sage: l
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
        return self.art_type(lines, basepoints, baseline=baseline)

    def build_set(self, s, baseline=0):
        r"""
        Return an character art output of a set.

        TESTS:

        When the constructor is passed a set, this method is called.  Since
        iteration over sets is non-deterministic so too is the results of this
        test::

            sage: ascii_art(set(DyckWords(3)))  # indirect doctest random
            {                                   /\   }
            {  /\      /\/\              /\    /  \  }
            { /  \/\, /    \, /\/\/\, /\/  \, /    \ }

        We can also call this method directly an pass an iterable that is not a
        set, but still obtain the same output formatting::

            sage: from sage.typeset.ascii_art import _ascii_art_factory as factory
            sage: factory.build_set(sorted(set(DyckWords(3))))
            {                                   /\   }
            {            /\    /\      /\/\    /  \  }
            { /\/\/\, /\/  \, /  \/\, /    \, /    \ }
        """
        comma = self.art_type([self.string_type(', ')], baseline=0)
        repr_elems = self.concatenate(s, comma)
        return self.build_container(
            repr_elems, self.left_curly_brace, self.right_curly_brace,
            baseline)

    def build_dict(self, d, baseline=0):
        r"""
        Return an character art output of a dictionary.

        TESTS::

            sage: from collections import OrderedDict
            sage: d = OrderedDict(enumerate(DyckWords(3)))
            sage: art = ascii_art(d)  # indirect doctest
            sage: art
            {                                             /\   }
            {                /\      /\        /\/\      /  \  }
            { 0:/\/\/\, 1:/\/  \, 2:/  \/\, 3:/    \, 4:/    \ }
            sage: art.get_breakpoints()
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
                (concat_no_breakpoint(k, v) for k, v in iteritems(d)),
                comma)
        return self.build_container(
                repr_elems, self.left_curly_brace, self.right_curly_brace,
                baseline)

    def build_list(self, l, baseline=0):
        r"""
        Return an character art output of a list.

        TESTS::

            sage: l = ascii_art(list(DyckWords(3)))  # indirect doctest
            sage: l
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
            repr_elems, self.left_square_bracket, self.right_square_bracket,
            baseline)

    def build_tuple(self, t, baseline=0):
        r"""
        Return an character art output of a tuple.

        TESTS::

            sage: ascii_art(tuple(DyckWords(3)))  # indirect doctest
            (                                   /\   )
            (            /\    /\      /\/\    /  \  )
            ( /\/\/\, /\/  \, /  \/\, /    \, /    \ )
        """
        comma = self.art_type([self.string_type(', ')],
                              baseline=0,
                              breakpoints=[1])
        repr_elems = self.concatenate(t, comma)
        return self.build_container(
            repr_elems, self.left_parenthesis, self.right_parenthesis,
            baseline)

    def concatenate(self, iterable, separator, empty=None, baseline=0):
        """
        Concatenate multiple character art instances

        The breakpoints are set as the breakpoints of the ``separator``
        together with the breakpoints of the objects in ``iterable``.
        If there is ``None``, the end of the separator is used.

        INPUT:

        - ``iterable`` -- iterable of character art

        - ``separable`` -- character art; the separator in-between the
          iterable

        - ``empty`` -- an optional character art which is returned if
          ``iterable`` is empty

        - ``baseline`` -- (default: 0) the baseline of the object

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
        result._baseline = baseline
        return result

    def parse_keywords(self, kwds):
        """
        Parse the keyword input given by the dict ``kwds``.

        INPUT:

        - ``kwds`` -- a dict

        OUTPUT:

        A triple:

        - the separator
        - the baseline
        - the baseline of the separator

        .. WARNING::

            The input is a dict, not a list of keyword arguments.

        .. NOTE::

            This will remove ``sep``/``separator`` and ``baseline``
            from ``kwds`` if they are specified.

        TESTS::

            sage: from sage.typeset.ascii_art import _ascii_art_factory as factory
            sage: d = {'sep': '2', 'baseline': 5}
            sage: factory.parse_keywords(d)
            ('2', 5, None)
            sage: d
            {}

        ::

            sage: d = {'foo': '2', 'baseline': 5}
            sage: factory.parse_keywords(d)
            (, 5, None)
            sage: d
            {'foo': '2'}

            sage: d = {'sep': '2', 'separator': '2'}
            sage: factory.parse_keywords(d)
            Traceback (most recent call last):
            ...
            ValueError: cannot specify both 'sep' and 'separator'
        """
        empty = self.build_empty()
        sep = kwds.pop("sep", empty)
        if sep == empty:
            sep = kwds.pop("separator", empty)
        elif "separator" in kwds:
            raise ValueError("cannot specify both 'sep' and 'separator'")
        return sep, kwds.pop("baseline", None), kwds.pop("sep_baseline", None)
