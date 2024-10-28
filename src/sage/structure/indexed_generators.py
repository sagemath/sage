# sage_setup: distribution = sagemath-objects
"""
Indexed Generators
"""
# ****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.structure.category_object import normalize_names


class IndexedGenerators():
    r"""nodetex
    Abstract base class for parents whose elements consist of generators
    indexed by an arbitrary set.

    Options controlling the printing of elements:

    - ``prefix`` -- string, prefix used for printing elements of this
      module (default: ``'x'``).  With the default, a monomial
      indexed by 'a' would be printed as ``x['a']``.

    - ``latex_prefix`` -- string or ``None``, prefix used in the `\LaTeX`
      representation of elements (default: ``None``); if this is
      anything except the empty string, it prints the index as a
      subscript.  If this is ``None``, it uses the setting for ``prefix``,
      so if ``prefix`` is set to "B", then a monomial indexed by 'a'
      would be printed as ``B_{a}``.  If this is the empty string, then
      don't print monomials as subscripts: the monomial indexed by 'a'
      would be printed as ``a``, or as ``[a]`` if ``latex_bracket`` is
      ``True``.

    - ``names`` -- dictionary with strings as values or list of strings (optional);
      a mapping from the indices of the generators to strings giving the
      generators explicit names. This is used instead of the print options
      ``prefix`` and ``bracket`` when ``names`` is specified.

    - ``latex_names`` -- dictionary with strings as values or list of strings
      (optional); same as ``names`` except using the `\LaTeX` representation

    - ``bracket`` -- ``None``, boolean, string, or list or tuple of
      strings (default: ``None``); if ``None``, use the value of the
      attribute ``self._repr_option_bracket``, which has default value
      ``True``.  (``self._repr_option_bracket`` is available for backwards
      compatibility.  Users should set ``bracket`` instead.  If
      ``bracket`` is set to anything except ``None``, it overrides
      the value of ``self._repr_option_bracket``.)  If ``False``, do not
      include brackets when printing elements: a monomial indexed by
      'a' would be printed as ``B'a'``, and a monomial indexed by
      (1,2,3) would be printed as ``B(1,2,3)``.  If ``True``, use "[" and
      "]" as brackets.  If it is one of "[", "(", or "{", use it and
      its partner as brackets.  If it is any other string, use it as
      both brackets.  If it is a list or tuple of strings, use the
      first entry as the left bracket and the second entry as the
      right bracket.

    - ``latex_bracket`` -- boolean, string, or list or tuple of strings
      (default: ``False``); if ``False``, do not include brackets in
      the LaTeX representation of elements.  This option is only
      relevant if ``latex_prefix`` is the empty string; otherwise,
      brackets are not used regardless.  If ``True``, use "\left[" and
      "\right]" as brackets.  If this is one of "[", "(", "\\{", "|",
      or "||", use it and its partner, prepended with "\left" and
      "\right", as brackets.  If this is any other string, use it as
      both brackets.  If this is a list or tuple of strings, use the
      first entry as the left bracket and the second entry as the
      right bracket.

    - ``scalar_mult`` -- string to use for scalar multiplication in
      the print representation (default: ``'*'``)

    - ``latex_scalar_mult`` -- string or ``None`` (default: ``None``);
      string to use for scalar multiplication in the latex
      representation.  If ``None``, use the empty string if ``scalar_mult``
      is set to "*", otherwise use the value of ``scalar_mult``.

    - ``tensor_symbol`` -- string or ``None`` (default: ``None``);
      string to use for tensor product in the print representation. If
      ``None``, use  ``sage.categories.tensor.symbol`` and
      ``sage.categories.tensor.unicode_symbol``.

    - ``sorting_key`` -- a key function (default: ``lambda x: x``);
      to use for sorting elements in the output of elements

    - ``sorting_reverse`` -- boolean (default: ``False``); if ``True``
      sort elements in reverse order in the output of elements

    - ``string_quotes`` -- boolean (default: ``True``); if ``True`` then
      display string indices with quotes

    - ``iterate_key`` -- boolean (default: ``False``); iterate through
      the elements of the key and print the result as comma separated
      objects for string output

    .. NOTE::

        These print options may also be accessed and modified using the
        :meth:`print_options` method, after the parent has been defined.

    EXAMPLES:

    We demonstrate a variety of the input options::

        sage: from sage.structure.indexed_generators import IndexedGenerators
        sage: I = IndexedGenerators(ZZ, prefix='A')
        sage: I._repr_generator(2)
        'A[2]'
        sage: I._latex_generator(2)
        'A_{2}'

        sage: I = IndexedGenerators(ZZ, bracket='(')
        sage: I._repr_generator(2)
        'x(2)'
        sage: I._latex_generator(2)
        'x_{2}'

        sage: I = IndexedGenerators(ZZ, prefix='', latex_bracket='(')
        sage: I._repr_generator(2)
        '[2]'
        sage: I._latex_generator(2)
        \left( 2 \right)

        sage: I = IndexedGenerators(ZZ, bracket=['|', '>'])
        sage: I._repr_generator(2)
        'x|2>'
    """
    def __init__(self, indices, prefix='x', **kwds):
        """
        Initialize ``self``.

        EXAMPLES:

        This is a mixin class, so don't need pickling equality::

            sage: I = sage.structure.indexed_generators.IndexedGenerators(ZZ)
            sage: TestSuite(I).run(skip='_test_pickling')
        """
        self._indices = indices

        # printing options for elements (set when initializing self).
        # This includes self._repr_option_bracket (kept for backwards
        # compatibility, declared to be True by default, needs to be
        # overridden explicitly).
        self._print_options = {'prefix': prefix,
                               'names': None,
                               'bracket': None,
                               'latex_bracket': False,
                               'latex_prefix': None,
                               'latex_names': None,
                               'scalar_mult': "*",
                               'latex_scalar_mult': None,
                               'tensor_symbol': None,
                               'string_quotes': True,
                               'sorting_key': lambda x: x,
                               'sorting_reverse': False,
                               'iterate_key': False}
        # 'bracket': its default value here is None, meaning that
        # the value of self._repr_option_bracket is used; the default
        # value of that attribute is True -- see immediately before
        # the method _repr_generator.  If 'bracket' is any value
        # except None, then it overrides the value of
        # self._repr_option_bracket.  Future users might consider
        # using 'bracket' instead of _repr_option_bracket.
        self.print_options(**kwds)

    def indices(self):
        """
        Return the indices of ``self``.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])                      # needs sage.modules
            sage: F.indices()                                                           # needs sage.modules
            {'a', 'b', 'c'}
        """
        return self._indices

    def prefix(self):
        """
        Return the prefix used when displaying elements of ``self``.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])                      # needs sage.modules
            sage: F.prefix()                                                            # needs sage.modules
            'B'

        ::

            sage: X = SchubertPolynomialRing(QQ)                                        # needs sage.combinat sage.modules
            sage: X.prefix()                                                            # needs sage.combinat sage.modules
            'X'
        """
        return self._print_options['prefix']

    def print_options(self, **kwds):
        """
        Return the current print options, or set an option.

        INPUT:

        All of the input is optional; if present, it should be
        in the form of keyword pairs, such as
        ``latex_bracket='('``.  The allowable keywords are:

        - ``prefix``
        - ``latex_prefix``
        - ``names``
        - ``latex_names``
        - ``bracket``
        - ``latex_bracket``
        - ``scalar_mult``
        - ``latex_scalar_mult``
        - ``tensor_symbol``
        - ``string_quotes``
        - ``sorting_key``
        - ``sorting_reverse``
        - ``iterate_key``

        See the documentation for :class:`IndexedGenerators` for
        descriptions of the effects of setting each of these options.

        OUTPUT: if the user provides any input, set the appropriate
        option(s) and return nothing.  Otherwise, return the
        dictionary of settings for print and LaTeX representations.

        EXAMPLES::

            sage: # needs sage.modules
            sage: F = CombinatorialFreeModule(ZZ, [1,2,3], prefix='x')
            sage: F.print_options()
            {...'prefix': 'x'...}
            sage: F.print_options(bracket='(')
            sage: F.print_options()
            {...'bracket': '('...}

        TESTS::

            sage: sorted(F.print_options().items())                                     # needs sage.modules
            [('bracket', '('), ('iterate_key', False),
             ('latex_bracket', False), ('latex_names', None),
             ('latex_prefix', None), ('latex_scalar_mult', None),
             ('names', None), ('prefix', 'x'),
             ('scalar_mult', '*'),
             ('sorting_key', <function ...<lambda> at ...>),
             ('sorting_reverse', False), ('string_quotes', True),
             ('tensor_symbol', None)]
            sage: F.print_options(bracket='[')  # reset                                 # needs sage.modules
        """
        # don't just use kwds.get(...) because I want to distinguish
        # between an argument like "option=None" and the option not
        # being there altogether.
        if kwds:
            for option in kwds:
                if option in self._print_options:
                    self._print_options[option] = kwds[option]
                else:
                    raise ValueError('{} is not a valid print option.'.format(option))
            return
        return self._print_options

    _repr_option_bracket = True

    def _parse_names(self, m, use_latex):
        """
        Return the name corresponding to ``m`` if it exists,
        otherwise return ``None``.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(ZZ, [1,2,3], names='a,b,c',               # needs sage.modules
            ....:                             latex_names='x,y,z')
            sage: F._parse_names(1, False)                                              # needs sage.modules
            'a'
            sage: F._parse_names(1, True)                                               # needs sage.modules
            'x'

            sage: F.print_options(latex_names=None)                                     # needs sage.modules
            sage: F._parse_names(1, True)                                               # needs sage.modules
            'a'

            sage: # needs sage.modules
            sage: F.print_options(latex_names={1:'x', 2:'y'}, names=None)
            sage: F._parse_names(1, False) is None
            True
            sage: F._parse_names(1, True)
            'x'
            sage: F._parse_names(3, True) is None
            True

            sage: # needs sage.modules
            sage: F.print_options(names={1:'a', 3:'c'}, latex_names=None)
            sage: F._parse_names(1, False)
            'a'
            sage: F._parse_names(1, True)
            'a'
            sage: F._parse_names(2, False) is None
            True
            sage: F._parse_names(2, True) is None
            True
        """
        names = self._print_options.get('names', None)
        if use_latex:
            latex_names = self._print_options.get('latex_names', None)
            if latex_names is not None:
                names = latex_names

        if names is not None:
            if isinstance(names, dict):
                try:
                    return names[m]
                except KeyError:
                    return None
            else: # treat it like a list
                try:
                    i = self._indices.rank(m)
                except (AttributeError, TypeError, KeyError, ValueError):
                    return None
                if i >= len(names):
                    return None
                try:
                    return names[i]
                except (AttributeError, TypeError, KeyError, ValueError):
                    return None
        return None

    def _repr_generator(self, m):
        """
        Return a string representing the generator indexed by ``m``.

        The output can be customized by setting any of the following
        options when initializing the parent:

        - ``prefix``
        - ``bracket``
        - ``scalar_mult``
        - ``names``
        - ``iterate_key``

        Alternatively, one can use the :meth:`print_options` method
        to achieve the same effect.  To modify the bracket setting,
        one can also set ``self._repr_option_bracket`` as long as one
        has *not* set the ``bracket`` option: if the
        ``bracket`` option is anything but ``None``, it overrides
        the value of ``self._repr_option_bracket``.

        See the documentation for :class:`CombinatorialFreeModule` for
        details on the initialization options.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])                      # needs sage.modules
            sage: e = F.basis()                                                         # needs sage.modules
            sage: e['a'] + 2*e['b']    # indirect doctest                               # needs sage.modules
            B['a'] + 2*B['b']

            sage: # needs sage.modules
            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'], prefix='F')
            sage: e = F.basis()
            sage: e['a'] + 2*e['b']    # indirect doctest
            F['a'] + 2*F['b']
            sage: F.print_options(string_quotes=False)
            sage: e['a'] + 2*e['b']
            F[a] + 2*F[b]

            sage: # needs sage.modules
            sage: F = CombinatorialFreeModule(QQ, ['aa', 'bb', 'cc'], prefix='F')
            sage: e = F.basis()
            sage: F.print_options(iterate_key=True)
            sage: e['aa'] + 2*e['bb']
            F['a', 'a'] + 2*F['b', 'b']
            sage: F.print_options(string_quotes=False)
            sage: e['aa'] + 2*e['bb']
            F[a, a] + 2*F[b, b]

            sage: # needs sage.combinat sage.modules
            sage: QS3 = CombinatorialFreeModule(QQ, Permutations(3), prefix="")
            sage: original_print_options = QS3.print_options()
            sage: a = 2*QS3([1,2,3])+4*QS3([3,2,1])
            sage: a                      # indirect doctest
            2*[[1, 2, 3]] + 4*[[3, 2, 1]]

            sage: QS3.print_options(bracket = False)                                    # needs sage.combinat sage.modules
            sage: a              # indirect doctest                                     # needs sage.combinat sage.modules
            2*[1, 2, 3] + 4*[3, 2, 1]

            sage: QS3.print_options(prefix='')                                          # needs sage.combinat sage.modules
            sage: a              # indirect doctest                                     # needs sage.combinat sage.modules
            2*[1, 2, 3] + 4*[3, 2, 1]

            sage: QS3.print_options(bracket="|", scalar_mult=" *@* ")                   # needs sage.combinat sage.modules
            sage: a              # indirect doctest                                     # needs sage.combinat sage.modules
            2 *@* |[1, 2, 3]| + 4 *@* |[3, 2, 1]|

            sage: QS3.print_options(bracket="|", scalar_mult="*", iterate_key=True)     # needs sage.combinat sage.modules
            sage: a              # indirect doctest                                     # needs sage.combinat sage.modules
            2*|1, 2, 3| + 4*|3, 2, 1|

            sage: QS3.print_options(**original_print_options)  # reset                  # needs sage.combinat sage.modules

        TESTS::

            sage: F = CombinatorialFreeModule(QQ, [('a', 'b'), ('c','d')])              # needs sage.modules
            sage: e = F.basis()                                                         # needs sage.modules
            sage: e[('a','b')] + 2*e[('c','d')]    # indirect doctest                   # needs sage.modules
            B[('a', 'b')] + 2*B[('c', 'd')]

            sage: F.<a,b,c> = CombinatorialFreeModule(QQ)                               # needs sage.modules
            sage: a + 2*b                                                               # needs sage.modules
            a + 2*b

            sage: # needs sage.modules
            sage: F = CombinatorialFreeModule(QQ, ZZ)
            sage: e = F.basis()
            sage: 3*e[1] + 2*e[-2]
            2*B[-2] + 3*B[1]
            sage: F.print_options(iterate_key=True)
            sage: 3*e[1] + 2*e[-2]
            2*B[-2] + 3*B[1]
        """
        ret = self._parse_names(m, False)
        if ret is not None:
            return ret
        bracket = self._print_options.get('bracket', None)
        bracket_d = {"{": "}", "[": "]", "(": ")"}
        if bracket is None:
            bracket = self._repr_option_bracket
        if bracket is True:
            left = "["
            right = "]"
        elif bracket is False:
            left = ""
            right = ""
        elif isinstance(bracket, (tuple, list)):
            left = bracket[0]
            right = bracket[1]
        elif bracket in bracket_d:
            left = bracket
            right = bracket_d[bracket]
        else:
            left = bracket
            right = bracket
        iterate_key = self._print_options.get('iterate_key', False)
        quotes = self._print_options.get('string_quotes', True)
        if iterate_key:
            try:
                m = iter(m)
            except TypeError:
                pass  # not iterable, so fallback to normal behavior
            else:
                if not quotes:
                    return self.prefix() + left + (', '.join(str(val) for val in m)) + right
                return self.prefix() + left + (', '.join(repr(val) for val in m)) + right
        if not quotes and isinstance(m, str):
            return self.prefix() + left + m + right
        return self.prefix() + left + repr(m) + right # mind the (m), to accept a tuple for m

    def _ascii_art_generator(self, m):
        r"""
        Return an ascii art representing the generator indexed by ``m``.

        TESTS::

            sage: # needs sage.combinat sage.modules
            sage: R = NonCommutativeSymmetricFunctions(QQ).R()
            sage: ascii_art(R[1,2,2,4])
            R
               ****
              **
             **
             *
            sage: Partitions.options(diagram_str='#', convention='french')
            sage: ascii_art(R[1,2,2,4])
            R
             #
             ##
              ##
               ####
            sage: Partitions.options._reset()

            sage: F.<a,b,c> = CombinatorialFreeModule(QQ)                               # needs sage.modules
            sage: ascii_art(a + 2*b)                                                    # needs sage.modules
            a + 2*b
        """
        from sage.typeset.ascii_art import AsciiArt, ascii_art
        ret = self._parse_names(m, False)
        if ret is not None:
            return ascii_art(ret)

        pref = AsciiArt([self.prefix()])
        if not pref:
            return ascii_art(m)
        r = pref * (AsciiArt([" " * len(pref)]) + ascii_art(m))
        r._baseline = r._h - 1
        return r

    def _unicode_art_generator(self, m):
        r"""
        Return an unicode art representing the generator indexed by ``m``.

        TESTS::

            sage: # needs sage.combinat
            sage: R = NonCommutativeSymmetricFunctions(QQ).R()                          # needs sage.modules
            sage: unicode_art(R[1,2,2,4])                                               # needs sage.modules
            R
               ┌┬┬┬┐
              ┌┼┼┴┴┘
             ┌┼┼┘
             ├┼┘
             └┘
            sage: Partitions.options.convention="french"
            sage: unicode_art(R[1,2,2,4])                                               # needs sage.modules
            R
             ┌┐
             ├┼┐
             └┼┼┐
              └┼┼┬┬┐
               └┴┴┴┘
            sage: Partitions.options._reset()

            sage: F.<a,b,c> = CombinatorialFreeModule(QQ)                               # needs sage.modules
            sage: unicode_art(a + 2*b)                                                  # needs sage.modules
            a + 2*b
        """
        from sage.typeset.unicode_art import UnicodeArt, unicode_art
        ret = self._parse_names(m, False)
        if ret is not None:
            return unicode_art(ret)

        pref = UnicodeArt([self.prefix()])
        if not pref:
            return unicode_art(m)
        r = pref * (UnicodeArt([" " * len(pref)]) + unicode_art(m))
        r._baseline = r._h - 1
        return r

    def _latex_generator(self, m):
        r"""
        Return a `\LaTeX` for the generator indexed by ``m``.

        The output can be customized by setting any of the following
        options when initializing the parent:

        - ``prefix``
        - ``latex_prefix``
        - ``latex_bracket``
        - ``names``
        - ``latex_names``

        (Alternatively, one can use the :meth:`print_options` method
        to achieve the same effect.)

        See the documentation for :class:`CombinatorialFreeModule` for
        details on the initialization options.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])                      # needs sage.modules
            sage: e = F.basis()                                                         # needs sage.modules
            sage: latex(e['a'] + 2*e['b'])    # indirect doctest                        # needs sage.modules
            B_{a} + 2 B_{b}

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'], prefix='C')          # needs sage.modules
            sage: e = F.basis()                                                         # needs sage.modules
            sage: latex(e['a'] + 2*e['b'])    # indirect doctest                        # needs sage.modules
            C_{a} + 2 C_{b}

            sage: # needs sage.combinat sage.modules
            sage: QS3 = CombinatorialFreeModule(QQ, Permutations(3),
            ....:                               prefix='', scalar_mult='*')
            sage: original_print_options = QS3.print_options()
            sage: a = 2*QS3([1,2,3])+4*QS3([3,2,1])
            sage: latex(a)                     # indirect doctest
            2 [1, 2, 3] + 4 [3, 2, 1]
            sage: QS3.print_options(latex_bracket=True)
            sage: latex(a)                     # indirect doctest
            2 \left[ [1, 2, 3] \right] + 4 \left[ [3, 2, 1] \right]
            sage: QS3.print_options(latex_bracket="(")
            sage: latex(a)                     # indirect doctest
            2 \left( [1, 2, 3] \right) + 4 \left( [3, 2, 1] \right)
            sage: QS3.print_options(latex_bracket=('\\myleftbracket',
            ....:                                  '\\myrightbracket'))
            sage: latex(a)                     # indirect doctest
            2 \myleftbracket [1, 2, 3] \myrightbracket + 4 \myleftbracket [3, 2, 1] \myrightbracket
            sage: QS3.print_options(**original_print_options) # reset

        TESTS::

            sage: # needs sage.modules
            sage: F = CombinatorialFreeModule(QQ, [('a', 'b'), (0,1,2)])
            sage: e = F.basis()
            sage: latex(e[('a','b')])    # indirect doctest
            B_{('a', 'b')}
            sage: latex(2*e[(0,1,2)])    # indirect doctest
            2 B_{\left(0, 1, 2\right)}
            sage: F = CombinatorialFreeModule(QQ, [('a', 'b'), (0,1,2)], prefix="")
            sage: e = F.basis()
            sage: latex(2*e[(0,1,2)])    # indirect doctest
            2 \left(0, 1, 2\right)

            sage: F.<a,b,c> = CombinatorialFreeModule(QQ, latex_names='x,y,z')          # needs sage.modules
            sage: latex(a + 2*b)                                                        # needs sage.modules
            x + 2 y
        """
        from sage.misc.latex import latex

        ret = self._parse_names(m, True)
        if ret is not None:
            return ret

        s = latex(m)
        if s.find('\\text{\\textt') != -1:
            # m contains "non-LaTeXed" strings, use string representation
            s = str(m)

        # dictionary with left-right pairs of "brackets".  put pairs
        # in here accept \\left and \\right as prefixes.
        bracket_d = {"{": "\\}", "[": "]", "(": ")", "\\{": "\\}",
                     "|": "|", "||": "||"}
        bracket = self._print_options.get('latex_bracket', False)
        if bracket is True:
            left = "\\left["
            right = "\\right]"
        elif bracket is False:
            left = ""
            right = ""
        elif isinstance(bracket, (tuple, list)):
            left = bracket[0]
            right = bracket[1]
        elif bracket in bracket_d:
            left = bracket
            right = bracket_d[bracket]
            if left == "{":
                left = "\\{"
            left = "\\left" + left
            right = "\\right" + right
        else:
            left = bracket
            right = bracket
        prefix = self._print_options.get('latex_prefix')
        if prefix is None:
            prefix = self._print_options.get('prefix')
        if prefix == "":
            return left + s + right
        return "%s_{%s}" % (prefix, s)


def split_index_keywords(kwds):
    """
    Split the dictionary ``kwds`` into two dictionaries, one containing
    keywords for :class:`IndexedGenerators`, and the other is everything else.

    OUTPUT:

    The dictionary containing only they keywords
    for :class:`IndexedGenerators`. This modifies the dictionary ``kwds``.

    .. WARNING::

        This modifies the input dictionary ``kwds``.

    EXAMPLES::

        sage: from sage.structure.indexed_generators import split_index_keywords
        sage: d = {'string_quotes': False, 'bracket': None, 'base': QQ}
        sage: split_index_keywords(d)
        {'bracket': None, 'string_quotes': False}
        sage: d
        {'base': Rational Field}
    """
    ret = {}
    for option in ['prefix', 'latex_prefix', 'bracket', 'latex_bracket',
                   'scalar_mult', 'latex_scalar_mult', 'tensor_symbol',
                   'sorting_key', 'sorting_reverse',
                   'string_quotes']:
        try:
            ret[option] = kwds.pop(option)
        except KeyError:
            pass
    return ret


def parse_indices_names(names, index_set, prefix, kwds=None):
    """
    Parse the names, index set, and prefix input, along with setting
    default values for keyword arguments ``kwds``.

    OUTPUT:

    The triple ``(N, I, p)``:

    - ``N`` is the tuple of variable names,
    - ``I`` is the index set, and
    - ``p`` is the prefix.

    This modifies the dictionary ``kwds``.

    .. NOTE::

        When the indices, names, or prefix have not been given, it
        should be passed to this function as ``None``.

    .. NOTE::

        For handling default prefixes, if the result will be ``None`` if
        it is not processed in this function.

    EXAMPLES::

        sage: from sage.structure.indexed_generators import parse_indices_names
        sage: d = {}
        sage: parse_indices_names('x,y,z', ZZ, None, d)
        (('x', 'y', 'z'), Integer Ring, None)
        sage: d
        {}
        sage: d = {}
        sage: parse_indices_names('x,y,z', None, None, d)
        (('x', 'y', 'z'), {'x', 'y', 'z'}, '')
        sage: d
        {'string_quotes': False}
        sage: d = {}
        sage: parse_indices_names(None, ZZ, None, d)
        (None, Integer Ring, None)
        sage: d
        {}

    ::

        sage: d = {'string_quotes':True, 'bracket':'['}
        sage: parse_indices_names(['a','b','c'], ZZ, 'x', d)
        (('a', 'b', 'c'), Integer Ring, 'x')
        sage: d
        {'bracket': '[', 'string_quotes': True}
        sage: parse_indices_names('x,y,z', None, 'A', d)
        (('x', 'y', 'z'), {'x', 'y', 'z'}, 'A')
        sage: d
        {'bracket': '[', 'string_quotes': True}
    """
    if index_set is None:
        if names is None:
            raise ValueError("either the indices or names must be given")

        if prefix is None:
            prefix = ''
        if kwds is None:
            kwds = {}
        kwds.setdefault('string_quotes', False)

    names, index_set = standardize_names_index_set(names, index_set, -1)

    return (names, index_set, prefix)


def standardize_names_index_set(names=None, index_set=None, ngens=None):
    """
    Standardize the ``names`` and ``index_set`` inputs.

    INPUT:

    - ``names`` -- (optional) the variable names
    - ``index_set`` -- (optional) the index set
    - ``ngens`` -- (optional) the number of generators

    If ``ngens`` is a negative number, then this does not check that
    the number of variable names matches the size of the index set.

    OUTPUT:

    A pair ``(names_std, index_set_std)``, where ``names_std`` is either
    ``None`` or a tuple of strings, and where ``index_set_std`` is a finite
    enumerated set.
    The purpose of ``index_set_std`` is to index the generators of some object
    (e.g., the basis of a module); the strings in ``names_std``, when they
    exist, are used for printing these indices. The ``ngens``

    If ``names`` contains exactly one name ``X`` and ``ngens`` is greater than
    1, then ``names_std`` are ``Xi`` for ``i`` in ``range(ngens)``.

    TESTS::

        sage: from sage.structure.indexed_generators import standardize_names_index_set
        sage: standardize_names_index_set('x,y')
        (('x', 'y'), {'x', 'y'})
        sage: standardize_names_index_set(['x','y'])
        (('x', 'y'), {'x', 'y'})
        sage: standardize_names_index_set(['x','y'], ['a','b'])
        (('x', 'y'), {'a', 'b'})
        sage: standardize_names_index_set('x,y', ngens=2)
        (('x', 'y'), {'x', 'y'})
        sage: standardize_names_index_set(index_set=['a','b'], ngens=2)
        (None, {'a', 'b'})
        sage: standardize_names_index_set('x', ngens=3)
        (('x0', 'x1', 'x2'), {'x0', 'x1', 'x2'})

        sage: standardize_names_index_set()
        Traceback (most recent call last):
        ...
        ValueError: the index_set, names, or number of generators must be specified
        sage: standardize_names_index_set(['x'], ['a', 'b'])
        Traceback (most recent call last):
        ...
        IndexError: the number of names must equal the size of the indexing set
        sage: standardize_names_index_set('x,y', ['a'])
        Traceback (most recent call last):
        ...
        IndexError: the number of names must equal the size of the indexing set
        sage: standardize_names_index_set('x,y,z', ngens=2)
        Traceback (most recent call last):
        ...
        IndexError: the number of names must equal the number of generators
        sage: standardize_names_index_set(index_set=['a'], ngens=2)
        Traceback (most recent call last):
        ...
        IndexError: the size of the indexing set must equal the number of generators
    """
    if names is not None:
        if ngens is None or ngens < 0:
            names = normalize_names(-1, names)
        else:
            names = normalize_names(ngens, names)

    if index_set is None:
        if names is None:
            # If neither is specified, we make range(ngens) the index set
            if ngens is None:
                raise ValueError("the index_set, names, or number of"
                                 " generators must be specified")
            index_set = tuple(range(ngens))
        else:
            # If only the names are specified, then we make the indexing set
            #   be the names
            index_set = tuple(names)

    from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
    if isinstance(index_set, dict):  # dict of {name: index} -- not likely to be used
        if names is not None:
            raise ValueError("cannot give index_set as a dict and names")
        names = normalize_names(-1, tuple(index_set.keys()))
        index_set = FiniteEnumeratedSet([index_set[n] for n in names])
    elif isinstance(index_set, str):
        index_set = FiniteEnumeratedSet(list(index_set))
    elif isinstance(index_set, (tuple, list)):
        index_set = FiniteEnumeratedSet(index_set)

    if ngens is None or ngens >= 0:
        if names is not None:
            if len(names) != index_set.cardinality():
                raise IndexError("the number of names must equal"
                                 " the size of the indexing set")
            if ngens is not None and len(names) != ngens:
                raise IndexError("the number of names must equal the"
                                 " number of generators")
        elif ngens is not None and index_set.cardinality() != ngens:
            raise IndexError("the size of the indexing set must equal"
                             " the number of generators")

    return (names, index_set)
