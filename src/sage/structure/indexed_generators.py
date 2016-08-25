# -*- coding: utf-8 -*-
"""
Indexed Generators
"""
#*****************************************************************************
#       Copyright (C) 2013      Travis Scrimshaw <tscrim at ucdavis.edu>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import Integer

class IndexedGenerators(object):
    r"""nodetex
    Abstract base class for parents whose elements consist of generators
    indexed by an arbitrary set.

    Options controlling the printing of elements:

    - ``prefix`` -- string, prefix used for printing elements of this
      module (optional, default 'x').  With the default, a monomial
      indexed by 'a' would be printed as ``x['a']``.

    - ``latex_prefix`` -- string or ``None``, prefix used in the `\LaTeX`
      representation of elements (optional, default ``None``). If this is
      anything except the empty string, it prints the index as a
      subscript.  If this is None, it uses the setting for ``prefix``,
      so if ``prefix`` is set to "B", then a monomial indexed by 'a'
      would be printed as ``B_{a}``.  If this is the empty string, then
      don't print monomials as subscripts: the monomial indexed by 'a'
      would be printed as ``a``, or as ``[a]`` if ``latex_bracket`` is
      True.

    - ``bracket`` -- ``None``, bool, string, or list or tuple of
      strings (optional, default ``None``): if ``None``, use the value of the
      attribute ``self._repr_option_bracket``, which has default value
      ``True``.  (``self._repr_option_bracket`` is available for backwards
      compatibility.  Users should set ``bracket`` instead.  If
      ``bracket`` is set to anything except ``None``, it overrides
      the value of ``self._repr_option_bracket``.)  If ``False``, do not
      include brackets when printing elements: a monomial indexed by
      'a' would be printed as ``B'a'``, and a monomial indexed by
      (1,2,3) would be printed as ``B(1,2,3)``.  If True, use "[" and
      "]" as brackets.  If it is one of "[", "(", or "{", use it and
      its partner as brackets.  If it is any other string, use it as
      both brackets.  If it is a list or tuple of strings, use the
      first entry as the left bracket and the second entry as the
      right bracket.

    - ``latex_bracket`` -- bool, string, or list or tuple of strings
      (optional, default False): if ``False``, do not include brackets in
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
      the print representation (optional, default "*")

    - ``latex_scalar_mult`` -- string or ``None`` (default: ``None``),
      string to use for scalar multiplication in the latex
      representation.  If None, use the empty string if ``scalar_mult``
      is set to "*", otherwise use the value of ``scalar_mult``.

    - ``tensor_symbol`` -- string or ``None`` (default: ``None``),
      string to use for tensor product in the print representation. If
      ``None``, use  ``sage.categories.tensor.symbol``.

    - ``generator_cmp`` -- deprecated

    - ``sorting_key`` -- a key function (default: ``lambda x: x``),
      to use for sorting elements in the output of elements

    - ``sorting_reverse`` -- bool (default: ``False``), if ``True`` 
      sort elements in reverse order in the output of elements

    - ``string_quotes`` -- bool (default: ``True``), if ``True`` then
      display string indices with quotes

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

        sage: I = IndexedGenerators(ZZ, prefix="", latex_bracket='(')
        sage: I._repr_generator(2)
        '[2]'
        sage: I._latex_generator(2)
        \left( 2 \right)

        sage: I = IndexedGenerators(ZZ, bracket=['|', '>'])
        sage: I._repr_generator(2)
        'x|2>'
    """
    def __init__(self, indices, prefix="x", **kwds):
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
                               'bracket': None,
                               'latex_bracket': False,
                               'latex_prefix': None,
                               'scalar_mult': "*",
                               'latex_scalar_mult': None,
                               'tensor_symbol': None,
                               'string_quotes': True,
                               'sorting_key': lambda x: x,
                               'sorting_reverse': False}
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

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F.indices()
            {'a', 'b', 'c'}
        """
        return self._indices

    def prefix(self):
        """
        Return the prefix used when displaying elements of self.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: F.prefix()
            'B'

        ::

            sage: X = SchubertPolynomialRing(QQ)
            sage: X.prefix()
            'X'
        """
        return self._print_options['prefix']

    def print_options(self, **kwds):
        """
        Return the current print options, or set an option.

        INPUT: all of the input is optional; if present, it should be
        in the form of keyword pairs, such as
        ``latex_bracket='('``.  The allowable keywords are:

        - ``prefix``
        - ``latex_prefix``
        - ``bracket``
        - ``latex_bracket``
        - ``scalar_mult``
        - ``latex_scalar_mult``
        - ``tensor_symbol``
        - ``string_quotes``
        - ``sorting_key``
        - ``sorting_reverse``

        See the documentation for :class:`IndexedGenerators` for
        descriptions of the effects of setting each of these options.

        OUTPUT: if the user provides any input, set the appropriate
        option(s) and return nothing.  Otherwise, return the
        dictionary of settings for print and LaTeX representations.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(ZZ, [1,2,3], prefix='x')
            sage: F.print_options()
            {...'prefix': 'x'...}
            sage: F.print_options(bracket='(')
            sage: F.print_options()
            {...'bracket': '('...}

        TESTS::

            sage: sorted(F.print_options().items())
            [('bracket', '('),
             ('latex_bracket', False), ('latex_prefix', None),
             ('latex_scalar_mult', None), ('prefix', 'x'),
             ('scalar_mult', '*'),
             ('sorting_key', <function <lambda> at ...>),
             ('sorting_reverse', False), ('string_quotes', True),
             ('tensor_symbol', None)]
            sage: F.print_options(bracket='[') # reset
            sage: F.print_options(generator_cmp=lambda x,y: (x < y) - (x > y))
            doctest:...: DeprecationWarning: Option generator_cmp is deprecated use sorting_key and sorting_reverse instead.
            See http://trac.sagemath.org/17229 for details.
        """
        # don't just use kwds.get(...) because I want to distinguish
        # between an argument like "option=None" and the option not
        # being there altogether.
        if kwds:
            if 'generator_cmp' in kwds:
                from sage.misc.superseded import deprecation
                deprecation(17229, "Option generator_cmp is deprecated use sorting_key and sorting_reverse instead.")
                from functools import cmp_to_key
                kwds['sorting_key'] = cmp_to_key(kwds['generator_cmp'])
                del kwds['generator_cmp']
            for option in kwds:
                # TODO: make this into a set and put it in a global variable?
                if option in ['prefix', 'latex_prefix', 'bracket', 'latex_bracket',
                              'scalar_mult', 'latex_scalar_mult', 'tensor_symbol',
                              'string_quotes',
                              'sorting_key', 'sorting_reverse'
                             ]:
                    self._print_options[option] = kwds[option]
                else:
                    raise ValueError('{} is not a valid print option.'.format(option))
            return
        return self._print_options

    _repr_option_bracket = True

    def _repr_generator(self, m):
        """
        Return a string representing the generator indexed by ``m``.

        The output can be customized by setting any of the following
        options when initializing the parent:

        - ``prefix``
        - ``bracket``
        - ``scalar_mult``

        Alternatively, one can use the :meth:`print_options` method
        to achieve the same effect.  To modify the bracket setting,
        one can also set ``self._repr_option_bracket`` as long as one
        has *not* set the ``bracket`` option: if the
        ``bracket`` option is anything but ``None``, it overrides
        the value of ``self._repr_option_bracket``.

        See the documentation for :class:`CombinatorialFreeModule` for
        details on the initialization options.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: e = F.basis()
            sage: e['a'] + 2*e['b']    # indirect doctest
            B['a'] + 2*B['b']

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'], prefix="F")
            sage: e = F.basis()
            sage: e['a'] + 2*e['b']    # indirect doctest
            F['a'] + 2*F['b']
            sage: F.print_options(string_quotes=False)
            sage: e['a'] + 2*e['b']
            F[a] + 2*F[b]

            sage: QS3 = CombinatorialFreeModule(QQ, Permutations(3), prefix="")
            sage: original_print_options = QS3.print_options()
            sage: a = 2*QS3([1,2,3])+4*QS3([3,2,1])
            sage: a                      # indirect doctest
            2*[[1, 2, 3]] + 4*[[3, 2, 1]]

            sage: QS3.print_options(bracket = False)
            sage: a              # indirect doctest
            2*[1, 2, 3] + 4*[3, 2, 1]

            sage: QS3.print_options(prefix='')
            sage: a              # indirect doctest
            2*[1, 2, 3] + 4*[3, 2, 1]

            sage: QS3.print_options(bracket="|", scalar_mult=" *@* ")
            sage: a              # indirect doctest
            2 *@* |[1, 2, 3]| + 4 *@* |[3, 2, 1]|

            sage: QS3.print_options(**original_print_options) # reset

        TESTS::

            sage: F = CombinatorialFreeModule(QQ, [('a', 'b'), ('c','d')])
            sage: e = F.basis()
            sage: e[('a','b')] + 2*e[('c','d')]    # indirect doctest
            B[('a', 'b')] + 2*B[('c', 'd')]
        """
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
        quotes = self._print_options.get('string_quotes', True)
        if not quotes and isinstance(m, str):
            return self.prefix() + left + m + right
        return self.prefix() + left + repr(m) + right # mind the (m), to accept a tuple for m

    def _ascii_art_generator(self, m):
        r"""
        Return an ascii art representing the generator indexed by ``m``.

        TESTS::

            sage: R = NonCommutativeSymmetricFunctions(QQ).R()
            sage: ascii_art(R[1,2,2,4])
            R
               ****
              **
             **
             *
            sage: Partitions.options(diagram_str="#", convention="french")
            sage: ascii_art(R[1,2,2,4])
            R
             #
             ##
              ##
               ####
            sage: Partitions.options._reset()
        """
        from sage.typeset.ascii_art import AsciiArt, ascii_art
        pref = AsciiArt([self.prefix()])
        r = pref * (AsciiArt([" "**Integer(len(pref))]) + ascii_art(m))
        r._baseline = r._h - 1
        return r

    def _unicode_art_generator(self, m):
        r"""
        Return an unicode art representing the generator indexed by ``m``.

        TESTS::

            sage: R = NonCommutativeSymmetricFunctions(QQ).R()
            sage: unicode_art(R[1,2,2,4])
            R
               ┌┬┬┬┐
              ┌┼┼┴┴┘
             ┌┼┼┘
             ├┼┘
             └┘
            sage: Partitions.options.convention="french"
            sage: unicode_art(R[1,2,2,4])
            R
             ┌┐
             ├┼┐
             └┼┼┐
              └┼┼┬┬┐
               └┴┴┴┘
            sage: Partitions.options._reset()
        """
        from sage.typeset.unicode_art import UnicodeArt, unicode_art
        pref = UnicodeArt([self.prefix()])
        r = pref * (UnicodeArt([" " ** Integer(len(pref))]) + unicode_art(m))
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

        (Alternatively, one can use the :meth:`print_options` method
        to achieve the same effect.)

        See the documentation for :class:`CombinatorialFreeModule` for
        details on the initialization options.

        EXAMPLES::

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'])
            sage: e = F.basis()
            sage: latex(e['a'] + 2*e['b'])    # indirect doctest
            B_{a} + 2B_{b}

            sage: F = CombinatorialFreeModule(QQ, ['a', 'b', 'c'], prefix="C")
            sage: e = F.basis()
            sage: latex(e['a'] + 2*e['b'])    # indirect doctest
            C_{a} + 2C_{b}

            sage: QS3 = CombinatorialFreeModule(QQ, Permutations(3), prefix="", scalar_mult="*")
            sage: original_print_options = QS3.print_options()
            sage: a = 2*QS3([1,2,3])+4*QS3([3,2,1])
            sage: latex(a)                     # indirect doctest
            2[1, 2, 3] + 4[3, 2, 1]
            sage: QS3.print_options(latex_bracket=True)
            sage: latex(a)                     # indirect doctest
            2\left[ [1, 2, 3] \right] + 4\left[ [3, 2, 1] \right]
            sage: QS3.print_options(latex_bracket="(")
            sage: latex(a)                     # indirect doctest
            2\left( [1, 2, 3] \right) + 4\left( [3, 2, 1] \right)
            sage: QS3.print_options(latex_bracket=('\\myleftbracket', '\\myrightbracket'))
            sage: latex(a)                     # indirect doctest
            2\myleftbracket [1, 2, 3] \myrightbracket + 4\myleftbracket [3, 2, 1] \myrightbracket
            sage: QS3.print_options(**original_print_options) # reset

        TESTS::

            sage: F = CombinatorialFreeModule(QQ, [('a', 'b'), (0,1,2)])
            sage: e = F.basis()
            sage: latex(e[('a','b')])    # indirect doctest
            B_{('a', 'b')}
            sage: latex(2*e[(0,1,2)])    # indirect doctest
            2B_{\left(0, 1, 2\right)}
            sage: F = CombinatorialFreeModule(QQ, [('a', 'b'), (0,1,2)], prefix="")
            sage: e = F.basis()
            sage: latex(2*e[(0,1,2)])    # indirect doctest
            2\left(0, 1, 2\right)
        """
        from sage.misc.latex import latex

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

