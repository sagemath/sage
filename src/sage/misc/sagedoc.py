# sage_setup: distribution = sagemath-repl
r"""
Format Sage documentation for viewing with IPython and the notebook

AUTHORS:

- William Stein (2005): initial version.
- Nick Alexander (2007): nodetex functions
- Nick Alexander (2008): search_src, search_def improvements
- Martin Albrecht (2008-03-21): parse LaTeX description environments in sagedoc
- John Palmieri (2009-04-11): fix for #5754 plus doctests
- Dan Drake (2009-05-21): refactor search_* functions, use system 'find' instead of sage -grep
- John Palmieri (2009-06-28): don't use 'find' -- use Python (os.walk, re.search) instead.
- Simon King (2011-09): Use os.linesep, avoid destruction of embedding information,
  enable nodetex in a docstring. Consequently use sage_getdoc.

TESTS:

Check that argspecs of extension function/methods appear correctly,
see :issue:`12849`::

    sage: from sage.env import SAGE_DOC
    sage: docfilename = os.path.join(SAGE_DOC, 'html', 'en', 'reference', 'calculus', 'sage', 'symbolic', 'expression.html')
    sage: with open(docfilename) as fobj:                                               # needs sagemath_doc_html
    ....:     for line in fobj:
    ....:         if "#sage.symbolic.expression.Expression.numerical_approx" in line:
    ....:             print(line)
    <span class="sig-name descname"><span class="pre">numerical_approx</span></span><span class="sig-paren">(</span><em class="sig-param"><span class="n"><span class="pre">prec</span></span><span class="o"><span class="pre">=</span></span><span class="default_value"><span class="pre">None</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">digits</span></span><span class="o"><span class="pre">=</span></span><span class="default_value"><span class="pre">None</span></span></em>, <em class="sig-param"><span class="n"><span class="pre">algorithm</span></span><span class="o"><span class="pre">=</span></span><span class="default_value"><span class="pre">None</span></span></em><span class="sig-paren">)</span>...

Check that sphinx is not imported at Sage start-up::

    sage: os.system("sage -c \"if 'sphinx' in sys.modules: sys.exit(1)\"")
    0
"""
# ****************************************************************************
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
import os
import re
import shutil
import sys
import pydoc
from sage.misc.temporary_file import tmp_dir
from sage.misc.viewer import browser
from sage.misc import sageinspect
import sage.version
from sage.env import SAGE_DOC, SAGE_SRC

# The detex function does two kinds of substitutions: math, which
# should only be done on the command line -- in the notebook, these
# should instead by taken care of by MathJax -- and nonmath, which
# should be done always.

# Math substitutions: don't forget the leading backslash '\\'. These
# are done using regular expressions, so it works best to also make
# the strings raw: r'\\blah'.
math_substitutes = [
    (r'\\to', '-->'),
    (r'\\rightarrow', '-->'),
    (r'\\leftarrow', '<--'),
    (r'\\leftrightarrow', '<->'),
    (r'\\longrightarrow', '--->'),
    (r'\\longleftarrow', '<---'),
    (r'\\longleftrightarrow', '<-->'),
    (r'\\Rightarrow', '==>'),
    (r'\\Leftarrow', '<=='),
    (r'\\Leftrightarrow', '<=>'),
    (r'\\Longrightarrow', '===>'),
    (r'\\Longleftarrow', '<==='),
    (r'\\Longleftrightarrow', '<==>'),
    (r'\\colon', ':'),
    (r'\\left', ''),
    (r'\\right', ''),
    (r'\\bigl', ''),
    (r'\\bigr', ''),
    (r'\\leq', '<='),
    (r'\\geq', '>='),
    (r'\\le', '<='),
    (r'\\ge', '>='),
    (r'\\cdots', '...'),
    (r'\\ldots', '...'),
    (r'\\dots', '...'),
    (r'\\cdot', ' *'),
    (r'\\ast', ' *'),
    (r' \\times', ' x'),
    (r'\\times', ' x'),
    (r'\\backslash', '\\'),
    (r'\\mapsto', ' |--> '),
    (r'\\longmapsto', ' |---> '),
    (r'\\lvert', '|'),
    (r'\\rvert', '|'),
    (r'\\mid', '|'),
    (r' \\circ', ' o'),
    (r'\\circ', ' o')
]
nonmath_substitutes = [
    ('\\_', '_'),
    ('\\item', '* '),
    ('<BLANKLINE>', ''),
    ('\\bf', ''),
    ('\\sage', 'Sage'),
    ('\\SAGE', 'Sage'),
    ('\\Sage', 'Sage'),
    ('\\rm', ''),
    ('backslash', '\\'),
    ('begin{enumerate}', ''),
    ('end{enumerate}', ''),
    ('begin{description}', ''),
    ('end{description}', ''),
    ('begin{itemize}', ''),
    ('end{itemize}', ''),
    ('begin{verbatim}', ''),
    ('end{verbatim}', ''),
    ('note{', 'NOTE: '),
]


def _rmcmd(s, cmd, left='', right=''):
    """
    Remove the LaTeX command ``cmd`` from the string ``s``.  This
    function is used by ``detex``.

    INPUT:

    - ``s`` -- string; string from which to remove the command

    - ``cmd`` -- string; command to be removed.  This should be a
      command which takes a single argument, like 'emph' or 'url'; the
      command is removed, but its argument is not.

    - ``left``, ``right`` -- string (default: ``''``); add these
      strings at the left and right ends of the command. See the
      examples.

    EXAMPLES::

        sage: from sage.misc.sagedoc import _rmcmd
        sage: _rmcmd('Check out \\url{http://www.sagemath.org}.', 'url')
        'Check out http://www.sagemath.org.'
        sage: _rmcmd('Text in \\emph{italics} looks like this.', 'emph', '*', '*')
        'Text in *italics* looks like this.'
        sage: _rmcmd('This is a \\very{silly} example.', 'very', right='!?')
        'This is a silly!? example.'
    """
    c = '\\%s{' % cmd
    while True:
        i = s.find(c)
        if i == -1:
            return s
        nesting = 1
        j = i + len(c) + 1
        while j < len(s) and nesting > 0:
            if s[j] == '{':
                nesting += 1
            elif s[j] == '}':
                nesting -= 1
            j += 1
        j -= 1  # j is position of closing '}'
        if j < len(s):
            s = s[:i] + left + s[i + len(c):j] + right + s[j + 1:]
        else:
            return s

# I wanted to be cool and use regexp's, but they aren't really
# useful, since really this is a parsing problem, because of
# nesting of commands, etc.   Since it doesn't have to be
# super super fast (it's a page of text scrolled to the user),
# the above works fine.

#
# import re
# def _rmcmd(s, cmd, left='', right=''):
#     c = '\\%s{.*}'%cmd
#     r = re.compile(c, re.DOTALL)
#     while True:
#         m = r.search(s)
#         if m is None: break
#         s = s[:m.start()] + left + s[m.start()+len(cmd)+1:m.end()-1] \
#             + right + s[m.end():]
#     return s


itempattern = re.compile(r"\\item\[?([^]]*)\]? *(.*)")
itemreplace = r"* \1 \2"


def detex(s, embedded=False):
    r"""nodetex
    This strips LaTeX commands from a string; it is used by the
    ``format`` function to process docstrings for display from the
    command line interface.

    INPUT:

    - ``s`` -- string
    - ``embedded`` -- boolean (default: ``False``)

    If ``embedded`` is ``False``, then do the replacements in both
    ``math_substitutes`` and ``nonmath_substitutes``.  If ``True``, then
    only do ``nonmath_substitutes``.

    OUTPUT: string

    EXAMPLES::

        sage: from sage.misc.sagedoc import detex
        sage: detex(r'Some math: `n \geq k`.  A website: \url{sagemath.org}.')
        'Some math: n >= k.  A website: sagemath.org.\n'
        sage: detex(r'More math: `x \mapsto y`.  {\bf Bold face}.')
        'More math: x  |-->  y.  { Bold face}.\n'
        sage: detex(r'`a, b, c, \ldots, z`')
        'a, b, c, ..., z\n'
        sage: detex(r'`a, b, c, \ldots, z`', embedded=True)
        '`a, b, c, \\ldots, z`'
        sage: detex(r'`\left(\lvert x\ast y \rvert\right]`')
        '(| x * y |]\n'
        sage: detex(r'`\left(\leq\le\leftarrow \rightarrow\unknownmacro\to`')
        '(<=<=<-- -->\\unknownmacro-->\n'
    """
    s = _rmcmd(s, 'url')
    s = _rmcmd(s, 'code')
    s = _rmcmd(s, 'class')
    s = _rmcmd(s, 'mbox')
    s = _rmcmd(s, 'text')
    s = _rmcmd(s, 'section')
    s = _rmcmd(s, 'subsection')
    s = _rmcmd(s, 'subsubsection')
    s = _rmcmd(s, 'note', 'NOTE: ', '')
    s = _rmcmd(s, 'emph', '*', '*')
    s = _rmcmd(s, 'textbf', '*', '*')

    s = re.sub(itempattern, itemreplace, s)

    for a, b in nonmath_substitutes:
        s = s.replace(a, b)
    if not embedded:  # not in the notebook
        s = _rmcmd(s, 'mathop')
        s = _rmcmd(s, 'mathrm')
        try:
            from .sphinxify import sphinxify
        except ImportError:
            s = s.replace('``', '"').replace('`', '') + '\n'
        else:
            s = sphinxify(s, format='text')
        # Do math substitutions. The strings to be replaced should be
        # TeX commands like "\\blah". Do a regular expression
        # replacement to replace "\\blah" but not "\\blahxyz", etc.:
        # test to make sure the next character is not a letter.
        for a, b in math_substitutes:
            s = re.sub(a + '([^a-zA-Z])', b + '\\1', s)
    return s


def skip_TESTS_block(docstring):
    r"""
    Remove blocks labeled "TESTS:" from ``docstring``.

    INPUT:

    - ``docstring`` -- string

    A "TESTS" block is a block starting "TESTS:" (or
    the same with two colons), on a line on its own, and ending either
    with a line indented less than "TESTS", or with a line with the
    same level of indentation -- not more -- matching one of the
    following:

    - a Sphinx directive of the form ".. foo:", optionally followed by
      other text.

    - text of the form "UPPERCASE:", optionally followed by other
      text.

    - lines which look like a reST header: one line containing
      anything, followed by a line consisting only of a string of
      hyphens, equal signs, or other characters which are valid
      markers for reST headers: ``- = ` : ' " ~ _ ^ * + # < >``.
      However, lines only containing double colons `::` do not
      end "TESTS" blocks.

    Return the string obtained from ``docstring`` by removing these
    blocks.

    EXAMPLES::

        sage: from sage.misc.sagedoc import skip_TESTS_block
        sage: start = ' Docstring\n\n'
        sage: test = ' TESTS: \n\n Here is a test::\n     sage: 2+2 \n     5 \n\n'
        sage: test2 = ' TESTS:: \n\n     sage: 2+2 \n     6 \n\n'

    Test lines starting with "REFERENCES:"::

        sage: refs = ' REFERENCES: \n text text \n'
        sage: skip_TESTS_block(start + test + refs).rstrip() == (start + refs).rstrip()
        True
        sage: skip_TESTS_block(start + test + test2 + refs).rstrip() == (start + refs).rstrip()
        True
        sage: skip_TESTS_block(start + test + refs + test2).rstrip() == (start + refs).rstrip()
        True

    Test Sphinx directives::

        sage: directive = ' .. todo:: \n     do some stuff \n'
        sage: skip_TESTS_block(start + test + refs + test2 + directive).rstrip() == (start + refs + directive).rstrip()
        True

    Test unindented lines::

        sage: unindented = 'NOT INDENTED\n'
        sage: skip_TESTS_block(start + test + unindented).rstrip() == (start + unindented).rstrip()
        True
        sage: skip_TESTS_block(start + test + unindented + test2 + unindented).rstrip() == (start + unindented + unindented).rstrip()
        True

    Test headers::

        sage: header = ' Header:\n ~~~~~~~~'
        sage: skip_TESTS_block(start + test + header) == start + header
        True

    Not a header because the characters on the second line must all be
    the same::

        sage: fake_header = ' Header:\n -=-=-=-=-='
        sage: skip_TESTS_block(start + test + fake_header).rstrip() == start.rstrip()
        True

    Not a header because it's indented compared to 'TEST' in the
    string ``test``::

        sage: another_fake = '\n    blah\n    ----'
        sage: skip_TESTS_block(start + test + another_fake).rstrip() == start.rstrip()
        True

    Double colons ``::`` are also not considered as headers (:issue:`27896`)::

        sage: colons = ' ::\n\n     sage: 2+2\n     4\n\n'
        sage: skip_TESTS_block(start + test2 + colons).rstrip() == start.rstrip()
        True
    """
    # tests_block: match a line starting with whitespace, then
    # "TEST" or "TESTS" followed by ":" or "::", then possibly
    # more whitespace, then the end of the line.
    tests_block = re.compile('([ ]*)TEST[S]?:[:]?[ ]*$')
    # end_of_block: match a line starting with whitespace, then Sphinx
    # directives of the form ".. foo:". This will match directive
    # names "foo" containing letters of either case, hyphens,
    # underscores.
    # Also match uppercase text followed by a colon, like
    # "REFERENCES:" or "ALGORITHM:".
    end_of_block = re.compile(r'[ ]*(\.\.[ ]+[-_A-Za-z]+|[A-Z]+):')
    # header: match a string of hyphens, or other characters which are
    # valid markers for reST headers: - = ` : ' " ~ _ ^ * + # < >
    # except for double colons ::
    header = re.compile(r'^[ ]*(?:([-=`\'"~_^*+#><])\1+|:|:::+)[ ]*$')
    s = ''
    skip = False
    previous = ''
    # indentation: amount of indentation at the start of 'TESTS:'.
    indentation = ''
    for l in docstring.split('\n'):
        if not skip:
            m = tests_block.match(l)
            if m:
                skip = True
                indentation = m.group(1)
            else:
                s += "\n"
                s += l
        else:
            if l and not l.isspace() and not l.startswith(indentation):
                # A non-blank line indented less than 'TESTS:'
                skip = False
                s += "\n"
                s += l
            elif end_of_block.match(l) and not tests_block.match(l):
                # A line matching end_of_block and indented the same as 'TESTS:'
                if l.startswith(indentation + " "):
                    continue
                skip = False
                s += "\n"
                s += l
            elif header.match(l):
                # A line matching header.
                if l.startswith(indentation + " "):
                    continue
                skip = False
                if previous:
                    s += "\n"
                    s += previous
                s += "\n"
                s += l
        previous = l
    return s[1:]  # Remove empty line from the beginning.


def process_dollars(s):
    r"""nodetex
    Replace dollar signs with backticks.

    More precisely, do a regular expression search.  Replace a plain
    dollar sign ($) by a backtick (`).  Replace an escaped dollar sign
    (\\$) by a dollar sign ($).  Don't change a dollar sign preceded or
    followed by a backtick (\`$ or \$`), because of strings like
    "``$HOME``".  Don't make any changes on lines starting with more
    spaces than the first nonempty line in ``s``, because those are
    indented and hence part of a block of code or examples.

    This also doesn't replaces dollar signs enclosed in curly braces,
    to avoid nested math environments.

    EXAMPLES::

        sage: from sage.misc.sagedoc import process_dollars
        sage: process_dollars('hello')
        'hello'
        sage: process_dollars('some math: $x=y$')
        doctest:warning...
        DeprecationWarning: using dollar signs to mark up math in Sage docstrings
        is deprecated; use backticks instead
        See https://github.com/sagemath/sage/issues/33973 for details.
        'some math: `x=y`'

    Replace \\$ with $, and don't do anything when backticks are involved::

        sage: process_dollars(r'a ``$REAL`` dollar sign: \$')
        'a ``$REAL`` dollar sign: $'

    Don't make any changes on lines indented more than the first
    nonempty line::

        sage: s = '\n first line\n     indented $x=y$'
        sage: s == process_dollars(s)
        True

    Don't replace dollar signs enclosed in curly braces::

        sage: process_dollars(r'f(n) = 0 \text{ if $n$ is prime}')
        'f(n) = 0 \\text{ if $n$ is prime}'

    This is not perfect::

        sage: process_dollars(r'$f(n) = 0 \text{ if $n$ is prime}$')
        '`f(n) = 0 \\text{ if $n$ is prime}$'

    The regular expression search doesn't find the last $.
    Fortunately, there don't seem to be any instances of this kind of
    expression in the Sage library, as of this writing.
    """
    if s.find("$") == -1:
        return s
    from sage.misc.superseded import deprecation
    # find how much leading whitespace s has, for later comparison:
    # ignore all $ on lines which start with more whitespace.
    whitespace = re.match(r'\s*\S', s.lstrip('\n'))
    whitespace = ' ' * (whitespace.end() - 1)  # leading whitespace
    # Indices will be a list of pairs of positions in s, to search between.
    # If the following search has no matches, then indices will be (0, len(s)).
    indices = [0]
    # This searches for "$blah$" inside a pair of curly braces --
    # don't change these, since they're probably coming from a nested
    # math environment.  So for each match, search to the left of its
    # start and to the right of its end, but not in between.
    for m in re.finditer(r"{[^{}$]*\$([^{}$]*)\$[^{}$]*}", s):
        indices[-1] = (indices[-1], m.start())
        indices.append(m.end())
    indices[-1] = (indices[-1], len(s))
    # regular expression for $ (not \$, `$, $`, and only on a line
    # with no extra leading whitespace).
    #
    # in detail:
    #   re.compile("^" # beginning of line
    #               + "(%s%)?" % whitespace
    #               + r"""(\S # non whitespace
    #                     .*?)? # non-greedy match any non-newline characters
    #                     (?<!`|\\)\$(?!`) # $ with negative lookbehind and lookahead
    #                  """, re.M | re.X)
    #
    # except that this doesn't work, so use the equivalent regular
    # expression without the 're.X' option.  Maybe 'whitespace' gets
    # eaten up by re.X?
    regexp = "^" + "(%s)?" % whitespace + r"(\S.*?)?(?<!`|\\)\$(?!`)"
    dollar = re.compile(regexp, re.M)
    # regular expression for \$
    slashdollar = re.compile(r"\\\$")
    for start, end in indices:
        while dollar.search(s, start, end):
            m = dollar.search(s, start, end)
            s = s[:m.end() - 1] + "`" + s[m.end():]
            deprecation(33973,
                        "using dollar signs to mark up math in Sage docstrings "
                        "is deprecated; use backticks instead")
        while slashdollar.search(s, start, end):
            m = slashdollar.search(s, start, end)
            s = s[:m.start()] + "$" + s[m.end():]
    return s


# When adding roles here, also add them to SAGE_ROOT/src/tox.ini [flake8]
# and document them in SAGE_ROOT/src/doc/en/developer/sage_manuals.rst
#
# Sage github issue shortcuts. For example, :issue:`7549`.
pythonversion = sys.version.split(' ')[0]
extlinks = {
    'python': (f'https://docs.python.org/release/{pythonversion}/%s', None),
    'issue': ('https://github.com/sagemath/sage/issues/%s', 'Issue #%s'),
    'sage_root': ('https://github.com/sagemath/sage/tree/develop/%s', 'SAGE_ROOT/%s'),
    'wikipedia': ('https://en.wikipedia.org/wiki/%s', 'Wikipedia article %s'),
    'arxiv': ('https://arxiv.org/abs/%s', 'arXiv %s'),
    'oeis': ('https://oeis.org/%s', 'OEIS sequence %s'),
    'doi': ('https://doi.org/%s', 'doi:%s'),
    'pari': ('https://pari.math.u-bordeaux.fr/dochtml/help/%s', 'pari:%s'),
    'mathscinet': ('https://www.ams.org/mathscinet-getitem?mr=%s', 'MathSciNet %s'),
    'common_lisp': ('https://www.lispworks.com/documentation/lw50/CLHS/Body/%s.htm', 'Common Lisp: %s'),
    'ecl': ('https://ecl.common-lisp.dev/static/manual/%s.html', 'ECL: %s'),
    'gap': ('https://docs.gap-system.org/doc/ref/%s_mj.html', 'GAP: %s'),
    'gap_package': ('https://docs.gap-system.org/pkg/%s', 'GAP package %s'),
    'giac_cascmd': ('https://www-fourier.ujf-grenoble.fr/~parisse/giac/doc/en/cascmd_en/%s.html', 'Giac: %s'),
    'giac_us': ('https://www-fourier.ujf-grenoble.fr/~parisse/giac_us.html#%s', 'Giac API: %s'),
    'maxima': ('https://maxima.sourceforge.io/docs/manual/maxima_singlepage.html#%s', 'Maxima: %s'),
    'meson': ('https://mesonbuild.com/%s', 'Meson: %s'),
    'polymake': ('https://polymake.org/doku.php/documentation/latest/%s', 'polymake: %s'),
    'ppl': ('https://www.bugseng.com/products/ppl/documentation/user/ppl-user-1.2-html/%s.html', 'PPL: %s'),
    'qepcad': ('https://www.usna.edu/CS/qepcadweb/B/%s.html', 'QEPCAD: %s'),
    'scip': ('https://scipopt.org/doc/html/%s.php', 'SCIP: %s'),
    'singular': ('https://www.singular.uni-kl.de/Manual/4-3-2/%s.htm', 'Singular: %s'),
    'soplex': ('https://soplex.zib.de/doc/html/%s.php', 'SoPlex: %s'),
}


def process_extlinks(s, embedded=False):
    r"""nodetex

    In docstrings at the command line, process markup related to the
    Sphinx extlinks extension. For example, replace ``:issue:`NUM```
    with ``https://github.com/sagemath/sage/issues/NUM``, and similarly with
    ``:python:TEXT`` and ``:wikipedia:TEXT``, looking up the url from
    the dictionary ``extlinks`` in ``sage_docbuild.conf``.
    If ``TEXT`` is of the form ``blah <LINK>``, then it uses ``LINK``
    rather than ``TEXT`` to construct the url.

    In the notebook, don't do anything: let sphinxify take care of it.

    INPUT:

    - ``s`` -- string, in practice a docstring
    - ``embedded`` -- boolean (default: ``False``)

    This function is called by :func:`format`, and if in the notebook,
    it sets ``embedded`` to be ``True``, otherwise ``False``.

    EXAMPLES::

        sage: from sage.misc.sagedoc import process_extlinks
        sage: process_extlinks('See :issue:`1234`, :wikipedia:`Wikipedia <Sage_(mathematics_software)>`, and :issue:`4321` ...')
        'See https://github.com/sagemath/sage/issues/1234, https://en.wikipedia.org/wiki/Sage_(mathematics_software), and https://github.com/sagemath/sage/issues/4321 ...'
        sage: process_extlinks('See :issue:`1234` for more information.', embedded=True)
        'See :issue:`1234` for more information.'
        sage: process_extlinks('see :python:`Implementing Descriptors <reference/datamodel.html#implementing-descriptors>` ...')
        'see https://docs.python.org/release/.../reference/datamodel.html#implementing-descriptors ...'
    """
    if embedded:
        return s
    for key in extlinks:
        while True:
            m = re.search(':%s:`([^`]*)`' % key, s)
            if not m:
                break
            link = m.group(1)
            m = re.search('.*<([^>]*)>', link)
            if m:
                link = m.group(1)
            s = re.sub(':%s:`([^`]*)`' % key,
                       extlinks[key][0].replace('%s', link),
                       s, count=1)
    return s


def process_mathtt(s):
    r"""nodetex
    Replace \\mathtt{BLAH} with BLAH in the command line.

    INPUT:

    - ``s`` -- string, in practice a docstring

    This function is called by :func:`format`.

    EXAMPLES::

        sage: from sage.misc.sagedoc import process_mathtt
        sage: process_mathtt(r'e^\mathtt{self}')
        'e^self'
    """
    while True:
        start = s.find("\\mathtt{")
        end = s.find("}", start)
        if start == -1 or end == -1:
            break
        s = s[:start] + s[start + 8:end] + s[end + 1:]
    return s


def process_optional_doctest_tags(s):
    r"""
    Remove ``# optional/needs`` doctest tags for present features from docstring ``s``.

    EXAMPLES:

        sage: from sage.misc.sagedoc import process_optional_doctest_tags
        sage: process_optional_doctest_tags("sage: # needs sage.rings.finite_rings\nsage: K.<x> = FunctionField(GF(5^2,'a')); K\nRational function field in x over Finite Field in a of size 5^2")  # needs sage.rings.finite_rings
        "sage: K.<x> = FunctionField(GF(5^2,'a')); K\nRational function field in x over Finite Field in a of size 5^2"
    """
    import io
    from sage.doctest.external import available_software
    from sage.doctest.parsing import parse_optional_tags, update_optional_tags

    start = 0
    with io.StringIO() as output:
        for m in re.finditer('( *sage: *.*#.*)\n', s):
            output.write(s[start:m.start(0)])
            line = m.group(1)
            tags = [tag for tag in parse_optional_tags(line)
                    if tag not in available_software]
            line = update_optional_tags(line, tags=tags)
            if not re.fullmatch(' *sage: *', line):
                print(line, file=output)
            start = m.end(0)
        output.write(s[start:])
        return output.getvalue()


def format(s, embedded=False):
    r"""noreplace
    Format Sage documentation ``s`` for viewing with IPython.

    This calls :func:`detex` on ``s`` to convert LaTeX commands to plain
    text, unless the directive ``nodetex`` is given in the first line
    of the string.

    Also, if ``s`` contains a string of the form ``<<<obj>>>``, then
    it replaces it with the docstring for ``obj``, unless the
    directive ``noreplace`` is given in the first line. If an error
    occurs under the attempt to find the docstring for ``obj``, then
    the substring ``<<<obj>>>`` is preserved.

    Directives must be separated by a comma.

    INPUT:

    - ``s`` -- string
    - ``embedded`` -- boolean (default: ``False``)

    OUTPUT: string

    Set ``embedded`` equal to ``True`` if formatting for use in the
    notebook; this just gets passed as an argument to :func:`detex`.

    .. SEEALSO::

        :func:`sage.misc.sageinspect.sage_getdoc` to get the formatted
        documentation of a given object.

    EXAMPLES::

        sage: from sage.misc.sagedoc import format
        sage: format('Let `A` be an `m` by `n` (0,1)-matrix. We identify `A` with a chessboard')
        'Let A be an m by n (0,1)-matrix. We identify A with a chessboard\n'

    If the first line of the string is 'nodetex', remove 'nodetex' but
    don't modify any TeX commands::

        sage: format("nodetex\n`x \\geq y`")
        '`x \\geq y`'

    Testing a string enclosed in triple angle brackets::

        sage: format('<<<identity_matrix')
        '<<<identity_matrix\n'
        sage: format('identity_matrix>>>')
        'identity_matrix>>>\n'
        sage: format('<<<identity_matrix>>>')                                           # needs sage.modules
        '...Definition: identity_matrix(...'
        sage: format('<<<identity_matrix>>>')[:28]                                      # needs sphinx
        'Definition: identity_matrix('

    TESTS:

    We check that the todo Sphinx extension is correctly activated::

        sage: sage.misc.sagedoc.format(sage.combinat.ranker.on_fly.__doc__)             # needs sphinx
        "...Return ...Todo: add tests as in combinat::rankers\n"

    In the following use case, the ``nodetex`` directive would have been ignored prior
    to :issue:`11815`::

        sage: cython_code = ["def testfunc(x):",
        ....: "    '''",
        ....: "    nodetex",
        ....: "    This is a doc string with raw latex",
        ....: "",
        ....: "    `x \\geq y`",
        ....: "    '''",
        ....: "    return -x"]
        sage: cython('\n'.join(cython_code))                                            # needs sage.misc.cython
        sage: from sage.misc.sageinspect import sage_getdoc
        sage: print(sage_getdoc(testfunc))                                              # needs sage.misc.cython
        <BLANKLINE>
            This is a doc string with raw latex
        <BLANKLINE>
            `x \geq y`
        <BLANKLINE>

    We check that the ``noreplace`` directive works, even combined with
    ``nodetex`` (see :issue:`11817`)::

        sage: print(format('''nodetex, noreplace\n<<<identity_matrix>>>`\\not= 0`'''))
        <<<identity_matrix>>>`\not= 0`

    If replacement is impossible, then no error is raised::

        sage: print(format('<<<bla\n<<<bla>>>\n<<<identity_matrix>>>'))
        <<<bla <<<bla>>>
        <BLANKLINE>
        Definition: identity_matrix(ring, n=0, sparse=False)
        <BLANKLINE>
        This function is available as identity_matrix(...) and
        matrix.identity(...).
        <BLANKLINE>
           Return the n x n identity matrix over the given ring.
        ...

    Check that backslashes are preserved in code blocks (:issue:`29140`)::

        sage: format('::\n'                                                             # needs sphinx
        ....:        '\n'
        ....:        r'    sage: print(r"\\\\.")' '\n'
        ....:        r'    \\\\.')
        '   sage: print(r"\\\\\\\\.")\n   \\\\\\\\.\n'
        sage: format(r'inline code ``\\\\.``')
        'inline code "\\\\\\\\."\n'
    """
    if not isinstance(s, str):
        raise TypeError("s must be a string")

    # Leading empty lines must be removed, since we search for directives
    # in the first line.
    s = s.lstrip(os.linesep)

    # parse directives at beginning of docstring
    # currently, only 'nodetex' and 'noreplace' are supported.
    # 'no' + 'doctest' may be supported eventually (don't type that as
    # one word, or the whole file will not be doctested).
    first_newline = s.find(os.linesep)
    if first_newline > -1:
        first_line = s[:first_newline]
    else:
        first_line = s
    # Moreover, we must strip blank space in order to get the directives
    directives = [d.strip().lower() for d in first_line.split(',')]

    if 'noreplace' in directives or 'nodetex' in directives:
        s = s[first_newline + len(os.linesep):]

    try:
        import sage.all
    except ImportError:
        pass

    docs = set()
    if 'noreplace' not in directives:
        i_0 = 0
        while True:
            i = s[i_0:].find("<<<")
            if i == -1:
                break
            j = s[i_0 + i + 3:].find('>>>')
            if j == -1:
                break
            obj = s[i_0 + i + 3:i_0 + i + 3 + j]
            if obj in docs:
                t = ''
            else:
                try:
                    x = eval('sage.all.%s' % obj, locals())
                except AttributeError:
                    # A pair <<<...>>> has been found, but the object not.
                    i_0 += i + 6 + j
                    continue
                except SyntaxError:
                    # This is a simple heuristics to cover the case of
                    # a non-matching set of <<< and >>>
                    i_0 += i + 3
                    continue
                t0 = sage.misc.sageinspect.sage_getdef(x, obj)
                t1 = sage.misc.sageinspect.sage_getdoc(x)
                t = 'Definition: ' + t0 + '\n\n' + t1
                docs.add(obj)
            s = s[:i_0 + i] + '\n' + t + s[i_0 + i + 6 + j:]
            i_0 += i

    if 'nodetex' not in directives:
        s = process_dollars(s)
        s = skip_TESTS_block(s)
        if not embedded:
            s = process_mathtt(s)
        s = process_extlinks(s, embedded=embedded)
        s = detex(s, embedded=embedded)

    if not embedded:
        s = process_optional_doctest_tags(s)

    return s


def format_src(s):
    """
    Format Sage source code ``s`` for viewing with IPython.

    If ``s`` contains a string of the form "<<<obj>>>", then it
    replaces it with the source code for "obj".

    INPUT:

    - ``s`` -- string

    OUTPUT: string

    EXAMPLES::

        sage: from sage.misc.sagedoc import format_src
        sage: format_src('unladen swallow')
        'unladen swallow'
        sage: format_src('<<<Sq>>>')[5:15]                                              # needs sage.combinat sage.modules
        'Sq(*nums):'
    """
    if not isinstance(s, str):
        raise TypeError("s must be a string")
    docs = set()

    try:
        import sage.all
    except ImportError:
        pass

    while True:
        i = s.find("<<<")
        if i == -1:
            break
        j = s[i + 3:].find('>>>')
        if j == -1:
            break
        obj = s[i + 3:i + 3 + j]
        if obj in docs:
            t = ''
        else:
            x = eval('sage.all.%s' % obj, locals())
            t = my_getsource(x)
            docs.add(obj)
        if t is None:
            print(x)
            t = ''
        s = s[:i] + '\n' + t + s[i + 6 + j:]

    return s


###############################

def _search_src_or_doc(what, string, extra1='', extra2='', extra3='',
                       extra4='', extra5='', **kwargs):
    r"""
    Search the Sage library or documentation for lines containing
    ``string`` and possibly some other terms. This function is used by
    :func:`search_src`, :func:`search_doc`, and :func:`search_def`.

    INPUT:

    - ``what`` -- either ``'src'`` or ``'doc'``, according to whether you
      are searching the documentation or source code
    - the rest of the input is the same as :func:`search_src`,
      :func:`search_doc`, and :func:`search_def`

    OUTPUT:

    If ``interact`` is ``False``, a string containing the results;
    otherwise, there is no output and the results are presented
    according to whether you are using the notebook or command-line
    interface. In the command-line interface, each line of the results
    has the form ``filename:num:line of code``, where ``num`` is the
    line number in ``filename`` and ``line of code`` is the line that
    matched your search terms.

    EXAMPLES::

        sage: from sage.misc.sagedoc import _search_src_or_doc
        sage: print(_search_src_or_doc('src', r'matrix\(',  # long time random
        ....:       'incidence_structures', 'self',
        ....:       '^combinat', interact=False))
        misc/sagedoc.py:        sage: _search_src_or_doc('src', 'matrix(', 'incidence_structures', 'self', '^combinat', interact=False)
        combinat/designs/incidence_structures.py:        M1 = self.incidence_matrix()
        combinat/designs/incidence_structures.py:        A = self.incidence_matrix()
        combinat/designs/incidence_structures.py:        M = transpose(self.incidence_matrix())
        combinat/designs/incidence_structures.py:    def incidence_matrix(self):
        combinat/designs/incidence_structures.py:        A = self.incidence_matrix()
        combinat/designs/incidence_structures.py:        A = self.incidence_matrix()
        combinat/designs/incidence_structures.py:        #A = self.incidence_matrix()

    TESTS:

    The examples are nice, but marking them "random" means we're not
    really testing if the function works, just that it completes. These
    tests aren't perfect, but are reasonable.

    ::

        sage: from sage.misc.sagedoc import _search_src_or_doc
        sage: len(_search_src_or_doc('src', r'matrix\(', 'incidence_structures', 'self', 'combinat', interact=False).splitlines()) > 1
        True
        sage: 'abvar/homology' in _search_src_or_doc('doc', 'homology', 'variety', interact=False)  # long time (4s on sage.math, 2012), needs sagemath_doc_html
        True
        sage: 'divisors' in _search_src_or_doc('src', '^ *def prime', interact=False)
        True

    When passing ``interactive=True``, in a terminal session this will pass the
    ``text/plain`` output to the configured pager, while in a notebook session
    it will display the ``text/html`` output in the notebook's pager.  However,
    in a non-interactive session (as in the doctests) it should just print the
    results to stdout::

        sage: from sage.misc.sagedoc import _search_src_or_doc
        sage: _search_src_or_doc('src',                        # long time
        ....:                    r'def _search_src_or_doc\(',
        ....:                    interact=True)
        misc/sagedoc.py:...:        def _search_src_or_doc(what, string, extra1='', extra2='', extra3='',
    """

    # process keyword arguments
    interact = kwargs.get('interact', True)
    path_re = kwargs.get('path_re', '')
    module = kwargs.get('module', 'sage')
    whole_word = kwargs.get('whole_word', False)
    ignore_case = kwargs.get('ignore_case', True)
    multiline = kwargs.get('multiline', False)

    # done processing keywords
    # define module, exts (file extension), title (title of search),
    # base_path (top directory in which to search)
    if what == 'src':
        base_path = SAGE_SRC
        if module.find('sage') == 0:
            module = module[4:].lstrip(".")  # remove 'sage' or 'sage.' from module
            base_path = os.path.join(base_path, 'sage')
        module = module.replace(".", os.sep)
        exts = ['py', 'pyx', 'pxd']
        title = 'Source Code'
    else:
        module = ''
        exts = ['html']
        title = 'Documentation'
        base_path = os.path.join(SAGE_DOC, 'html')
    if not os.path.exists(base_path):
        print("""Warning: the Sage documentation is not available""")

    strip = len(base_path)
    results = []
    # in regular expressions, '\bWORD\b' matches 'WORD' but not
    # 'SWORD' or 'WORDS'.  so if the user requests a whole_word
    # search, append and prepend '\b' to each string.
    regexp = string
    extra_regexps = extras = [extra1, extra2, extra3, extra4, extra5]
    if whole_word:
        regexp = r'\b' + regexp + r'\b'
        extra_regexps = [r'\b%s\b' % e for e in extra_regexps]
    if ignore_case:
        # 'flags' is a flag passed to re.search. use bit-wise or "|" to combine flags.
        flags = re.IGNORECASE
    else:
        flags = 0

    # done with preparation; ready to start search
    for dirpath, dirs, files in os.walk(os.path.join(base_path, module)):
        try:
            dirs.remove('_static')
        except ValueError:
            pass
        for f in files:
            if not f.startswith('.') and re.search(r"\.(" + "|".join(exts) + ")$", f):
                filename = os.path.join(dirpath, f)
                if re.search(path_re, filename):
                    if multiline:
                        with open(filename) as fobj:
                            line = fobj.read()
                        if re.search(regexp, line, flags):
                            match_list = line
                        else:
                            match_list = None
                        for extra in extra_regexps:
                            if extra and match_list:
                                if not re.search(extra, match_list):
                                    match_list = None
                        if match_list:
                            results.append(filename[strip:].lstrip("/") + '\n')
                    else:
                        with open(filename) as fobj:
                            match_list = [(lineno, line)
                                          for lineno, line in enumerate(fobj)
                                          if re.search(regexp, line, flags)]
                        for extra in extra_regexps:
                            if extra:
                                match_list = [s for s in match_list
                                              if re.search(extra, s[1], re.MULTILINE | flags)]
                        for num, line in match_list:
                            results.append('{}:{}:{}'.format(
                                filename[strip:].lstrip('/'), num + 1, line))

    text_results = ''.join(results).rstrip()

    if not interact:
        return text_results

    html_results = format_search_as_html(title, results, [string] + extras)
    # potentially used below

    # Pass through the IPython pager in a mime bundle
    from IPython.core.page import page
    if not isinstance(text_results, str):
        text_results = text_results.decode('utf-8', 'replace')

    page({
        'text/plain': text_results,
        # 'text/html': html_results
        # don't return HTML results since they currently are not
        # correctly formatted for Jupyter use
    })


def search_src(string, extra1='', extra2='', extra3='', extra4='',
               extra5='', **kwds):
    r"""
    Search Sage library source code for lines containing ``string``.
    The search is case-insensitive by default.

    INPUT:

    - ``string`` -- string to find in the Sage source code

    - ``extra1``, ..., ``extra5`` -- additional strings to require when
      searching.  Lines must match all of these, as well as ``string``

    - ``whole_word`` -- (default: ``False``) if ``True``, search for
      ``string`` and ``extra1`` (etc.) as whole words only.  This
      assumes that each of these arguments is a single word, not a
      regular expression, and it might have unexpected results if used
      with regular expressions.

    - ``ignore_case`` -- boolean (default: ``True``); if ``False``, perform a
      case-sensitive search

    - ``multiline`` -- (default: ``False``) if ``True``, search more
      than one line at a time.  In this case, print any matching file
      names, but don't print line numbers.

    - ``interact`` -- boolean (default: ``True``); if ``False``, return
      a string with all the matches. Otherwise, this function returns
      ``None``, and the results are displayed appropriately, according
      to whether you are using the notebook or the command-line
      interface. You should not ordinarily need to use this.

    - ``path_re`` -- (default: ``''``) regular expression which
      the filename (including the path) must match

    - ``module`` -- (default: ``'sage'``) the module in which to
      search.  The default is 'sage', the entire Sage library.  If
      ``module`` doesn't start with "sage", then the links in the
      notebook output may not function.

    OUTPUT:

    If ``interact`` is ``False``, then return a string with all of
    the matches, separated by newlines.  On the other hand, if
    ``interact`` is ``True`` (the default), there is no output.  Instead:
    at the command line, the search results are printed on the screen
    in the form ``filename:line_number:line of text``, showing the
    filename in which each match occurs, the line number where it
    occurs, and the actual matching line.  (If ``multiline`` is ``True``,
    then only the filename is printed for each match.)  The file paths
    in the output are relative to ``$SAGE_SRC``.  In the
    notebook, each match produces a link to the actual file in which
    it occurs.

    The ``string`` and ``extraN`` arguments are treated as regular
    expressions, as is ``path_re``, and errors will be raised if they
    are invalid. The matches will be case-insensitive unless
    ``ignore_case`` is ``False``.

    .. NOTE::

        The ``extraN`` parameters are present only because
        ``search_src(string, *extras, interact=False)``
        is not parsed correctly by Python 2.6; see http://bugs.python.org/issue1909.

    EXAMPLES:

    First note that without using ``interact=False``, this function
    produces no output, while with ``interact=False``, the output is a
    string.  These examples almost all use this option, so that they
    have something to which to compare their output.

    You can search for "matrix" by typing ``search_src("matrix")``.
    This particular search will produce many results::

        sage: len(search_src("matrix", interact=False).splitlines())  # random # long time
        9522

    You can restrict to the Sage calculus code with
    ``search_src("matrix", module="sage.calculus")``, and this
    produces many fewer results::

        sage: len(search_src("matrix", module="sage.calculus", interact=False).splitlines())  # random
        26

    Note that you can do tab completion on the ``module`` string.
    Another way to accomplish a similar search::

        sage: len(search_src("matrix", path_re="calc",                                  # needs sage.modules
        ....:                interact=False).splitlines()) > 15
        True

    The following produces an error because the string 'fetch(' is a
    malformed regular expression::

        sage: try:
        ....:     print(search_src(" fetch(", "def", interact=False))
        ....: except Exception as e:
        ....:     print(e)
        missing ), unterminated subpattern at position 6

    To fix this, *escape* the parenthesis with a backslash::

        sage: print(search_src(r" fetch\(", "def", interact=False))  # random # long time
        matrix/matrix0.pyx:    cdef fetch(self, key):
        matrix/matrix0.pxd:    cdef fetch(self, key)

        sage: print(search_src(r" fetch\(", "def", "pyx", interact=False))  # random # long time
        matrix/matrix0.pyx:    cdef fetch(self, key):

    As noted above, the search is case-insensitive, but you can make it
    case-sensitive with the 'ignore_case' key word::

        sage: s = search_src('Matrix', path_re='matrix', interact=False); s.find('x') > 0
        True

        sage: s = search_src('MatRiX', path_re='matrix', interact=False); s.find('x') > 0
        True

        sage: s = search_src('MatRiX', path_re='matrix',
        ....:                interact=False, ignore_case=False); s.find('x') > 0
        False

    Searches are by default restricted to single lines, but this can
    be changed by setting ``multiline`` to be True.  In the following,
    since ``search_src(string, interact=False)`` returns a string with
    one line for each match, counting the length of
    ``search_src(string, interact=False).splitlines()`` gives the
    number of matches. ::

        sage: len(search_src('log', 'derivative', interact=False).splitlines()) < 40
        True
        sage: len(search_src('log', 'derivative',
        ....:                interact=False, multiline=True).splitlines()) > 70
        True

    A little recursive narcissism: let's do a doctest that searches for
    this function's doctests. Note that you can't put "sage:" in the
    doctest string because it will get replaced by the Python ">>>"
    prompt.

    ::

        sage: print(search_src(r'^ *sage[:] .*search_src\(', interact=False))  # long time
        misc/sagedoc.py:... len(search_src("matrix", interact=False).splitlines())...
        misc/sagedoc.py:... len(search_src("matrix", module="sage.calculus", interact=False).splitlines())...
        misc/sagedoc.py:... len(search_src("matrix", path_re="calc"...
        misc/sagedoc.py:... print(search_src(" fetch(", "def", interact=False))...
        misc/sagedoc.py:... print(search_src(r" fetch\(", "def", interact=False))...
        misc/sagedoc.py:... print(search_src(r" fetch\(", "def", "pyx", interact=False))...
        misc/sagedoc.py:... s = search_src('Matrix', path_re='matrix', interact=False); s.find('x') > 0...
        misc/sagedoc.py:... s = search_src('MatRiX', path_re='matrix', interact=False); s.find('x') > 0...
        misc/sagedoc.py:... s = search_src('MatRiX', path_re='matrix',...
        misc/sagedoc.py:... len(search_src('log', 'derivative', interact=False).splitlines()) < 40...
        misc/sagedoc.py:... len(search_src('log', 'derivative'...
        misc/sagedoc.py:... print(search_src(r'^ *sage[:] .*search_src\(', interact=False))...
        misc/sagedoc.py:... len(search_src("matrix", interact=False).splitlines()) > 9000...
        misc/sagedoc.py:... print(search_src('matrix', 'column', 'row', 'sub',...
        misc/sagedoc.py:... sage: results = search_src('format_search_as_html',...

    TESTS:

    As of this writing, there are about 9500 lines in the Sage library that
    contain "matrix"; it seems safe to assume we will continue to have
    over 9000 such lines::

        sage: len(search_src("matrix", interact=False).splitlines()) > 9000 # long time
        True

    Check that you can pass 5 parameters::

        sage: print(search_src('matrix', 'column', 'row', 'sub', 'start', 'index', interact=False)) # random # long time
        matrix/matrix0.pyx:598:        Get The 2 x 2 submatrix of M, starting at row index and column
        matrix/matrix0.pyx:607:        Get the 2 x 3 submatrix of M starting at row index and column index
        matrix/matrix0.pyx:924:        Set the 2 x 2 submatrix of M, starting at row index and column
        matrix/matrix0.pyx:933:        Set the 2 x 3 submatrix of M starting at row index and column
    """
    return _search_src_or_doc('src', string, extra1=extra1, extra2=extra2,
                              extra3=extra3, extra4=extra4, extra5=extra5,
                              **kwds)


def search_doc(string, extra1='', extra2='', extra3='', extra4='',
               extra5='', **kwds):
    r"""
    Search Sage HTML documentation for lines containing ``string``. The
    search is case-insensitive by default.

    The file paths in the output are relative to ``$SAGE_DOC``.

    INPUT: same as for :func:`search_src`.

    OUTPUT: same as for :func:`search_src`

    EXAMPLES:

    See the documentation for :func:`search_src` for more examples. ::

        sage: search_doc('creates a polynomial', path_re='tutorial', interact=False) # random
        html/en/tutorial/tour_polynomial.html:<p>This creates a polynomial ring and tells Sage to use (the string)

    If you search the documentation for 'tree', then you will get too
    many results, because many lines in the documentation contain the
    word 'toctree'.  If you use the ``whole_word`` option, though, you
    can search for 'tree' without returning all of the instances of
    'toctree'.  In the following, since ``search_doc('tree',
    interact=False)`` returns a string with one line for each match,
    counting the length of ``search_doc('tree',
    interact=False).splitlines()`` gives the number of matches. ::

        sage: # long time, needs sagemath_doc_html
        sage: N = len(search_doc('tree', interact=False).splitlines())
        sage: L = search_doc('tree', whole_word=True, interact=False).splitlines()
        sage: len(L) < N
        True
        sage: import re
        sage: tree_re = re.compile(r'(^|\W)tree(\W|$)', re.I)
        sage: all(tree_re.search(l) for l in L)
        True
    """
    return _search_src_or_doc('doc', string, extra1=extra1, extra2=extra2,
                              extra3=extra3, extra4=extra4, extra5=extra5,
                              **kwds)


def search_def(name, extra1='', extra2='', extra3='', extra4='',
               extra5='', **kwds):
    r"""
    Search Sage library source code for function definitions containing
    ``name``. The search is case-insensitive by default.

    INPUT: same as for :func:`search_src`.

    OUTPUT: same as for :func:`search_src`

    .. NOTE::

        The regular expression used by this function only finds function
        definitions that are preceded by spaces, so if you use tabs on a
        "def" line, this function will not find it. As tabs are not
        allowed in Sage library code, this should not be a problem.

    EXAMPLES:

    See the documentation for :func:`search_src` for more examples. ::

        sage: print(search_def("fetch", interact=False))  # random # long time
        matrix/matrix0.pyx:    cdef fetch(self, key):
        matrix/matrix0.pxd:    cdef fetch(self, key)

        sage: print(search_def("fetch", path_re="pyx", interact=False))  # random # long time
        matrix/matrix0.pyx:    cdef fetch(self, key):
    """
    # since we convert name to a regular expression, we need to do the
    # 'whole_word' conversion here, rather than pass it on to
    # _search_src_or_doc.
    if 'whole_word' in kwds and kwds['whole_word']:
        name = r'\b' + name + r'\b'
        if extra1:
            extra1 = r'\b' + extra1 + r'\b'
        if extra2:
            extra2 = r'\b' + extra2 + r'\b'
        if extra3:
            extra3 = r'\b' + extra3 + r'\b'
        if extra4:
            extra4 = r'\b' + extra4 + r'\b'
        if extra5:
            extra5 = r'\b' + extra5 + r'\b'
        kwds['whole_word'] = False

    return _search_src_or_doc('src', '^ *[c]?def.*%s' % name, extra1=extra1,
                              extra2=extra2, extra3=extra3, extra4=extra4,
                              extra5=extra5, **kwds)


def format_search_as_html(what, results, search):
    r"""
    Format the output from ``search_src``, ``search_def``, or
    ``search_doc`` as html, for use in the notebook.

    INPUT:

    - ``what`` -- string; what was searched (source code or
      documentation)
    - ``results`` -- string or list; the results of the search as a string or list of
      search results
    - ``search`` -- string or list; what was being searched for, either as a
      string which is taken verbatim, or a list of multiple search terms if
      there were more than one

    This function parses ``results``: each line should have either the form
    ``FILENAME`` or ``FILENAME: string`` where FILENAME is the file in which
    the string that matched the search was found.  If FILENAME ends in '.html',
    then this is part of the documentation; otherwise, it is in the source
    code.  In either case, an appropriate link is created.

    EXAMPLES::

        sage: from sage.misc.sagedoc import format_search_as_html
        sage: format_search_as_html('Source', 'algebras/steenrod_algebra_element.py:        an antihomomorphism: if we call the antipode `c`, then', 'antipode antihomomorphism')
        '<html><font color="black"><h2>Search Source: "antipode antihomomorphism"</h2></font><font color="darkpurple"><ol><li><a href="/src/algebras/steenrod_algebra_element.py" target="_blank"><tt>algebras/steenrod_algebra_element.py</tt></a>\n</ol></font></html>'
        sage: format_search_as_html('Other', 'html/en/reference/sage/algebras/steenrod_algebra_element.html:an antihomomorphism: if we call the antipode <span class="math">c</span>, then', 'antipode antihomomorphism')
        '<html><font color="black"><h2>Search Other: "antipode antihomomorphism"</h2></font><font color="darkpurple"><ol><li><a href="/doc/live/reference/sage/algebras/steenrod_algebra_element.html" target="_blank"><tt>reference/sage/algebras/steenrod_algebra_element.html</tt></a>\n</ol></font></html>'

    TESTS:

        Test that results from a ``search_src`` with ``multiline=True`` works
        reasonably::

            sage: results = search_src('format_search_as_html',  # long time
            ....:                      multiline=True, interact=False)
            sage: format_search_as_html('Source', results,  # long time
            ....:                       'format_search_as_html')
            '<html><font color="black"><h2>Search Source: "format_search_as_html"</h2></font><font color="darkpurple"><ol><li><a href="/src/misc/sagedoc.py" target="_blank"><tt>misc/sagedoc.py</tt></a>\n</ol></font></html>'
    """

    if not isinstance(search, list):
        search = [search]

    s = [
        '<html>',
        '<font color="black">',
        '<h2>Search {}: {}</h2>'.format(
            what, ', '.join('"{}"'.format(s) for s in search if s.strip())),
        '</font>',
        '<font color="darkpurple">',
        '<ol>'
    ]

    append = s.append

    if not isinstance(results, list):
        results = results.splitlines()

    files = set()
    for L in results:
        filename = L.strip().split(':', 1)[0]
        if filename:
            files.add(filename)
    files = sorted(files)
    for F in files:
        if F.endswith('.html'):
            F = F.split('/', 2)[2]
            url = '/doc/live/' + F
        else:
            # source code
            url = '/src/' + F
        append('<li><a href="%s" target="_blank"><tt>%s</tt></a>\n' % (url, F))
    append('</ol>')
    append('</font>')
    append('</html>')
    return ''.join(s)


#######################################
#    Add detex'ing of documentation
#######################################


def my_getsource(obj, oname=''):
    """
    Retrieve the source code for ``obj``.

    INPUT:

    - ``obj`` -- a Sage object, function, etc.

    - ``oname`` -- string (optional); a name under which the object is
      known. Currently ignored by Sage

    OUTPUT:

    Its documentation (string)

    EXAMPLES::

        sage: from sage.misc.sagedoc import my_getsource
        sage: s = my_getsource(identity_matrix)                                         # needs sage.modules
        sage: s[15:34]                                                                  # needs sage.modules
        'def identity_matrix'
    """
    try:
        s = sageinspect.sage_getsource(obj)
        return format_src(s)
    except Exception as msg:
        print('Error getting source:', msg)
        return None


class _sage_doc:
    """
    Open Sage documentation in a web browser, from either the
    command-line or the notebook.

    - Type "browse_sage_doc.DOCUMENT()" to open the named document --
      for example, "browse_sage_doc.tutorial()" opens the tutorial.
      Available documents are

      - tutorial: the Sage tutorial
      - reference: the Sage reference manual
      - constructions: "how do I construct ... in Sage?"
      - developer: the Sage developer's guide.

    - Type "browse_sage_doc(OBJECT, output=FORMAT, view=BOOL)" to view
      the documentation for OBJECT, as in
      "browse_sage_doc(identity_matrix, 'html').  ``output`` can be
      either 'html' or 'rst': the form of the output.  ``view`` is
      only relevant if ``output`` is ``html``; in this case, if
      ``view`` is ``True`` (its default value), then open up the
      documentation in a web browser.  Otherwise, just output the
      documentation as a string.

    EXAMPLES::

        sage: browse_sage_doc._open("reference", testing=True)[0]                       # needs sagemath_doc_html
        'http://localhost:8000/doc/live/reference/index.html'
        sage: browse_sage_doc(identity_matrix, 'rst')[-107:-47]                         # needs sage.modules
        '...Full MatrixSpace of 3 by 3 sparse matrices...'
    """
    def __init__(self):
        """
        EXAMPLES::

            sage: browse_sage_doc._base_url
            'http://localhost:8000/doc/live/'
        """
        self._base_url = "http://localhost:8000/doc/live/"
        self._base_path = os.path.join(SAGE_DOC, "html", "en")

    def __call__(self, obj, output='html', view=True):
        r"""
        Return the documentation for ``obj``.

        INPUT:

        - ``obj`` -- a Sage object
        - ``output`` -- 'html', 'rst', or 'text': return documentation in this form
        - ``view`` -- only has an effect if output is 'html': in this
          case, if ``view`` is ``True``, display the documentation in
          a web browser.  Otherwise, return the documentation as a
          string.

        EXAMPLES::

            sage: browse_sage_doc(identity_matrix, 'rst')                               # needs sage.modules
            "...**File:**...**Type:**...**Definition:** identity_matrix..."
            sage: identity_matrix.__doc__ in browse_sage_doc(identity_matrix, 'rst')    # needs sage.modules
            True
            sage: browse_sage_doc(identity_matrix, 'html', False)                       # needs sagemath_doc_html sphinx
            '...div...File:...Type:...Definition:...identity_matrix...'

        In the 'text' version, double colons have been replaced with
        single ones (among other things)::

            sage: '::' in browse_sage_doc(identity_matrix, 'rst')                       # needs sage.modules
            True
            sage: '::' in browse_sage_doc(identity_matrix, 'text')                      # needs sphinx
            False
        """
        if output != 'html' and view:
            view = False

        s = ''
        newline = "\n\n"  # blank line to start new paragraph

        try:
            filename = sageinspect.sage_getfile(obj)
            s += '**File:** %s' % filename
            s += newline
        except TypeError:
            pass

        obj_name = ''
        locs = sys._getframe(1).f_locals
        for var in locs:
            if id(locs[var]) == id(obj):
                obj_name = var

        s += '**Type:** %s' % type(obj)
        s += newline
        s += '**Definition:** %s' % sageinspect.sage_getdef(obj, obj_name)
        s += newline
        s += '**Docstring:**'
        s += newline
        s += sageinspect.sage_getdoc(obj, obj_name, embedded=True)

        # now s should be the reST version of the docstring
        if output == 'html':
            try:
                from .sphinxify import sphinxify
            except ImportError:
                from html import escape
                html = escape(s)
            else:
                html = sphinxify(s)
            if view:
                path = os.path.join(tmp_dir(), "temp.html")
                filed = open(path, 'w')

                static_path = os.path.join(SAGE_DOC, "html", "en", "_static")
                if os.path.exists(static_path):
                    title = obj_name + ' - Sage ' + sage.version.version + ' Documentation'
                    template = """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
  "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">

<html xmlns="http://www.w3.org/1999/xhtml">
  <head>
    <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
    <title>%(title)s</title>
    <link rel="stylesheet" href="%(static_path)s/default.css" type="text/css" />
    <link rel="stylesheet" href="%(static_path)s/pygments.css" type="text/css" />
    <style type="text/css">
      <!--
        div.body {
          margin: 1.0em;
          padding: 1.0em;
        }
        div.bodywrapper {
          margin: 0;
        }
      -->
    </style>
    <script type="text/javascript">
      var DOCUMENTATION_OPTIONS = {
        URL_ROOT:    '',
        VERSION:     '%(version)s',
        COLLAPSE_MODINDEX: false,
        FILE_SUFFIX: '.html',
        HAS_SOURCE:  false
      };
    </script>
    <script type="text/javascript" src="%(static_path)s/jquery.js"></script>
    <script type="text/javascript" src="%(static_path)s/doctools.js"></script>
    <link rel="shortcut icon" href="%(static_path)s/favicon.ico" />
    <link rel="icon" href="%(static_path)s/sageicon.png" type="image/x-icon" />
  </head>
  <body>
    <div class="document">
      <div class="documentwrapper">
        <div class="bodywrapper">
          <div class="body">
            %(html)s
          </div>
        </div>
      </div>
    </div>
  </body>
</html>"""
                    html = template % {'html': html,
                                       'static_path': static_path,
                                       'title': title,
                                       'version': sage.version.version}

                filed.write(html)
                filed.close()
                os.system(browser() + " " + path)
            else:
                return html
        elif output == 'rst':
            return s
        elif output == 'text':
            try:
                from .sphinxify import sphinxify
            except ImportError:
                return s
            else:
                return sphinxify(s, format='text')
        else:
            raise ValueError("output type {} not recognized".format(output))

    def _open(self, name, testing=False):
        """
        Open the document ``name`` in a web browser.  This constructs
        the appropriate URL and/or path name and passes it to the web
        browser.

        INPUT:

        - ``name`` -- string; name of the documentation

        - ``testing`` -- boolean (default: ``False``); if ``True``,
          then just return the URL and path-name for this document
          (don't open the web browser)

        EXAMPLES::

            sage: browse_sage_doc._open("reference", testing=True)[0]                   # needs sagemath_doc_html
            'http://localhost:8000/doc/live/reference/index.html'
            sage: browse_sage_doc._open("tutorial", testing=True)[1]                    # needs sagemath_doc_html
            '.../html/en/tutorial/index.html'
        """
        url = self._base_url + os.path.join(name, "index.html")
        path = os.path.join(self._base_path, name, "index.html")
        if not os.path.exists(path):
            raise OSError("""The document '{0}' does not exist.  Please build it
with 'sage -docbuild {0} html' and try again.""".format(name))

        if testing:
            return (url, path)

        os.system(browser() + " " + path)

    def tutorial(self):
        """
        The Sage tutorial.  To get started with Sage, start here.

        EXAMPLES::

            sage: tutorial()  # indirect doctest, not tested
        """
        self._open("tutorial")

    def reference(self):
        """
        The Sage reference manual.

        EXAMPLES::

            sage: reference()  # indirect doctest, not tested
            sage: manual()  # indirect doctest, not tested
        """
        self._open("reference")

    manual = reference

    def developer(self):
        """
        The Sage developer's guide.  Learn to develop programs for Sage.

        EXAMPLES::

            sage: developer()  # indirect doctest, not tested
        """
        self._open("developer")

    def constructions(self):
        """
        Sage constructions.  Attempts to answer the question "How do I
        construct ... in Sage?"

        EXAMPLES::

            sage: constructions()  # indirect doctest, not tested
        """
        self._open("constructions")


browse_sage_doc = _sage_doc()
tutorial = browse_sage_doc.tutorial
reference = browse_sage_doc.reference
manual = browse_sage_doc.reference
developer = browse_sage_doc.developer
constructions = browse_sage_doc.constructions

python_help = pydoc.help


def help(module=None):
    """
    If there is an argument ``module``, print the Python help message
    for ``module``.  With no argument, print a help message about
    getting help in Sage.

    EXAMPLES::

        sage: help()
        Welcome to Sage ...
    """
    if module is not None:
        python_help(module)
    else:
        print("""Welcome to Sage {}!

To view the Sage tutorial in your web browser, type "tutorial()", and
to view the (very detailed) Sage reference manual, type "manual()".
For help on any Sage function, for example "matrix_plot", type
"matrix_plot?" to see a help message, type "help(matrix_plot)" to see
a very similar message, type "browse_sage_doc(matrix_plot)" to view a
help message in a web browser, and type "matrix_plot??" to look at the
function's source code.

(When you type something like "matrix_plot?", "help(matrix_plot)", or
"matrix_plot??", Sage may start a paging program to display the
requested message. Type a space to scroll to the next page, type "h"
to get help on the paging program, and type "q" to quit it and return
to the "sage:" prompt.)

For license information for Sage and its components, read the file
"COPYING.txt" in the top-level directory of the Sage installation,
or type "license()".

To enter Python's interactive online help utility, type "python_help()".
To get help on a Python function, module or package, type "help(MODULE)" or
"python_help(MODULE)".""".format(sage.version.version))
