"""
HTML Fragments

This module defines a HTML fragment class, which holds a piece of
HTML. This is primarily used in browser-based notebooks, though it
might be useful for creating static pages as well.

This module defines :class:`MathJax`, an object of which performs the task of
producing an HTML representation of any object. The produced HTML is
renderable in a browser-based notebook with the help of MathJax.
"""

#*****************************************************************************
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.latex import latex
from sage.misc.sage_eval import sage_eval
from sage.structure.sage_object import SageObject


class HtmlFragment(str, SageObject):
    r"""
    A HTML fragment.

    This is a piece of HTML, usually not a complete document.  For
    example, just a ``<div>...</div>`` piece and not the entire
    ``<html>...</html>``.

    EXAMPLES::

        sage: from sage.misc.html import HtmlFragment
        sage: HtmlFragment('<b>test</b>')
        <b>test</b>

    .. automethod:: _rich_repr_
    """

    def _rich_repr_(self, display_manager, **kwds):
        """
        Rich Output Magic Method

        See :mod:`sage.repl.rich_output` for details.

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: h = sage.misc.html.HtmlFragment('<b>old</b>')
            sage: h._rich_repr_(dm)  # the doctest backend does not support html
            OutputPlainText container
        """
        OutputHtml = display_manager.types.OutputHtml
        if OutputHtml in display_manager.supported_output():
            return OutputHtml(self)
        else:
            return display_manager.types.OutputPlainText(self)


def math_parse(s):
    r"""
    Transform the string ``s`` with TeX maths to an HTML string renderable by
    MathJax.

    INPUT:

    - ``s`` -- a string

    OUTPUT:

    A :class:`HtmlFragment` instance.

    Specifically this method does the following:

    * Replace all ``\$text\$``\'s by ``\(text\)``
    * Replace all ``\$\$text\$\$``\'s by ``\[text\]``
    * Replace all ``\\\$``\'s by ``\$``\'s. Note that this has precedence over
      the above two cases.

    EXAMPLES::

        sage: print(sage.misc.html.math_parse('This is $2+2$.'))
        This is \(2+2\).
        sage: print(sage.misc.html.math_parse('This is $$2+2$$.'))
        This is \[2+2\].
        sage: print(sage.misc.html.math_parse('This is \\[2+2\\].'))
        This is \[2+2\].
        sage: print(sage.misc.html.math_parse(r'\$2+2\$ is rendered to $2+2$.'))
        <span>$</span>2+2<span>$</span> is rendered to \(2+2\).

    """
    # Below t always has the "parsed so far" version of s, and s is
    # just the part of the original input s that hasn't been parsed.
    t = ''
    while True:
        i = s.find('$')
        if i == -1:
            # No dollar signs -- definitely done.
            return HtmlFragment(t + s)
        elif i > 0 and s[i-1] == '\\':
            # A dollar sign with a backslash right before it, so this is a
            # normal dollar sign. If processEscapes is enabled in MathJax, "\$"
            # will do the job. But as we do not assume that, we use the span
            # tag safely.
            t += s[:i-1] + '<span>$</span>'
            s = s[i+1:]
            continue
        elif i+1 < len(s) and s[i+1] == '$':
            # Found a math environment. Double dollar sign so display mode.
            disp = True
        else:
            # Found math environment. Single dollar sign so default mode.
            disp = False

        # Now find the matching $ sign and form the html string.

        if disp:
            j = s[i+2:].find('$$')
            if j == -1:
                j = len(s)
                s += '$$'
            else:
                j += i + 2
            txt = s[i+2:j]
        else:
            j = s[i+2:].find('$')
            if j == -1:
                j = len(s)
                s += '$'
            else:
                j += i + 2
            txt = s[i+1:j]

        if disp:
            t += s[:i] + r'\[{0}\]'.format(' '.join(txt.splitlines()))
        else:
            t += s[:i] + r'\({0}\)'.format(' '.join(txt.splitlines()))
        s = s[j+1:]
        if disp:
            s = s[1:]
    return HtmlFragment(t)


class MathJaxExpr:
    """
    An arbitrary MathJax expression that can be nicely concatenated.

    EXAMPLES::

        sage: from sage.misc.html import MathJaxExpr
        sage: MathJaxExpr("a^{2}") + MathJaxExpr("x^{-1}")
        a^{2}x^{-1}
    """
    def __init__(self, y):
        """
        Initialize a MathJax expression.

        INPUT:

        - ``y`` - a string

        Note that no error checking is done on the type of ``y``.

        EXAMPLES::

            sage: from sage.misc.html import MathJaxExpr
            sage: jax = MathJaxExpr(3); jax  # indirect doctest
            3
            sage: TestSuite(jax).run(skip ="_test_pickling")
        """
        self.__y = y

    def __repr__(self):
        """
        Print representation.

        EXAMPLES::

            sage: from sage.misc.html import MathJaxExpr
            sage: jax = MathJaxExpr('3')
            sage: jax.__repr__()
            '3'
        """
        return str(self.__y)

    def __add__(self, y):
        """
        'Add' MathJaxExpr ``self`` to ``y``.  This concatenates them
        (assuming that they're strings).

        EXAMPLES::

            sage: from sage.misc.html import MathJaxExpr
            sage: j3 = MathJaxExpr('3')
            sage: jx = MathJaxExpr('x')
            sage: j3 + jx
            3x
        """
        return MathJaxExpr(self.__y + y)

    def __radd__(self, y):
        """
        'Add' MathJaxExpr ``y`` to ``self``.  This concatenates them
        (assuming that they're strings).

        EXAMPLES::

            sage: from sage.misc.html import MathJaxExpr
            sage: j3 = MathJaxExpr('3')
            sage: jx = MathJaxExpr('x')
            sage: j3.__radd__(jx)
            x3
        """
        return MathJaxExpr(y + self.__y)


class MathJax:
    r"""
    Render LaTeX input using MathJax.  This returns a :class:`MathJaxExpr`.

    EXAMPLES::

        sage: from sage.misc.html import MathJax
        sage: MathJax()(3)
        <html>\[\newcommand{\Bold}[1]{\mathbf{#1}}3\]</html>
        sage: MathJax()(ZZ)
        <html>\[\newcommand{\Bold}[1]{\mathbf{#1}}\Bold{Z}\]</html>
    """

    def __call__(self, x, combine_all=False):
        r"""
        Render LaTeX input using MathJax.  This returns a :class:`MathJaxExpr`.

        INPUT:

        - ``x`` - a Sage object

        - ``combine_all`` - boolean (Default: ``False``): If ``combine_all`` is
          ``True`` and the input is a tuple, then it does not return a tuple
          and instead returns a string with all the elements separated by
          a single space.

        OUTPUT:

        A :class:`MathJaxExpr`

        EXAMPLES::

            sage: from sage.misc.html import MathJax
            sage: MathJax()(3)
            <html>\[\newcommand{\Bold}[1]{\mathbf{#1}}3\]</html>
            sage: str(MathJax().eval(ZZ['x'], mode='display')) == str(MathJax()(ZZ['x']))
            True
        """
        return self.eval(x, combine_all=combine_all)

    def eval(self, x, globals=None, locals=None, mode='display', combine_all=False):
        r"""
        Render LaTeX input using MathJax.  This returns a :class:`MathJaxExpr`.

        INPUT:

        - ``x`` - a Sage object

        -  ``globals`` - a globals dictionary

        -  ``locals`` - extra local variables used when
           evaluating Sage code in ``x``.

        - ``mode`` - string (optional, default ``'display'``):
           ``'display'`` for displaymath, ``'inline'`` for inline
           math, or ``'plain'`` for just the LaTeX code without the
           surrounding html and script tags.

        - ``combine_all`` - boolean (Default: ``False``): If ``combine_all`` is
          ``True`` and the input is a tuple, then it does not return a tuple
          and instead returns a string with all the elements separated by
          a single space.

        OUTPUT:

        A :class:`MathJaxExpr`

        EXAMPLES::

            sage: from sage.misc.html import MathJax
            sage: MathJax().eval(3, mode='display')
            <html>\[\newcommand{\Bold}[1]{\mathbf{#1}}3\]</html>
            sage: MathJax().eval(3, mode='inline')
            <html>\(\newcommand{\Bold}[1]{\mathbf{#1}}3\)</html>
            sage: MathJax().eval(type(3), mode='inline')
            <html>\(\newcommand{\Bold}[1]{\mathbf{#1}}\verb|&lt;class|\verb| |\verb|'sage.rings.integer.Integer'>|\)</html>
        """
        # Get a regular LaTeX representation of x
        x = latex(x, combine_all=combine_all)

        # The "\text{\texttt{...}}" blocks are reformed to be renderable by MathJax.
        # These blocks are produced by str_function() defined in sage.misc.latex.
        prefix = r"\text{\texttt{"
        parts = x.split(prefix)
        for i, part in enumerate(parts):
            if i == 0:
                continue    # Nothing to do with the head part
            n = 1
            for closing, c in enumerate(part):
                if c == "{" and part[closing - 1] != "\\":
                    n += 1
                if c == "}" and part[closing - 1] != "\\":
                    n -= 1
                if n == -1:
                    break
            # part should end in "}}", so omit the last two characters
            # from y
            y = part[:closing-1]
            for delimiter in r"""|"'`#%&,.:;?!@_~^+-/\=<>()[]{}0123456789E""":
                if delimiter not in y:
                    break
            if delimiter == "E":
                # y is too complicated
                delimiter = "|"
                y = "(complicated string)"
            wrapper = r"\verb" + delimiter + "%s" + delimiter
            y = y.replace("{ }", " ").replace("{-}", "-")
            for c in r"#$%&\^_{}~":
                char_wrapper = r"{\char`\%s}" % c
                y = y.replace(char_wrapper, c)
            subparts = []
            nspaces = 0
            for subpart in y.split(" "):
                if subpart == "":
                    nspaces += 1
                    continue
                if nspaces > 0:
                    subparts.append(wrapper % (" " * nspaces))
                nspaces = 1
                subparts.append(wrapper % subpart)
            subparts.append(part[closing + 1:])
            parts[i] = "".join(subparts)

        from sage.misc.latex_macros import sage_configurable_latex_macros
        from sage.misc.latex import _Latex_prefs
        latex_string = ''.join(
            sage_configurable_latex_macros +
            [_Latex_prefs._option['macros']] +
            parts
        )
        mathjax_string = latex_string.replace('<', '&lt;')
        if mode == 'display':
            html = r'<html>\[{0}\]</html>'
        elif mode == 'inline':
            html = r'<html>\({0}\)</html>'
        elif mode == 'plain':
            return mathjax_string
        else:
            raise ValueError("mode must be either 'display', 'inline', or 'plain'")
        return MathJaxExpr(html.format(mathjax_string))


class HTMLFragmentFactory(SageObject):

    def _repr_(self):
        """
        Return string representation

        OUTPUT:

        String.

        EXAMPLES::

            sage: html
            Create HTML output (see html? for details)
        """
        return 'Create HTML output (see html? for details)'

    def __call__(self, obj, concatenate=True, strict=False):
        r"""
        Construct a HTML fragment

        INPUT:

        - ``obj`` -- anything. An object for which you want an HTML
          representation.

        - ``concatenate`` -- if ``True``, combine HTML representations of
          elements of the container ``obj``

        - ``strict`` -- if ``True``, construct an HTML representation of
          ``obj`` even if ``obj`` is a string

        OUTPUT:

        A :class:`HtmlFragment` instance.

        EXAMPLES::

            sage: h = html('<hr>');  pretty_print(h)
            <hr>
            sage: type(h)
            <class 'sage.misc.html.HtmlFragment'>

            sage: html(1/2)
            <html>\[\newcommand{\Bold}[1]{\mathbf{#1}}\frac{1}{2}\]</html>

            sage: html('<a href="http://sagemath.org">sagemath</a>')
            <a href="http://sagemath.org">sagemath</a>

            sage: html('<a href="http://sagemath.org">sagemath</a>', strict=True)
            <html>\[\newcommand{\Bold}[1]{\mathbf{#1}}\verb|&lt;a|\verb| |\verb|href="http://sagemath.org">sagemath&lt;/a>|\]</html>
        """
        # string obj is interpreted as an HTML in not strict mode
        if isinstance(obj, str) and not strict:
            return HtmlFragment(math_parse(obj))

        # prefer dedicated _html_() method
        try:
            result = obj._html_()
            return HtmlFragment(result)
        except AttributeError:
            pass

        # otherwise convert latex to html
        if concatenate:
            if isinstance(obj, (tuple, list)):
                obj = tuple(obj)
            result = MathJax().eval(obj, mode='display', combine_all=True)
        else:
            result = MathJax().eval(obj, mode='display', combine_all=False)
        return HtmlFragment(result)

    def eval(self, s, locals=None):
        r"""
        Evaluate embedded <sage> tags

        INPUT:

        - ``s`` -- string.

        - ``globals`` -- dictionary. The global variables when
          evaluating ``s``. Default: the current global variables.

        OUTPUT:

        A :class:`HtmlFragment` instance.

        EXAMPLES::

            sage: a = 123
            sage: html.eval('<sage>a</sage>')
            \(123\)
            sage: html.eval('<sage>a</sage>', locals={'a': 456})
            \(456\)
        """
        if locals is None:
            from sage.repl.user_globals import get_globals
            locals = get_globals()
        s = str(s)
        s = math_parse(s)
        t = ''
        while s:
            i = s.find('<sage>')
            if i == -1:
                 t += s
                 break
            j = s.find('</sage>')
            if j == -1:
                 t += s
                 break
            t += s[:i] + r'\({}\)'.format(latex(sage_eval(s[6+i:j], locals=locals)))
            s = s[j+7:]
        return HtmlFragment(t)

    def iframe(self, url, height=400, width=800):
        r"""
        Generate an iframe HTML fragment

        INPUT:

        - ``url`` -- string. A url, either with or without URI scheme
          (defaults to "http"), or an absolute file path.

        - ``height`` -- the number of pixels for the page height.
          Defaults to 400.

        - ``width`` -- the number of pixels for the page width.
          Defaults to 800.

        OUTPUT:

        A :class:`HtmlFragment` instance.

        EXAMPLES::

            sage: pretty_print(html.iframe("sagemath.org"))
            <iframe height="400" width="800"
            src="http://sagemath.org"></iframe>
            sage: pretty_print(html.iframe("http://sagemath.org",30,40))
            <iframe height="30" width="40"
            src="http://sagemath.org"></iframe>
            sage: pretty_print(html.iframe("https://sagemath.org",30))
            <iframe height="30" width="800"
            src="https://sagemath.org"></iframe>
            sage: pretty_print(html.iframe("/home/admin/0/data/filename"))
            <iframe height="400" width="800"
            src="file:///home/admin/0/data/filename"></iframe>
            sage: pretty_print(html.iframe('data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAA'
            ....: 'AUAAAAFCAYAAACNbyblAAAAHElEQVQI12P4//8/w38GIAXDIBKE0DHxgljNBA'
            ....: 'AO9TXL0Y4OHwAAAABJRU5ErkJggg=="'))
            <iframe height="400" width="800"
            src="http://data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAUAAAAFCAYAAACNbyblAAAAHElEQVQI12P4//8/w38GIAXDIBKE0DHxgljNBAAO9TXL0Y4OHwAAAABJRU5ErkJggg==""></iframe>
        """
        if url.startswith('/'):
            url = 'file://{0}'.format(url)
        elif '://' not in url:
            url = 'http://{0}'.format(url)
        return HtmlFragment('<iframe height="{0}" width="{1}" src="{2}"></iframe>'
                            .format(height, width, url))


html = HTMLFragmentFactory()
html.__doc__ = HTMLFragmentFactory.__call__.__doc__


def pretty_print_default(enable=True):
    r"""
    Enable or disable default pretty printing.

    Pretty printing means rendering things in HTML and by MathJax so that a
    browser-based frontend can render real math.

    This function is pretty useless without the notebook, it should not
    be in the global namespace.

    INPUT:

    -  ``enable`` -- bool (optional, default ``True``).  If ``True``, turn on
       pretty printing; if ``False``, turn it off.

    EXAMPLES::

        sage: pretty_print_default(True)
        sage: 'foo'  # the doctest backend does not support html
        'foo'
        sage: pretty_print_default(False)
        sage: 'foo'
        'foo'
    """
    from sage.repl.rich_output import get_display_manager
    dm = get_display_manager()
    dm.preferences.text = 'latex' if enable else None
