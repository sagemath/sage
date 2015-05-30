"""
HTML typesetting for the notebook
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

def math_parse(s):
    r"""
    Turn the HTML-ish string s that can have \$\$ and \$'s in it into
    pure HTML.  See below for a precise definition of what this means.

    INPUT:

    - ``s`` -- a string

    OUTPUT:

    - a string.

    Do the following:

    * Replace all ``\$ text \$``\'s by
      ``<script type="math/tex"> text </script>``
    * Replace all ``\$\$ text \$\$``\'s by
      ``<script type="math/tex; mode=display"> text </script>``
    * Replace all ``\ \$``\'s by ``\$``\'s.  Note that in
      the above two cases nothing is done if the ``\$``
      is preceeded by a backslash.
    * Replace all ``\[ text \]``\'s by
      ``<script type="math/tex; mode=display"> text </script>``

    EXAMPLES::

        sage: sage.misc.html.math_parse('This is $2+2$.')
        'This is <script type="math/tex">2+2</script>.'
        sage: sage.misc.html.math_parse('This is $$2+2$$.')
        'This is <script type="math/tex; mode=display">2+2</script>.'
        sage: sage.misc.html.math_parse('This is \\[2+2\\].')
        'This is <script type="math/tex; mode=display">2+2</script>.'
        sage: sage.misc.html.math_parse(r'This is \[2+2\].')
        'This is <script type="math/tex; mode=display">2+2</script>.'

    TESTS::

        sage: sage.misc.html.math_parse(r'This \$\$is $2+2$.')
        'This $$is <script type="math/tex">2+2</script>.'
    """
    # first replace \\[ and \\] by <script type="math/tex; mode=display">
    # and </script>, respectively.
    while True:
        i = s.find('\\[')
        if i == -1:
            break
        else:
            s = s[:i] + '<script type="math/tex; mode=display">' + s[i+2:]
            j = s.find('\\]')
            if j == -1:  # missing right-hand delimiter, so add one
                s = s + '</script>'
            else:
                s = s[:j] + '</script>' + s[j+2:]

    # Below t always has the "parsed so far" version of s, and s is
    # just the part of the original input s that hasn't been parsed.
    t = ''
    while True:
        i = s.find('$')
        if i == -1:
            # No dollar signs -- definitely done.
            return t + s
        elif i > 0 and s[i-1] == '\\':
            # A dollar sign with a backslash right before it, so
            # we ignore it by sticking it in the parsed string t
            # and skip to the next iteration.
            t += s[:i-1] + '$'
            s = s[i+1:]
            continue
        elif i+1 < len(s) and s[i+1] == '$':
            # Found a math environment. Double dollar sign so display mode.
            disp = '; mode=display'
        else:
            # Found math environment. Single dollar sign so default mode.
            disp = ''

        # Now find the matching $ sign and form the html string.

        if len(disp) > 0:
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

        t += s[:i] + '<script type="math/tex%s">%s</script>'%(disp,
                      ' '.join(txt.splitlines()))
        s = s[j+1:]
        if len(disp) > 0:
            s = s[1:]
    return t

class HTMLExpr(str):
    r"""
    A class for HTML expression
    """
    def __repr__(self):
        return str(self)

class HTML:
    def __call__(self, s, globals=None, locals=None):
        """
        Display the given HTML expression in the notebook.

        INPUT:

        - ``s`` -- a string

        OUTPUT:

        - prints a code that embeds HTML in the output.

        By default in the notebook an output cell has two parts, first a plain
        text preformat part, then second a general HTML part (not pre).  If
        you call html(s) at any point then that adds something that will be
        displayed in the preformated part in html.

        EXAMPLES::

            sage: html('<a href="http://sagemath.org">sagemath</a>')
            <html><font color='black'><a href="http://sagemath.org">sagemath</a></font></html>
            sage: html('<hr>')
            <html><font color='black'><hr></font></html>
        """
        return HTMLExpr(self.eval(s, globals, locals))

    def eval(self, s, globals=None, locals=None):
        r"""
        Return an html representation for an object ``s``.

        If ``s`` has a method ``_html_()``, call that. Otherwise, call
        :func:`math_parse` on ``str(s)``, evaluate any variables in
        the result, and add some html preamble and postamble.

        In any case, *print* the resulting html string. This method
        always *returns* an empty string.

        EXAMPLES::

            sage: html.eval('<hr>')
            <html><font color='black'><hr></font></html>
            ''
        """
        if hasattr(s, '_html_'):
            s._html_()
            return ''
        if globals is None:
            globals = {}
        if locals is None:
            locals = {}
        s = str(s)
        s = math_parse(s)
        t = ''
        while len(s) > 0:
            i = s.find('<sage>')
            if i == -1:
                 t += s
                 break
            j = s.find('</sage>')
            if j == -1:
                 t += s
                 break
            t += s[:i] + '<script type="math/tex">%s</script>'%\
                     latex(sage_eval(s[6+i:j], locals=locals))
            s = s[j+7:]
        print("<html><font color='black'>{}</font></html>".format(t))
        return ''

    def table(self, x, header = False):
        r"""
        Print a nested list as a HTML table.  Strings of html
        will be parsed for math inside dollar and double-dollar signs.
        2D graphics will be displayed in the cells.  Expressions will
        be latexed.


        INPUT:

        - ``x`` -- a list of lists (i.e., a list of table rows)
        - ``header`` -- a row of headers.  If ``True``, then the first
          row of the table is taken to be the header.

        EXAMPLES::

            sage: html.table([(i, j, i == j) for i in [0..1] for j in [0..1]])
            <html>
            <div class="notruncate">
            <table class="table_form">
            <tbody>
            <tr class ="row-a">
            <td><script type="math/tex">0</script></td>
            <td><script type="math/tex">0</script></td>
            <td><script type="math/tex">\mathrm{True}</script></td>
            </tr>
            <tr class ="row-b">
            <td><script type="math/tex">0</script></td>
            <td><script type="math/tex">1</script></td>
            <td><script type="math/tex">\mathrm{False}</script></td>
            </tr>
            <tr class ="row-a">
            <td><script type="math/tex">1</script></td>
            <td><script type="math/tex">0</script></td>
            <td><script type="math/tex">\mathrm{False}</script></td>
            </tr>
            <tr class ="row-b">
            <td><script type="math/tex">1</script></td>
            <td><script type="math/tex">1</script></td>
            <td><script type="math/tex">\mathrm{True}</script></td>
            </tr>
            </tbody>
            </table>
            </div>
            </html>

            sage: html.table([(x,n(sin(x), digits=2)) for x in [0..3]], header = ["$x$", "$\sin(x)$"])
            <html>
            <div class="notruncate">
            <table class="table_form">
            <tbody>
            <tr>
            <th><script type="math/tex">x</script></th>
            <th><script type="math/tex">\sin(x)</script></th>
            </tr>
            <tr class ="row-a">
            <td><script type="math/tex">0</script></td>
            <td><script type="math/tex">0.00</script></td>
            </tr>
            <tr class ="row-b">
            <td><script type="math/tex">1</script></td>
            <td><script type="math/tex">0.84</script></td>
            </tr>
            <tr class ="row-a">
            <td><script type="math/tex">2</script></td>
            <td><script type="math/tex">0.91</script></td>
            </tr>
            <tr class ="row-b">
            <td><script type="math/tex">3</script></td>
            <td><script type="math/tex">0.14</script></td>
            </tr>
            </tbody>
            </table>
            </div>
            </html>

        """
        from table import table
        table(x, header_row=header)._html_()

    def iframe(self, url, height=400, width=800):
        r"""
        Put an existing web page into a worksheet.

        INPUT:

        - ``url`` -- a url string, either with or without URI scheme
          (defaults to "http").
        - ``height`` -- the number of pixels for the page height.
          Defaults to 400.
        - ``width`` -- the number of pixels for the page width.
          Defaults to 800.

        OUTPUT:

            Opens the url in a worksheet. If the url is a regular web page it
            will appear in the worksheet. This was originally intended to bring
            GeoGebra worksheets into Sage, but it can be used for many other
            purposes.

        EXAMPLES::

            sage: html.iframe("sagemath.org")
            <html><font color='black'><iframe height="400" width="800"
            src="http://sagemath.org"></iframe></font></html>
            sage: html.iframe("http://sagemath.org",30,40)
            <html><font color='black'><iframe height="30" width="40"
            src="http://sagemath.org"></iframe></font></html>
            sage: html.iframe("https://sagemath.org",30)
            <html><font color='black'><iframe height="30" width="800"
            src="https://sagemath.org"></iframe></font></html>
            sage: html.iframe("/home/admin/0/data/filename")
            <html><font color='black'><iframe height="400" width="800"
            src="/home/admin/0/data/filename"></iframe></font></html>
            sage: html.iframe('data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAA'
            ... 'AUAAAAFCAYAAACNbyblAAAAHElEQVQI12P4//8/w38GIAXDIBKE0DHxgljNBA'
            ... 'AO9TXL0Y4OHwAAAABJRU5ErkJggg=="')
            <html><font color='black'><iframe height="400" width="800"
            src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAUAAAAFCAYAAACNbyblAAAAHElEQVQI12P4//8/w38GIAXDIBKE0DHxgljNBAAO9TXL0Y4OHwAAAABJRU5ErkJggg==""></iframe></font></html>

        AUTHOR:

        - Bruce Cohen (2011-06-14)
        """
        if ":" not in url and not url.startswith('/'):
            url = "http://" + url
        string = ( '<iframe height="%d" width="%d" src="%s"></iframe>' %
                    (height, width, url) )
        return html(string)

html = HTML()
# Ensure that html appear in the sphinx doc as a function
# so that the link :func:`html` is correctly set up.
html.__doc__ = HTML.__call__.__doc__
