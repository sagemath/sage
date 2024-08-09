# sage_setup: distribution = sagemath-repl
r"""
Rich Output for the Browser
"""

import re

from sage.repl.rich_output.output_basic import OutputBase
from sage.repl.rich_output.buffer import OutputBuffer

# regex to match "<html>\[...\]</html>" or "<html>\(...\)</html>"
latex_re = re.compile(r'<html>(?P<mathstart>\\\[|\\\()(?P<latex>.*)(?P<mathend>\\\]|\\\))</html>',
                      flags=re.DOTALL)

class OutputHtml(OutputBase):

    def __init__(self, html):
        """
        HTML Output.

        INPUT:

        - ``html`` --
          :class:`~sage.repl.rich_output.buffer.OutputBuffer`. Alternatively, a
          string (bytes) can be passed directly which will then be converted
          into an :class:`~sage.repl.rich_output.buffer.OutputBuffer`. String
          containing the html fragment code. Excludes the surrounding
          ``<body>`` and ``<html>`` tag.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputHtml
            sage: OutputHtml('<div>Foo<b>B</b>ar</div>')
            OutputHtml container
        """
        self.html = OutputBuffer(html)

        # if the html is a simple wrapper of latex for mathjax rendering, then
        # the latex string is saved for possible latex output such as Jupyter's
        # pdf export of a notebook
        m = latex_re.match(html)
        if m:
            mathjax_string = m.group('latex')
            latex_string = mathjax_string.replace('&lt;', '<')
            if m.group('mathstart') == r'\[' and m.group('mathend') == r'\]':
                self.latex = OutputBuffer('$$' + latex_string + '$$')
            else:
                self.latex = OutputBuffer('$' + latex_string + '$')
        else:
            self.latex = None

    @classmethod
    def example(cls):
        r"""
        Construct a sample Html output container.

        This static method is meant for doctests, so they can easily
        construct an example.

        OUTPUT: an instance of :class:`OutputHtml`

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputHtml
            sage: OutputHtml.example()
            OutputHtml container
            sage: OutputHtml.example().html.get_str()
            '<div>Hello World!</div>'
        """
        return cls('<div>Hello World!</div>')

    def print_to_stdout(self):
        r"""
        Write the data to stdout.

        This is just a convenience method to help with debugging.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputHtml
            sage: rich_output = OutputHtml.example()
            sage: rich_output.print_to_stdout()
            <div>Hello World!</div>
        """
        print(self.html.get_unicode())

    def with_html_tag(self):
        r"""
        Return the HTML code surrounded by ``<html>`` tag.

        This is just a convenience method.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputHtml
            sage: rich_output = OutputHtml.example()
            sage: rich_output.print_to_stdout()
            <div>Hello World!</div>
            sage: rich_output.with_html_tag()
            '<html><div>Hello World!</div></html>'
        """
        return '<html>{0}</html>'.format(self.html.get_unicode())
