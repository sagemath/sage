# -*- encoding: utf-8 -*-
r"""
Rich Output for the Browser
"""

from sage.repl.rich_output.output_basic import OutputBase
from sage.repl.rich_output.buffer import OutputBuffer


class OutputHtml(OutputBase):

    def __init__(self, html):
        """
        HTML Output
        
        INPUT:

        - ``html`` --
          :class:`~sage.repl.rich_output.buffer.OutputBuffer`. Alternatively,
          a string (bytes) can be passed directly which will then be
          converted into an
          :class:`~sage.repl.rich_output.buffer.OutputBuffer`. String
          containing the html fragment code. Excludes the surrounding
          ``<body>`` and ``<html>`` tag.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputHtml
            sage: OutputHtml('<div>Foo<b>B</b>ar</div>')
            OutputHtml container
        """
        self.html = OutputBuffer(html)

    @classmethod
    def example(cls):
        r"""
        Construct a sample Html output container

        This static method is meant for doctests, so they can easily
        construt an example.

        OUTPUT:

        An instance of :class:`OutputHtml`.
        
        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputHtml
            sage: OutputHtml.example()
            OutputHtml container
            sage: OutputHtml.example().html.get()
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
        print(self.html.get())

    def with_html_tag(self):
        r"""
        Return the HTML code surrounded by ``<html>`` tag

        This is just a convenience method.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputHtml
            sage: rich_output = OutputHtml.example()
            sage: rich_output.print_to_stdout()
            <div>Hello World!</div>
            sage: rich_output.with_html_tag()
            '<html><div>Hello World!</div></html>'
        """
        return '<html>{0}</html>'.format(self.html.get())
        
    
