# sage_setup: distribution = sagemath-repl
r"""
Basic Output Types

The Sage rich representation system requires a special container class
to hold the data for each type of rich output. They all inherit from
:class:`OutputBase`, though a more typical example is
:class:`OutputPlainText`. Some output classes consist of more than one
data buffer, for example jmol or certain animation formats. The output
class is independent of user preferences and of the display
backend.

The display backends can define derived classes to attach
backend-specific display functionality to, for example how to launch a
viewer. But they must not change how the output container is
created. To enforce this, the Sage ``_rich_repr_`` magic method will
only ever see the output class defined here. The display manager will
promote it to a backend-specific subclass if necessary prior to
displaying it.

To create new types of output, you must create your own subclass of
:class:`OutputBase` and register it in
:mod:`sage.repl.rich_output.output_catalog`.

.. warning::

    All rich output data in subclasses of :class:`OutputBase` must be
    contained in :class:`~sage.repl.rich_output.buffer.OutputBuffer`
    instances. You must never reference any files on the local file
    system, as there is no guarantee that the notebook server and the
    worker process are on the same computer. Or even share a common
    file system.
"""

#*****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.structure.sage_object import SageObject
from sage.repl.rich_output.buffer import OutputBuffer


class OutputBase(SageObject):
    """
    Base class for all rich output containers.
    """

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.repl.rich_output.output_basic import OutputBase
            sage: output_base = OutputBase()
            sage: output_base._repr_()
            'OutputBase container'
        """
        return '{0} container'.format(self.__class__.__name__)

    @classmethod
    def example(cls):
        """
        Construct a sample instance.

        This static method is meant for doctests, so they can easily
        construct an example.

        OUTPUT: an instance of the :class:`OutputBase` subclass

        EXAMPLES::

            sage: from sage.repl.rich_output.output_basic import OutputBase
            sage: OutputBase.example()
            Traceback (most recent call last):
            ...
            NotImplementedError: derived classes must implement this class method
        """
        raise NotImplementedError('derived classes must implement this class method')


class OutputPlainText(OutputBase):

    def __init__(self, plain_text):
        """
        Plain Text Output.

        INPUT:

        - ``plain_text`` --
          :class:`~sage.repl.rich_output.buffer.OutputBuffer`. Alternatively,
          a bytes (string in Python 2.x) or string (unicode in Python
          2.x) can be passed directly which will then be converted
          into an
          :class:`~sage.repl.rich_output.buffer.OutputBuffer`. The
          plain text output.

        This should always be exactly the same as the (non-rich)
        output from the ``_repr_`` method. Every backend object must
        support plain text output as fallback.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputPlainText
            sage: OutputPlainText('foo')
            OutputPlainText container
        """
        self.text = OutputBuffer(plain_text)

    @classmethod
    def example(cls):
        """
        Construct a sample plain text output container.

        This static method is meant for doctests, so they can easily
        construct an example.

        OUTPUT: an instance of :class:`OutputPlainText`

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputPlainText
            sage: OutputPlainText.example()
            OutputPlainText container
            sage: OutputPlainText.example().text.get_str()
            'Example plain text output'
        """
        return cls('Example plain text output')

    def print_to_stdout(self):
        """
        Write the data to stdout.

        This is just a convenience method to help with debugging.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputPlainText
            sage: plain_text = OutputPlainText.example()
            sage: plain_text.print_to_stdout()
            Example plain text output
        """
        print(self.text.get_str())


class OutputAsciiArt(OutputBase):

    def __init__(self, ascii_art):
        """
        ASCII Art Output.

        INPUT:

        - ``ascii_art`` --
          :class:`~sage.repl.rich_output.buffer.OutputBuffer`. Alternatively,
          a string (bytes) can be passed directly which will then be
          converted into an
          :class:`~sage.repl.rich_output.buffer.OutputBuffer`. Ascii
          art rendered into a string.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputAsciiArt
            sage: OutputAsciiArt(':-}')
            OutputAsciiArt container
        """
        self.ascii_art = OutputBuffer(ascii_art)

    @classmethod
    def example(cls):
        r"""
        Construct a sample ascii art output container.

        This static method is meant for doctests, so they can easily
        construct an example.

        OUTPUT: an instance of :class:`OutputAsciiArt`

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputAsciiArt
            sage: OutputAsciiArt.example()
            OutputAsciiArt container
            sage: OutputAsciiArt.example().ascii_art.get_str()
            '[                        *   *   *    * ]\n[      **   **   *    *  *   *  *    *  ]\n[ ***, * , *  , **, ** , *, * , * , *   ]'
        """
        return cls('[                        *   *   *    * ]\n'
                   '[      **   **   *    *  *   *  *    *  ]\n'
                   '[ ***, * , *  , **, ** , *, * , * , *   ]')

    def print_to_stdout(self):
        """
        Write the data to stdout.

        This is just a convenience method to help with debugging.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputAsciiArt
            sage: ascii_art = OutputAsciiArt.example()
            sage: ascii_art.print_to_stdout()
            [                        *   *   *    * ]
            [      **   **   *    *  *   *  *    *  ]
            [ ***, * , *  , **, ** , *, * , * , *   ]
        """
        print(self.ascii_art.get_str())


class OutputUnicodeArt(OutputBase):

    def __init__(self, unicode_art):
        """
        Unicode Art Output.

        Similar to :class:`OutputAsciiArt` but using the entire
        unicode range.

        INPUT:

        - ``unicode_art`` --
          :class:`~sage.repl.rich_output.buffer.OutputBuffer`. Alternatively,
          a string (unicode in Python 2.x) can be passed directly
          which will then be converted into an
          :class:`~sage.repl.rich_output.buffer.OutputBuffer`. Unicode
          art rendered into a string.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputUnicodeArt
            sage: OutputUnicodeArt(u':-}')
            OutputUnicodeArt container
        """
        # Internally, all buffers store bytes. Unicode is always utf-8
        # encoded.
        if not isinstance(unicode_art, bytes):
            unicode_art = unicode_art.encode('utf-8')
        self.unicode_art = OutputBuffer(unicode_art)

    @classmethod
    def example(cls):
        r"""
        Construct a sample unicode art output container.

        This static method is meant for doctests, so they can easily
        construct an example.

        OUTPUT: an instance of :class:`OutputUnicodeArt`

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputUnicodeArt
            sage: OutputUnicodeArt.example()
            OutputUnicodeArt container
            sage: print(OutputUnicodeArt.example().unicode_art.get_unicode())
            ⎛-11   0   1⎞
            ⎜  3  -1   0⎟
            ⎝ -1  -1   0⎠
        """
        return cls('⎛-11   0   1⎞\n'
                   '⎜  3  -1   0⎟\n'
                   '⎝ -1  -1   0⎠')

    def print_to_stdout(self):
        """
        Write the data to stdout.

        This is just a convenience method to help with debugging.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputUnicodeArt
            sage: unicode_art = OutputUnicodeArt.example()
            sage: unicode_art.print_to_stdout()
            ⎛-11   0   1⎞
            ⎜  3  -1   0⎟
            ⎝ -1  -1   0⎠
        """
        print(self.unicode_art.get_unicode())


class OutputLatex(OutputBase):

    def __init__(self, latex):
        """
        LaTeX Output.

        .. NOTE::

            The LaTeX commands will only use a subset of LaTeX that
            can be displayed by MathJax.

        INPUT:

        - ``latex`` --
          :class:`~sage.repl.rich_output.buffer.OutputBuffer`. Alternatively,
          a string (bytes) can be passed directly which will then be
          converted into an
          :class:`~sage.repl.rich_output.buffer.OutputBuffer`. String
          containing the latex equation code. Excludes the surrounding
          dollar signs / LaTeX equation environment.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputLatex
            sage: OutputLatex(latex(sqrt(x)))                                           # needs sage.symbolic
            OutputLatex container
        """
        self.latex = OutputBuffer(latex)

    def display_equation(self):
        r"""
        Return the LaTeX code for a display equation.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputLatex
            sage: rich_output = OutputLatex('1')
            sage: rich_output.latex
            buffer containing 1 bytes
            sage: rich_output.latex.get_str()
            '1'
            sage: rich_output.display_equation()
            '\\begin{equation}\n1\n\\end{equation}'
        """
        return '\n'.join([r'\begin{equation}', self.latex.get_str(),
                          r'\end{equation}'])

    def inline_equation(self):
        r"""
        Return the LaTeX code for an inline equation.

        OUTPUT: string

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputLatex
            sage: rich_output = OutputLatex('1')
            sage: rich_output.latex
            buffer containing 1 bytes
            sage: rich_output.latex.get_str()
            '1'
            sage: rich_output.inline_equation()
            '\\begin{math}\n1\n\\end{math}'
        """
        return '\n'.join([r'\begin{math}', self.latex.get_str(), r'\end{math}'])

    @classmethod
    def example(cls):
        r"""
        Construct a sample LaTeX output container.

        This static method is meant for doctests, so they can easily
        construct an example.

        OUTPUT: an instance of :class:`OutputLatex`

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputLatex
            sage: OutputLatex.example()
            OutputLatex container
            sage: OutputLatex.example().latex.get_str()
            '\\newcommand{\\Bold}[1]{\\mathbf{#1}}\\int \\sin\\left(x\\right)\\,{d x}'
        """
        return cls(r'\newcommand{\Bold}[1]{\mathbf{#1}}'
                   r'\int \sin\left(x\right)\,{d x}')

    def print_to_stdout(self):
        r"""
        Write the data to stdout.

        This is just a convenience method to help with debugging.

        EXAMPLES::

            sage: from sage.repl.rich_output.output_catalog import OutputLatex
            sage: rich_output = OutputLatex.example()
            sage: rich_output.print_to_stdout()
            \newcommand{\Bold}[1]{\mathbf{#1}}\int \sin\left(x\right)\,{d x}
        """
        print(self.latex.get_str())
