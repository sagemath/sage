r"""
Standalone LaTeX Document class and TikzPicture

This module contains two Python classes. Firstly, it contains a class
:class:`Standalone` to represent a LaTeX file using the standalone__
document class.

__ http://www.ctan.org/pkg/standalone

From its documentation:

    *The standalone bundle allows users to easily place picture environments
    or other material in own source files and compile these on their own or as
    part of a main document. A special standalone class is provided for use
    with such files, which by default crops the resulting output file to the
    content. The standalone package enables the user to simply load the
    standalone files using ``\input`` inside a main document.*

Secondly, it contains a class :class:`TikzPicture` which inherits from
:class:`Standalone` that represents a LaTeX file using the standalone
document class and containing a tikzpicture.

A Python Module for PGF/Tikz pictures. A TikzPicture object is created from
a string starting with ``r'\begin{tikzpicture}'`` and ending with
``r'\end{tikzpicture}'``.

The module allows to convert a standalone LaTeX document class file,
including tikzpictures, to an image. It allows conversion to pdf, png and
svg formats. It also show them automatically in Jupyter using rich
representation.

According to wikipedia, `PGF/TikZ`__ is a pair of languages for producing
vector graphics (e.g., technical illustrations and drawings) from a
geometric/algebraic description, with standard features including the
drawing of points, lines, arrows, paths, circles, ellipses and polygons.

__ https://www.ctan.org/pkg/pgf

EXAMPLES:

Standalone LaTeX document class
-------------------------------

First *Hello World* example::

    sage: from sage.misc.latex_standalone import Standalone
    sage: Standalone('Hello World')
    \documentclass{standalone}
    \begin{document}
    Hello World
    \end{document}

Loading a few latex packages::

    sage: Standalone('Hello World', usepackage=['amsmath', 'amsfont'])
    \documentclass{standalone}
    \usepackage{amsmath}
    \usepackage{amsfont}
    \begin{document}
    Hello World
    \end{document}

Setting few standalone options (see documentation of standalone for a
complete list)::

    sage: Standalone('Hello World', standalone_config=["border=4mm", "beamer=true"])
    \documentclass{standalone}
    \standaloneconfig{border=4mm}
    \standaloneconfig{beamer=true}
    \begin{document}
    Hello World
    \end{document}

Adding your own list of macros::

    sage: Standalone('Hello World', macros=[r'\newcommand{\ZZ}{\mathbb{Z}}'])
    \documentclass{standalone}
    \newcommand{\ZZ}{\mathbb{Z}}
    \begin{document}
    Hello World
    \end{document}

It provides conversion to images of different format::

    sage: from sage.misc.latex_standalone import Standalone
    sage: s = Standalone('Hello World')
    sage: _ = s.pdf()    # not tested
    sage: _ = s.png()    # not tested
    sage: _ = s.svg()    # not tested
    sage: s              # not tested, in Jupyter, this shows the image directly below the cell

TikzPicture
-----------

This module also contains a class :class:`TikzPicture` which inherits from
:class:`Standalone` to represent more specifically a tikzpicture which is
within a standalone document class.

First construct a string describing a tikzpicture::

    sage: lines = []
    sage: lines.append(r'\begin{tikzpicture}')
    sage: lines.append(r'\draw[very thick,orange,->] (0,0) -- (1,1);')
    sage: lines.append(r'\end{tikzpicture}')
    sage: s = '\n'.join(lines)
    sage: print(s)
    \begin{tikzpicture}
    \draw[very thick,orange,->] (0,0) -- (1,1);
    \end{tikzpicture}

One may provide it as input to ``TikzPicture``::

    sage: from sage.misc.latex_standalone import TikzPicture
    sage: t = TikzPicture(s)

In the terminal, the following shows the content of the standalone
document class tex file which contains the tikzpicture. In Jupyter, it
shows the picture itself::

    sage: t
    \documentclass[tikz]{standalone}
    \begin{document}
    \begin{tikzpicture}
    \draw[very thick,orange,->] (0,0) -- (1,1);
    \end{tikzpicture}
    \end{document}

As it is the case for :class:`Standalone`, the constructor of
``TikzPicture`` has many arguments allowing for example to add some more
``\usepackage`` lines::

    sage: t = TikzPicture(s, usepackage=['amsmath'])
    sage: t
    \documentclass[tikz]{standalone}
    \usepackage{amsmath}
    \begin{document}
    \begin{tikzpicture}
    \draw[very thick,orange,->] (0,0) -- (1,1);
    \end{tikzpicture}
    \end{document}

Moreover, it allows to load some tikz libraries::

    sage: t = TikzPicture(s, usetikzlibrary=['arrows'])
    sage: t
    \documentclass[tikz]{standalone}
    \usetikzlibrary{arrows}
    \begin{document}
    \begin{tikzpicture}
    \draw[very thick,orange,->] (0,0) -- (1,1);
    \end{tikzpicture}
    \end{document}

The following example illustrates that it works when providing the
tikzpicture code generated by Sage from some polyhedron::

    sage: from sage.misc.latex_standalone import TikzPicture
    sage: V = [[1,0,1], [1,0,0], [1,1,0], [0,0,-1],
    ....:      [0,1,0], [-1,0,0], [0,1,1], [0,0,1], [0,-1,0]]
    sage: P = Polyhedron(vertices=V).polar()                                            # needs sage.geometry.polyhedron
    sage: s1 = P.projection().tikz([674,108,-731],112, output_type='LatexExpr')         # needs sage.geometry.polyhedron sage.plot
    sage: t1 = TikzPicture(s1)                                                          # needs sage.geometry.polyhedron sage.plot

Open the image in a viewer (the returned value is a string giving the
absolute path to the file in some temporary directory)::

    sage: path_to_file = t1.pdf()               # not tested                            # needs sage.geometry.polyhedron sage.plot

Instead, you may save a pdf of the tikzpicture into a file of your choice
(but this does not open the viewer)::

    sage: _ = t1.pdf('tikz_polytope.pdf')       # not tested                            # needs sage.geometry.polyhedron sage.plot

Opening the image in a viewer can be turned off::

    sage: _ = t1.pdf(view=False)                # long time (2s), optional - latex, needs sage.geometry.polyhedron sage.plot

The same can be done with png format (translated from pdf with convert
command which needs the installation of imagemagick)::

    sage: _ = t1.png(view=False)                # long time (2s), optional - imagemagick latex, needs sage.geometry.polyhedron sage.plot

The string representation gives the header (5 lines) and tail (5 lines) of
the tikzpicture. In Jupyter, it will instead use rich representation and
show the image directly below the cell in png or svg format::

    sage: t1                                                                            # needs sage.geometry.polyhedron sage.plot
    \documentclass[tikz]{standalone}
    \begin{document}
    \begin{tikzpicture}%
            [x={(0.249656cm, -0.577639cm)},
            y={(0.777700cm, -0.358578cm)},
            z={(-0.576936cm, -0.733318cm)},
    ...
    \node[vertex] at (0.00000, -1.00000, 0.00000)     {};
    \node[vertex] at (-0.50000, -0.50000, -0.50000)     {};
    %%
    %%
    \end{tikzpicture}
    \end{document}

Use ``print(t)`` to see the complete content of the file::

    sage: print(t1)                             # not tested                            # needs sage.geometry.polyhedron sage.plot

Adding a border in the options avoids cropping the vertices of a graph::

    sage: # needs sage.graphs
    sage: g = graphs.PetersenGraph()
    sage: s2 = latex(g)   # takes 3s but the result is cached           # optional - latex
    sage: t2 = TikzPicture(s2, standalone_config=["border=4mm"],        # optional - latex
    ....:                  usepackage=['tkz-graph'])
    sage: _ = t2.pdf()                          # not tested

The current latex representation of a transducer is a tikzpicture using
the tikz library automata. The string can be used as input::

    sage: # needs sage.graphs sage.modules
    sage: s3 = latex(transducers.GrayCode())
    sage: t3 = TikzPicture(s3, usetikzlibrary=['automata'])
    sage: _ = t3.pdf(view=False)                # long time (2s)        # optional - latex

AUTHORS:

- Sébastien Labbé, initial version in slabbe-0.2.spkg, nov 2015.
- Sébastien Labbé, inclusion into SageMath from slabbe-0.6.2, July 2021.
"""

# ****************************************************************************
#       Copyright (C) 2015-2022 Sébastien Labbé <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from subprocess import run
import os

from sage.structure.sage_object import SageObject
from sage.misc.superseded import experimental


class Standalone(SageObject):
    r"""
    LaTeX standalone document class.

    INPUT:

    - ``content`` -- string; the content to be added in the document
      between lines ``r'\begin{document}'`` and ``r'\end{document}'``
    - ``document_class_options`` -- list of strings (default: ``[]``);
      latex document class standalone options. Such options appear on the
      line ``\documentclass[...]{standalone}`` between the brackets.
    - ``standalone_config`` -- list of strings (default: ``[]``);
      standalone configuration options. Such options are defined with
      ``\standaloneconfig{...}``.
    - ``usepackage`` -- list of strings (default: ``[]``); latex packages
    - ``macros`` -- list of strings (default: ``[]``); stuff you need for the picture
    - ``use_sage_preamble`` -- boolean (default: ``False``); whether to include sage
      latex preamble and sage latex macros, that is, the content of
      :func:`sage.misc.latex.extra_preamble()`,
      :func:`sage.misc.latex.extra_macros()` and
      :func:`sage.misc.latex_macros.sage_latex_macros()`

    EXAMPLES::

        sage: from sage.misc.latex_standalone import Standalone
        sage: content = "\\section{Intro}\nTest\n"
        sage: t = Standalone(content)
        sage: t
        \documentclass{standalone}
        \begin{document}
        \section{Intro}
        Test
        \end{document}

    ::

        sage: t = Standalone(content, standalone_config=["border=4mm"],
        ....:                usepackage=['amsmath'])
        sage: t
        \documentclass{standalone}
        \standaloneconfig{border=4mm}
        \usepackage{amsmath}
        \begin{document}
        \section{Intro}
        Test
        \end{document}
    """
    def __init__(self, content, document_class_options=None,
                 standalone_config=None, usepackage=None, macros=None,
                 use_sage_preamble=False):
        r"""
        See :class:`Standalone` for full information.

        EXAMPLES::

            sage: from sage.misc.latex_standalone import Standalone
            sage: content = "\\section{Intro}\n\nTest\n"
            sage: t = Standalone(content)
        """
        self._content = content
        self._document_class_options = [] if document_class_options is None else list(document_class_options)
        self._standalone_config = [] if standalone_config is None else standalone_config
        self._usepackage = [] if usepackage is None else usepackage
        self._macros = [] if macros is None else macros
        if use_sage_preamble:
            from sage.misc.latex import _Latex_prefs
            for key in ['preamble', 'macros']:
                s = _Latex_prefs._option[key]
                if s:
                    self._macros.append(s)
            from sage.misc.latex_macros import sage_latex_macros
            self._macros.extend(sage_latex_macros())

    def _latex_file_header_lines(self):
        r"""
        EXAMPLES::

            sage: latex.extra_preamble('')
            sage: from sage.misc.latex_standalone import Standalone
            sage: s = "\\begin{tikzpicture}\n\\draw (0,0) -- (1,1);\n\\end{tikzpicture}"
            sage: A = ['tikz']
            sage: B = ["border=4mm"]
            sage: C = ['amsmath']
            sage: t = Standalone(s, document_class_options=A, standalone_config=B, usepackage=C)
            sage: t._latex_file_header_lines()[:6]
            ['\\documentclass[tikz]{standalone}',
             '\\standaloneconfig{border=4mm}',
             '\\usepackage{amsmath}']
        """
        lines = []
        if self._document_class_options:
            options = ','.join(self._document_class_options)
            lines.append(r"\documentclass[{}]{{standalone}}".format(options))
        else:
            lines.append(r"\documentclass{standalone}")
        lines.extend(r"\standaloneconfig{{{}}}".format(config)
                     for config in self._standalone_config)
        lines.extend(r"\usepackage{{{}}}".format(package)
                     for package in self._usepackage)
        lines.extend(self._macros)
        return lines

    def _repr_(self):
        r"""
        Return a string representation of the Standalone file.

        It contains the first few and last few lines of the content.

        .. NOTE::

            Use ``print(t)`` or ``str(t)`` to show or get the full content.

        EXAMPLES:

        When the content has 10 lines or less, it shows it all::

            sage: from sage.misc.latex_standalone import Standalone
            sage: s = "\\begin{tikzpicture}\n\\draw (0,0) -- (1,1);\n\\end{tikzpicture}"
            sage: t = Standalone(s, document_class_options=['tikz'], standalone_config=["border=4mm"], usepackage=['amsmath'])
            sage: t
            \documentclass[tikz]{standalone}
            \standaloneconfig{border=4mm}
            \usepackage{amsmath}
            \begin{document}
            \begin{tikzpicture}
            \draw (0,0) -- (1,1);
            \end{tikzpicture}
            \end{document}

        When the content more than 10 lines, it shows the head (first 5
        lines) and tail (last 5 lines) of the content together with the
        complete header of the standalone latex document. The number of
        missing lines and the number of characters in the content is
        written allowing to detect a change if needed::

            sage: lines = []
            sage: lines.append(r'\begin{tikzpicture}')
            sage: lines.append(r'\draw[->] (-.5,0) -- (20,0);')
            sage: lines.extend(r'\draw({i},-.5) -- ({i},.5);'.format(i=i) for i in range(20))
            sage: lines.append(r'\end{tikzpicture}')
            sage: t = Standalone('\n'.join(lines), document_class_options=['tikz'])
            sage: t
            \documentclass[tikz]{standalone}
            \begin{document}
            \begin{tikzpicture}
            \draw[->] (-.5,0) -- (20,0);
            \draw(0,-.5) -- (0,.5);
            \draw(1,-.5) -- (1,.5);
            \draw(2,-.5) -- (2,.5);
            ---
            13 lines not printed (566 characters in total).
            Use print to see the full content.
            ---
            \draw(16,-.5) -- (16,.5);
            \draw(17,-.5) -- (17,.5);
            \draw(18,-.5) -- (18,.5);
            \draw(19,-.5) -- (19,.5);
            \end{tikzpicture}
            \end{document}
        """
        lines = self._latex_file_header_lines()
        lines.append(r"\begin{document}")
        L = self._content.splitlines()
        if len(L) <= 10:
            lines.extend(L)
        else:
            lines.extend(L[:5])
            lines.append('---')
            lines.append('{} lines not printed ({} characters in total).'.format(len(L) - 10,
                                                           len(self._content)))
            lines.append('Use print to see the full content.')
            lines.append('---')
            lines.extend(L[-5:])
        lines.append(r"\end{document}")
        return '\n'.join(lines)

    def _rich_repr_(self, display_manager, **kwds):
        r"""
        Rich Output Magic Method.

        See :mod:`sage.repl.rich_output` for details.

        EXAMPLES::

            sage: from sage.repl.rich_output import get_display_manager
            sage: dm = get_display_manager()
            sage: dm.is_in_terminal()
            False

            sage: from sage.misc.latex_standalone import TikzPicture
            sage: lines = []
            sage: lines.append(r'\begin{tikzpicture}')
            sage: lines.append(r'\draw[very thick,orange,->] (0,0) -- (1,1);')
            sage: lines.append(r'\end{tikzpicture}')
            sage: s = '\n'.join(lines)
            sage: t = TikzPicture(s)
            sage: t._rich_repr_(dm)      # random result is Text in doctest
            OutputImagePng container

        Using vector svg instead of png::

            sage: dm.preferences.graphics = 'vector'
            sage: t._rich_repr_(dm)      # random result is Text in doctest
            OutputImageSvg container
            sage: dm.preferences.graphics = 'raster'
        """
        # Do not use rich output in the terminal
        if display_manager.is_in_terminal():
            return
        # Do not use rich output if not in IPython notebook (Jupyter)
        from sage.repl.rich_output.backend_ipython import BackendIPythonNotebook
        if not isinstance(display_manager._backend, BackendIPythonNotebook):
            return

        types = display_manager.types
        prefer_raster = (
            ('png', types.OutputImagePng),
        )
        prefer_vector = (
            ('svg', types.OutputImageSvg),
            ('pdf', types.OutputImagePdf),
        )
        graphics = display_manager.preferences.graphics
        if graphics == 'disable':
            return
        elif graphics == 'raster' or graphics is None:
            preferred = prefer_raster + prefer_vector
        elif graphics == 'vector':
            preferred = prefer_vector + prefer_raster
        else:
            raise ValueError('unknown graphics output preference')

        for format, output_container in preferred:
            if output_container in display_manager.supported_output():
                filename = getattr(self, format)(view=False, **kwds)
                from sage.repl.rich_output.buffer import OutputBuffer
                buf = OutputBuffer.from_file(filename)
                return output_container(buf)

    def __str__(self):
        r"""
        Return the complete string of the standalone document class file.

        EXAMPLES::

            sage: from sage.misc.latex_standalone import Standalone
            sage: s = "\\begin{tikzpicture}\n\\draw (0,0) -- (1,1);\n\\end{tikzpicture}"
            sage: t = Standalone(s, document_class_options=['tikz'])
            sage: print(t)
            \RequirePackage{luatex85}
            \documentclass[tikz]{standalone}
            \begin{document}
            \begin{tikzpicture}
            \draw (0,0) -- (1,1);
            \end{tikzpicture}
            \end{document}
        """
        lines = []
        # LuaLaTeX, TeXLive 2016, standalone: undefined control sequence
        # https://tex.stackexchange.com/questions/315025
        # fixed in 2018, meanwhile, we add the fix here
        lines.append(r"\RequirePackage{luatex85}")
        lines.extend(self._latex_file_header_lines())
        lines.append(r"\begin{document}")
        lines.append(self._content)
        lines.append(r"\end{document}")
        return '\n'.join(lines)

    def content(self):
        r"""
        Return the content of the standalone document class file.

        EXAMPLES::

            sage: from sage.misc.latex_standalone import Standalone
            sage: t = Standalone('Hello World')
            sage: t.content()
            'Hello World'

        ::

            sage: from sage.misc.latex_standalone import TikzPicture
            sage: s = "\\begin{tikzpicture}\n\\draw (0,0) -- (1,1);\n\\end{tikzpicture}"
            sage: t = TikzPicture(s)
            sage: print(t.content())
            \begin{tikzpicture}
            \draw (0,0) -- (1,1);
            \end{tikzpicture}
        """
        return self._content

    def add_document_class_option(self, option):
        r"""
        Add a document class option.

        INPUT:

        - ``option`` -- string

        EXAMPLES::

            sage: from sage.misc.latex_standalone import Standalone
            sage: t = Standalone('Hello World')
            sage: t.add_document_class_option('beamer')
            sage: t
            \documentclass[beamer]{standalone}
            \begin{document}
            Hello World
            \end{document}
        """
        self._document_class_options.append(option)

    def add_standalone_config(self, config):
        r"""
        Add a standalone config.

        INPUT:

        - ``config`` -- string

        EXAMPLES::

            sage: from sage.misc.latex_standalone import Standalone
            sage: t = Standalone('Hello World')
            sage: t.add_standalone_config("border=4mm")
            sage: t
            \documentclass{standalone}
            \standaloneconfig{border=4mm}
            \begin{document}
            Hello World
            \end{document}
        """
        self._standalone_config.append(config)

    def add_usepackage(self, package):
        r"""
        Add a ``usepackage`` line.

        INPUT:

        - ``package`` -- string, name of package

        EXAMPLES::

            sage: from sage.misc.latex_standalone import Standalone
            sage: t = Standalone('Hello World')
            sage: t.add_usepackage('amsmath')
            sage: t
            \documentclass{standalone}
            \usepackage{amsmath}
            \begin{document}
            Hello World
            \end{document}
        """
        self._usepackage.append(package)

    def add_macro(self, macro):
        r"""
        Add a macro.

        INPUT:

        - ``macro`` -- string, newcommand line

        EXAMPLES::

            sage: from sage.misc.latex_standalone import Standalone
            sage: t = Standalone('Hello World')
            sage: t.add_macro(r'\newcommand{\ZZ}{\mathbb{Z}}')
            sage: t
            \documentclass{standalone}
            \newcommand{\ZZ}{\mathbb{Z}}
            \begin{document}
            Hello World
            \end{document}
        """
        self._macros.append(macro)

    def pdf(self, filename=None, view=True, program=None):
        r"""
        Compile the latex code with pdflatex and create a pdf file.

        INPUT:

        - ``filename`` -- string (default: ``None``); the output filename.
          If ``None``, it saves the file in a temporary directory.

        - ``view`` -- boolean (default: ``True``); whether to open the file in a
          pdf viewer. This option is ignored and automatically set to
          ``False`` if ``filename`` is not ``None``.

        - ``program`` -- string (default: ``None``); ``'pdflatex'`` or
          ``'lualatex'``. If ``None``, it uses ``'lualatex'`` if it is
          available, otherwise ``'pdflatex'``.

        OUTPUT: string, path to pdf file

        EXAMPLES::

            sage: from sage.misc.latex_standalone import Standalone
            sage: t = Standalone('Hello World')
            sage: _ = t.pdf(view=False)                 # long time (1s)    # optional - latex

        Same for instances of :class:`TikzPicture`::

            sage: from sage.misc.latex_standalone import TikzPicture
            sage: s = "\\begin{tikzpicture}\n\\draw (0,0) -- (1,1);\n\\end{tikzpicture}"
            sage: t = TikzPicture(s)
            sage: _ = t.pdf(view=False)                 # not tested

        A filename may be provided where to save the file, in which case
        the viewer does not open the file::

            sage: from sage.misc.temporary_file import tmp_filename
            sage: filename = tmp_filename('temp','.pdf')
            sage: path_to_file = t.pdf(filename)        # long time (1s)    # optional - latex
            sage: path_to_file[-4:]                     # long time (fast)  # optional - latex
            '.pdf'

        The filename may contain spaces::

            sage: filename = tmp_filename('filename with spaces','.pdf')
            sage: path_to_file = t.pdf(filename)        # long time (1s)    # optional - latex

        TESTS:

        We test the behavior when a wrong tex string is provided::

            sage: s = "\\begin{tikzpicture}\n\\draw (0,0) -- (1,1);\n\\end{tikzpicture}"
            sage: s_missing_last_character = s[:-1]
            sage: t = TikzPicture(s_missing_last_character)
            sage: _ = t.pdf()                                               # optional - latex
            Traceback (most recent call last):
            ...
            CalledProcessError: Command '['...latex', '-interaction=nonstopmode',
            'tikz_...tex']' returned nonzero exit status 1.
        """
        from sage.features.latex import lualatex, pdflatex

        # Set default program
        if program is None:
            if lualatex().is_present():
                program = 'lualatex'
            else:
                program = 'pdflatex'

        # Check availability of programs
        if program == 'pdflatex':
            pdflatex().require()
        elif program == 'lualatex':
            lualatex().require()
        else:
            raise ValueError("program(={}) should be pdflatex or lualatex".format(program))

        # set up filenames
        from sage.misc.temporary_file import tmp_filename
        temp_filename_tex = tmp_filename('tikz_', '.tex')
        with open(temp_filename_tex, 'w') as f:
            f.write(str(self))
        base, temp_filename_tex = os.path.split(temp_filename_tex)
        temp_filename, ext = os.path.splitext(temp_filename_tex)

        # running pdflatex or lualatex
        cmd = [program, '-interaction=nonstopmode', temp_filename_tex]
        result = run(cmd, cwd=base, capture_output=True, text=True)

        # If a problem with the tex source occurs, provide the log
        if result.returncode != 0:
            print("Command \n"
                  "   '{}'\n"
                  "returned nonzero exit status {}.\n"
                  "Here is the content of the stderr:{}\n"
                  "Here is the content of the stdout:"
                  "{}\n".format(' '.join(result.args),
                                result.returncode,
                                result.stderr.strip(),
                                result.stdout.strip()))
        result.check_returncode()
        temp_filename_pdf = os.path.join(base, temp_filename + '.pdf')

        # move the pdf into the good location
        if filename:
            filename = os.path.abspath(filename)
            import shutil
            shutil.move(temp_filename_pdf, filename)
            return filename

        # open the tmp pdf
        elif view:
            from sage.misc.viewer import pdf_viewer
            cmd = pdf_viewer().split()
            cmd.append(temp_filename_pdf)
            # we use check_call as opposed to run, because
            # it gives the sage prompt back to the user
            # see https://stackoverflow.com/a/71342967
            # run(cmd, cwd=base, capture_output=True, check=True)
            from subprocess import check_call, PIPE
            check_call(cmd, cwd=base, stdout=PIPE, stderr=PIPE)

        return temp_filename_pdf

    def dvi(self, filename=None, view=True, program='latex'):
        r"""
        Compile the latex code with latex and create a dvi file.

        INPUT:

        - ``filename`` -- string (default: ``None``); the output filename.
          If ``None``, it saves the file in a temporary directory

        - ``view`` -- boolean (default: ``True``); whether to open the file in a
          dvi viewer. This option is ignored and automatically set to
          ``False`` if ``filename`` is not ``None``.

        - ``program`` -- string (default: ``'latex'``); ``'latex'``

        OUTPUT: string, path to dvi file

        EXAMPLES::

            sage: from sage.misc.latex_standalone import Standalone
            sage: t = Standalone('Hello World')
            sage: _ = t.dvi(view=False)                 # long time (1s)    # optional - latex

        Same for instances of :class:`TikzPicture`::

            sage: from sage.misc.latex_standalone import TikzPicture
            sage: s = "\\begin{tikzpicture}\n\\draw (0,0) -- (1,1);\n\\end{tikzpicture}"
            sage: t = TikzPicture(s)
            sage: _ = t.dvi(view=False)                 # not tested

        A filename may be provided where to save the file, in which case
        the viewer does not open the file::

            sage: from sage.misc.temporary_file import tmp_filename
            sage: filename = tmp_filename('temp','.dvi')
            sage: path_to_file = t.dvi(filename)        # long time (1s)    # optional - latex
            sage: path_to_file[-4:]                     # long time (fast)  # optional - latex
            '.dvi'

        The filename may contain spaces::

            sage: filename = tmp_filename('filename with spaces','.dvi')
            sage: path_to_file = t.dvi(filename)        # long time (1s)    # optional - latex

        TESTS:

        We test the behavior when a wrong tex string is provided::

            sage: s = "\\begin{tikzpicture}\n\\draw (0,0) -- (1,1);\n\\end{tikzpicture}"
            sage: s_missing_last_character = s[:-1]
            sage: t = TikzPicture(s_missing_last_character)
            sage: _ = t.dvi()                                               # optional - latex
            Traceback (most recent call last):
            ...
            CalledProcessError: Command '['latex', '-interaction=nonstopmode',
            'tikz_...tex']' returned nonzero exit status 1.

        We test the behavior when a wrong value is provided::

            sage: t = Standalone('Hello World')
            sage: _ = t.dvi(program='lates')
            Traceback (most recent call last):
            ...
            ValueError: program(=lates) should be latex
        """
        from sage.features.latex import latex

        # Set default program
        if program is None:
            program = 'latex'

        # Check availability of programs
        if program == 'latex':
            latex().require()
        else:
            raise ValueError("program(={}) should be latex".format(program))

        # set up filenames
        from sage.misc.temporary_file import tmp_filename
        temp_filename_tex = tmp_filename('tikz_', '.tex')
        with open(temp_filename_tex, 'w') as f:
            f.write(str(self))
        base, temp_filename_tex = os.path.split(temp_filename_tex)
        temp_filename, ext = os.path.splitext(temp_filename_tex)

        # running pdflatex or lualatex
        cmd = [program, '-interaction=nonstopmode', temp_filename_tex]
        result = run(cmd, cwd=base, capture_output=True, text=True)

        # If a problem with the tex source occurs, provide the log
        if result.returncode != 0:
            print("Command \n"
                  "   '{}'\n"
                  "returned nonzero exit status {}.\n"
                  "Here is the content of the stderr:{}\n"
                  "Here is the content of the stdout:"
                  "{}\n".format(' '.join(result.args),
                                result.returncode,
                                result.stderr.strip(),
                                result.stdout.strip()))
        result.check_returncode()
        temp_filename_dvi = os.path.join(base, temp_filename + '.dvi')

        # move the pdf into the good location
        if filename:
            filename = os.path.abspath(filename)
            import shutil
            shutil.move(temp_filename_dvi, filename)
            return filename

        # open the tmp dvi
        elif view:
            from sage.misc.viewer import dvi_viewer
            cmd = dvi_viewer().split()
            cmd.append(temp_filename_dvi)
            # we use check_call as opposed to run, because
            # it gives the sage prompt back to the user
            # see https://stackoverflow.com/a/71342967
            # run(cmd, cwd=base, capture_output=True, check=True)
            from subprocess import check_call, PIPE
            check_call(cmd, cwd=base, stdout=PIPE, stderr=PIPE)

        return temp_filename_dvi

    def png(self, filename=None, density=150, view=True):
        r"""
        Compile the latex code with pdflatex and converts to a png file.

        INPUT:

        - ``filename`` -- string (default: ``None``); the output filename.
          If ``None``, it saves the file in a temporary directory.

        - ``density`` -- integer (default: ``150``); horizontal and vertical
          density of the image

        - ``view`` -- boolean (default: ``True``); whether to open the file in a
          png viewer. This option is ignored and automatically set to
          ``False`` if ``filename`` is not ``None``.

        OUTPUT: string, path to png file

        EXAMPLES::

            sage: from sage.misc.latex_standalone import Standalone
            sage: t = Standalone('Hello World')
            sage: _ = t.png(view=False)             # long time (1s)    # optional - latex imagemagick

        Same for instances of :class:`TikzPicture`::

            sage: from sage.misc.latex_standalone import TikzPicture
            sage: s = "\\begin{tikzpicture}\n\\draw (0,0) -- (1,1);\n\\end{tikzpicture}"
            sage: t = TikzPicture(s)
            sage: _ = t.png(view=False)             # not tested

        ::

            sage: from sage.misc.temporary_file import tmp_filename
            sage: filename = tmp_filename('temp','.png')
            sage: path_to_file = t.png(filename)    # long time (1s)    # optional - latex imagemagick
            sage: path_to_file[-4:]                 # long time (fast)  # optional - latex imagemagick
            '.png'
        """
        from sage.features.imagemagick import ImageMagick
        ImageMagick().require()

        temp_filename_pdf = self.pdf(filename=None, view=False)
        temp_filename, ext = os.path.splitext(temp_filename_pdf)
        temp_filename_png = temp_filename + '.png'

        # convert to png
        cmd = ['convert', '-density',
               '{0}x{0}'.format(density), '-trim', temp_filename_pdf,
               temp_filename_png]
        result = run(cmd, capture_output=True, text=True)

        # If a problem occurs, provide the log
        if result.returncode != 0:
            print("Command \n"
                  "   '{}'\n"
                  "returned nonzero exit status {}.\n"
                  "Here is the content of the stderr:{}\n"
                  "Here is the content of the stdout:"
                  "{}\n".format(' '.join(result.args),
                                result.returncode,
                                result.stderr.strip(),
                                result.stdout.strip()))
        result.check_returncode()

        # move the png into the good location
        if filename:
            filename = os.path.abspath(filename)
            import shutil
            shutil.move(temp_filename_png, filename)
            return filename

        # open the tmp png
        elif view:
            from sage.misc.viewer import png_viewer
            cmd = png_viewer().split()
            cmd.append(temp_filename_png)
            # we use check_call as opposed to run, because
            # it gives the sage prompt back to the user
            # see https://stackoverflow.com/a/71342967
            # run(cmd, capture_output=True, check=True)
            from subprocess import check_call, PIPE
            check_call(cmd, stdout=PIPE, stderr=PIPE)

        return temp_filename_png

    def svg(self, filename=None, view=True, program='pdftocairo'):
        r"""
        Compile the latex code with pdflatex and converts to a svg file.

        INPUT:

        - ``filename`` -- string (default: ``None``); the output filename.
          If ``None``, it saves the file in a temporary directory.

        - ``view`` -- boolean (default: ``True``); whether to open the file in
          a browser. This option is ignored and automatically set to
          ``False`` if ``filename`` is not ``None``.

        - ``program`` -- string (default: ``'pdftocairo'``); ``'pdftocairo'`` or
          ``'pdf2svg'``

        OUTPUT: string, path to svg file

        EXAMPLES::

            sage: from sage.misc.latex_standalone import Standalone
            sage: t = Standalone('Hello World')
            sage: _ = t.svg(view=False)     # not tested

        Same for instances of :class:`TikzPicture`::

            sage: from sage.misc.latex_standalone import TikzPicture
            sage: s = "\\begin{tikzpicture}\n\\draw (0,0) -- (1,1);\n\\end{tikzpicture}"
            sage: t = TikzPicture(s)
            sage: _ = t.svg(view=False)     # not tested

        ::

            sage: from sage.misc.temporary_file import tmp_filename
            sage: filename = tmp_filename('temp', '.svg')
            sage: path_to_file = t.svg(filename,            # long time (1s)    # optional - latex pdf2svg
            ....:                      program='pdf2svg')
            sage: path_to_file[-4:]                         # long time (fast)  # optional - latex pdf2svg
            '.svg'
            sage: path_to_file = t.svg(filename,            # long time (1s)    # optional - latex pdftocairo
            ....:                      program='pdftocairo')
            sage: path_to_file[-4:]                         # long time (fast)  # optional - latex pdftocairo
            '.svg'
        """
        # set the temporary filenames
        temp_filename_pdf = self.pdf(filename=None, view=False)
        temp_filename, ext = os.path.splitext(temp_filename_pdf)
        temp_filename_svg = temp_filename + '.svg'

        # set the command
        if program == 'pdftocairo':
            from sage.features.poppler import pdftocairo
            pdftocairo().require()
            cmd = ['pdftocairo', '-svg', temp_filename_pdf, temp_filename_svg]
        elif program == 'pdf2svg':
            from sage.features.pdf2svg import pdf2svg
            pdf2svg().require()
            cmd = ['pdf2svg', temp_filename_pdf, temp_filename_svg]
        else:
            raise ValueError("program(={}) should be 'pdftocairo' or"
                             " 'pdf2svg'".format(program))

        # convert to svg
        result = run(cmd, capture_output=True, text=True)

        # If a problem occurs, provide the log
        if result.returncode != 0:
            print("Command \n"
                  "   '{}'\n"
                  "returned nonzero exit status {}.\n"
                  "Here is the content of the stderr:{}\n"
                  "Here is the content of the stdout:"
                  "{}\n".format(' '.join(result.args),
                                result.returncode,
                                result.stderr.strip(),
                                result.stdout.strip()))
        result.check_returncode()

        # move the svg into the good location
        if filename:
            filename = os.path.abspath(filename)
            import shutil
            shutil.move(temp_filename_svg, filename)
            return filename

        # open the tmp svg
        elif view:
            from sage.misc.viewer import browser
            cmd = browser().split()
            cmd.append(temp_filename_svg)
            # we use check_call as opposed to run, because
            # it gives the sage prompt back to the user
            # see https://stackoverflow.com/a/71342967
            # run(cmd, capture_output=True, check=True)
            from subprocess import check_call, PIPE
            check_call(cmd, stdout=PIPE, stderr=PIPE)

        return temp_filename_svg

    def eps(self, filename=None, view=True, program='dvips'):
        r"""
        Compile the latex code with pdflatex and converts to a eps file.

        INPUT:

        - ``filename`` -- string (default: ``None``); the output filename.
          If ``None``, it saves the file in a temporary directory

        - ``view`` -- boolean (default: ``True``); whether to open the file in
          a browser. This option is ignored and automatically set to
          ``False`` if ``filename`` is not ``None``.

        - ``program`` -- string (default: ``'dvips'``);
          ``'pdftocairo'`` or ``'dvips'``

        OUTPUT: string, path to eps file

        EXAMPLES::

            sage: from sage.misc.latex_standalone import Standalone
            sage: t = Standalone('Hello World')
            sage: _ = t.eps(view=False)     # not tested

        Same for instances of :class:`TikzPicture`::

            sage: from sage.misc.latex_standalone import TikzPicture
            sage: s = "\\begin{tikzpicture}\n\\draw (0,0) -- (1,1);\n\\end{tikzpicture}"
            sage: t = TikzPicture(s)
            sage: _ = t.eps(view=False)     # not tested

        We test the creation of the files::

            sage: from sage.misc.temporary_file import tmp_filename
            sage: filename = tmp_filename('temp', '.eps')
            sage: path_to_file = t.eps(filename,            # long time (1s)    # optional - latex dvips
            ....:                      program='dvips')
            sage: path_to_file[-4:]                         # long time (fast)  # optional - latex dvips
            '.eps'
            sage: path_to_file = t.eps(filename,            # long time (1s)    # optional - latex pdftocairo
            ....:                      program='pdftocairo')
            sage: path_to_file[-4:]                         # long time (fast)  # optional - latex pdftocairo
            '.eps'

        TESTS:

        We test the behavior when a wrong value is provided::

            sage: t = Standalone('Hello World')
            sage: _ = t.eps(program='convert')
            Traceback (most recent call last):
            ...
            ValueError: program(=convert) should be 'pdftocairo' or 'dvips'
        """

        if program == 'pdftocairo':
            from sage.features.poppler import pdftocairo
            pdftocairo().require()
            # set the temporary filenames
            temp_filename_pdf = self.pdf(filename=None, view=False)
            temp_filename, ext = os.path.splitext(temp_filename_pdf)
            temp_filename_eps = temp_filename + '.eps'
            # set the command
            cmd = ['pdftocairo', '-eps', temp_filename_pdf, temp_filename_eps]
        elif program == 'dvips':
            from sage.features.latex import dvips
            dvips().require()
            # set the temporary filenames
            temp_filename_dvi = self.dvi(filename=None, view=False)
            temp_filename, ext = os.path.splitext(temp_filename_dvi)
            temp_filename_eps = temp_filename + '.eps'
            # set the command
            cmd = ['dvips', '-E', '-o', temp_filename_eps, temp_filename_dvi]
        else:
            raise ValueError("program(={}) should be 'pdftocairo' or"
                             " 'dvips'".format(program))

        # convert to eps
        result = run(cmd, capture_output=True, text=True)

        # If a problem occurs, provide the log
        if result.returncode != 0:
            print("Command \n"
                  "   '{}'\n"
                  "returned nonzero exit status {}.\n"
                  "Here is the content of the stderr:{}\n"
                  "Here is the content of the stdout:"
                  "{}\n".format(' '.join(result.args),
                                result.returncode,
                                result.stderr.strip(),
                                result.stdout.strip()))
        result.check_returncode()

        # move the eps into the good location
        if filename:
            filename = os.path.abspath(filename)
            import shutil
            shutil.move(temp_filename_eps, filename)
            return filename

        # open the tmp eps
        elif view:
            from sage.misc.viewer import viewer
            cmd = viewer().split()
            cmd.append(temp_filename_eps)
            # we use check_call as opposed to run, because
            # it gives the sage prompt back to the user
            # see https://stackoverflow.com/a/71342967
            # run(cmd, capture_output=True, check=True)
            from subprocess import check_call, PIPE
            check_call(cmd, stdout=PIPE, stderr=PIPE)

        return temp_filename_eps

    def tex(self, filename=None, content_only=False, include_header=None):
        r"""
        Writes the latex code to a file.

        INPUT:

        - ``filename`` -- string (default: ``None``); the output filename.
          If ``None``, it saves the file in a temporary directory.
        - ``content_only`` -- boolean (default: ``False``); whether to include
          the header latex part. If ``True``, it prints only the
          content to the file.

        OUTPUT: string, path to tex file

        EXAMPLES::

            sage: from sage.misc.latex_standalone import Standalone
            sage: t = Standalone('Hello World')
            sage: _ = t.tex()
            sage: _ = t.tex(content_only=True)

        Same for instances of :class:`TikzPicture`::

            sage: from sage.misc.latex_standalone import TikzPicture
            sage: s = "\\begin{tikzpicture}\n\\draw (0,0) -- (1,1);\n\\end{tikzpicture}"
            sage: t = TikzPicture(s)
            sage: _ = t.tex()
            sage: _ = t.tex(content_only=True)

        Write to a given filename::

            sage: from sage.misc.temporary_file import tmp_filename
            sage: filename = tmp_filename('temp','.tex')
            sage: path_to_file = t.tex(filename)
            sage: path_to_file[-4:]
            '.tex'
        """
        if filename is None:
            from sage.misc.temporary_file import tmp_filename
            filename = tmp_filename('tikz_', '.tex')
        else:
            filename = os.path.abspath(filename)

        if include_header is not None:
            content_only = not include_header
            from sage.misc.superseded import deprecation
            deprecation(20343, "When merging this code from slabbe into "
                    "SageMath the argument include_header=False was "
                    "replaced by content_only=True. Please update your code "
                    "before include_header option gets removed from SageMath.")

        if content_only:
            output = self.content()
        else:
            output = str(self)

        with open(filename, 'w') as f:
            f.write(output)

        return filename

    def save(self, filename, **kwds):
        r"""
        Save the graphics to an image file.

        INPUT:

        - ``filename`` -- string. The filename and the image format
          given by the extension, which can be one of the following:

            * ``.pdf``,
            * ``.png``,
            * ``.svg``,
            * ``.eps``,
            * ``.dvi``,
            * ``.sobj`` (for a Sage object you can load later),
            * empty extension will be treated as ``.sobj``.

        All other keyword arguments will be passed to the plotter.

        OUTPUT: none

        .. NOTE::

            This method follows the signature of the method
            :meth:`sage.plot.Graphics.save` in order to be compatible with
            with sagetex. In particular so that ``\sageplot{t}`` written
            in a ``tex`` file works when ``t`` is an instance of
            :class:`Standalone` or :class:`TikzPicture`.

        EXAMPLES::

            sage: from sage.misc.temporary_file import tmp_filename
            sage: from sage.misc.latex_standalone import Standalone
            sage: t = Standalone('Hello World')
            sage: filename = tmp_filename('temp','.pdf')
            sage: t.save(filename)                          # long time (1s)    # optional - latex
            sage: filename = tmp_filename('temp','.eps')
            sage: t.save(filename)                          # long time (1s)    # optional - latex dvips
        """
        ext = os.path.splitext(filename)[1].lower()
        if ext == '' or ext == '.sobj':
            raise NotImplementedError()
        elif ext == '.pdf':
            self.pdf(filename, **kwds)
        elif ext == '.png':
            self.png(filename, **kwds)
        elif ext == '.svg':
            self.svg(filename, **kwds)
        elif ext == '.eps':
            self.eps(filename, **kwds)
        elif ext == '.dvi':
            self.dvi(filename, **kwds)
        else:
            raise ValueError("allowed file extensions for images are "
                             ".pdf, .png, .svg, .eps, .dvi!")


class TikzPicture(Standalone):
    r"""
    A TikzPicture embedded in a LaTeX standalone document class.

    INPUT:

    - ``content`` -- string, tikzpicture code starting with ``r'\begin{tikzpicture}'``
      and ending with ``r'\end{tikzpicture}'``
    - ``standalone_config`` -- list of strings (default: ``[]``);
      latex document class standalone configuration options
    - ``usepackage`` -- list of strings (default: ``[]``); latex
      packages
    - ``usetikzlibrary`` -- list of strings (default: ``[]``); tikz libraries
      to use
    - ``macros`` -- list of strings (default: ``[]``); stuff you need for the picture
    - ``use_sage_preamble`` -- boolean (default: ``False``); whether to include sage
      latex preamble and sage latex macros, that is, the content of
      :func:`sage.misc.latex.extra_preamble()`,
      :func:`sage.misc.latex.extra_macros()` and
      :func:`sage.misc.latex_macros.sage_latex_macros()`

    EXAMPLES:

    Create your own tikz string from scratch and provide it::

        sage: from sage.misc.latex_standalone import TikzPicture
        sage: lines = []
        sage: lines.append(r'\begin{tikzpicture}')
        sage: lines.append(r'\draw[very thick,orange,->] (0,0) -- (1,1);')
        sage: lines.append(r'\end{tikzpicture}')
        sage: s = '\n'.join(lines)
        sage: t = TikzPicture(s)
        sage: t
        \documentclass[tikz]{standalone}
        \begin{document}
        \begin{tikzpicture}
        \draw[very thick,orange,->] (0,0) -- (1,1);
        \end{tikzpicture}
        \end{document}

    Then use it by exporting the tikzpicture to other formats, all of the
    below methods return a string providing the path to the filename, which
    is by default in a temporary folder::

        sage: # not tested
        sage: _ = t.pdf()
        sage: _ = t.png()
        sage: _ = t.svg()
        sage: _ = t.tex()
        sage: _ = t.pdf(filename='abc.pdf')

    Here we create a tikzpicture for the latex representation of a graph.
    This is using tkz-graph tex library::

        sage: # needs sage.graphs
        sage: g = graphs.PetersenGraph()
        sage: s = latex(g)                                              # optional - latex
        sage: t = TikzPicture(s, standalone_config=["border=4mm"],      # optional - latex
        ....:                 usepackage=['tkz-graph'])
        sage: _ = t.pdf(view=False)                     # long time (2s), optional - latex latex_package_tkz_graph

    Here are standalone configurations, packages, tikz libraries and macros
    that can be set::

        sage: options = ['preview', 'border=4mm', 'beamer', 'float']
        sage: usepackage = ['nicefrac', 'amsmath', 'pifont', 'tikz-3dplot',
        ....:    'pgfplots']
        sage: tikzlib = ['arrows', 'snakes', 'backgrounds', 'patterns',
        ....:      'matrix', 'shapes', 'fit', 'calc', 'shadows', 'plotmarks',
        ....:      'positioning', 'pgfplots.groupplots', 'mindmap']
        sage: macros = [r'\newcommand{\ZZ}{\mathbb{Z}}']
        sage: s = "\\begin{tikzpicture}\n\\draw (0,0) -- (1,1);\n\\end{tikzpicture}"
        sage: t = TikzPicture(s, standalone_config=options, usepackage=usepackage,
        ....:        usetikzlibrary=tikzlib, macros=macros)
        sage: _ = t.pdf(view=False)                     # long time (2s), optional - latex
    """
    def __init__(self, content, standalone_config=None, usepackage=None,
            usetikzlibrary=None, macros=None, use_sage_preamble=False):
        r"""
        See :class:`TikzPicture` for full information.

        EXAMPLES::

            sage: from sage.misc.latex_standalone import TikzPicture
            sage: s = "\\begin{tikzpicture}\n\\draw (0,0) -- (1,1);\n\\end{tikzpicture}"
            sage: t = TikzPicture(s)
        """
        Standalone.__init__(self, content, document_class_options=['tikz'],
            standalone_config=standalone_config, usepackage=usepackage,
            macros=macros, use_sage_preamble=use_sage_preamble)

        self._usetikzlibrary = [] if usetikzlibrary is None else usetikzlibrary

    def _latex_file_header_lines(self):
        r"""
        EXAMPLES::

            sage: latex.extra_preamble('')
            sage: from sage.misc.latex_standalone import TikzPicture
            sage: s = "\\begin{tikzpicture}\n\\draw (0,0) -- (1,1);\n\\end{tikzpicture}"
            sage: t = TikzPicture(s, standalone_config=["border=4mm"], usepackage=['tkz-graph'])
            sage: t._latex_file_header_lines()[:6]
            ['\\documentclass[tikz]{standalone}',
             '\\standaloneconfig{border=4mm}',
             '\\usepackage{tkz-graph}']
        """
        lines = Standalone._latex_file_header_lines(self)
        for library in self._usetikzlibrary:
            lines.append(r"\usetikzlibrary{{{}}}".format(library))
        return lines

    def add_usetikzlibrary(self, library):
        r"""
        Add a ``usetikzlibrary`` line.

        INPUT:

        - ``library`` -- string, name of library

        EXAMPLES::

            sage: from sage.misc.latex_standalone import TikzPicture
            sage: s = "\\begin{tikzpicture}\n\\draw (0,0) -- (1,1);\n\\end{tikzpicture}"
            sage: t = TikzPicture(s)
            sage: t.add_usetikzlibrary('arrows')
            sage: t
            \documentclass[tikz]{standalone}
            \usetikzlibrary{arrows}
            \begin{document}
            \begin{tikzpicture}
            \draw (0,0) -- (1,1);
            \end{tikzpicture}
            \end{document}
        """
        self._usetikzlibrary.append(library)

    @classmethod
    def from_dot_string(cls, dotdata, prog='dot'):
        r"""
        Convert a graph to a tikzpicture using graphviz and dot2tex.

        .. NOTE::

            Prerequisite: dot2tex optional Sage package and graphviz must be
            installed.

        INPUT:

        - ``dotdata`` -- dot format string
        - ``prog`` -- string (default: ``'dot'``); the program used for the
          layout corresponding to one of the software of the graphviz
          suite: ``'dot'``, ``'neato'``, ``'twopi'``, ``'circo'`` or ``'fdp'``

        EXAMPLES::

            sage: # needs sage.graphs
            sage: from sage.misc.latex_standalone import TikzPicture
            sage: G = graphs.PetersenGraph()
            sage: dotdata = G.graphviz_string()
            sage: tikz = TikzPicture.from_dot_string(dotdata)   # long time (3s), optional - dot2tex graphviz
            sage: _ = tikz.pdf()      # not tested

        ::

            sage: # needs sage.graphs
            sage: dotdata = G.graphviz_string(labels='latex')
            sage: tikz = TikzPicture.from_dot_string(dotdata)   # long time (3s), optional - dot2tex graphviz
            sage: _ = tikz.pdf()      # not tested

        ::

            sage: # needs sage.combinat sage.graphs sage.groups
            sage: W = CoxeterGroup(["A",2])
            sage: G = W.cayley_graph()
            sage: dotdata = G.graphviz_string()
            sage: tikz = TikzPicture.from_dot_string(dotdata)   # long time (3s), optional - dot2tex graphviz
            sage: _ = tikz.pdf()      # not tested

        ::

            sage: # needs sage.combinat sage.graphs sage.groups
            sage: dotdata = G.graphviz_string(labels='latex')
            sage: tikz = TikzPicture.from_dot_string(dotdata)   # long time (3s), optional - dot2tex graphviz
            sage: _ = tikz.pdf()      # not tested
        """
        from sage.features import PythonModule
        PythonModule("dot2tex").require()
        from sage.features.graphviz import Graphviz
        Graphviz().require()

        import dot2tex
        tikz = dot2tex.dot2tex(dotdata,
                               format='tikz',
                               autosize=True,
                               crop=True,
                               figonly='True',
                               prog=prog).strip()
        return TikzPicture(tikz, standalone_config=["border=4mm"],
                           usetikzlibrary=['shapes'])

    @classmethod
    @experimental(issue_number=20343)
    def from_graph(cls, graph, merge_multiedges=True,
            merge_label_function=tuple, **kwds):
        r"""
        Convert a graph to a tikzpicture using graphviz and dot2tex.

        .. NOTE::

            Prerequisite: dot2tex optional Sage package and graphviz must be
            installed.

        .. WARNING::

            This method might be deleted in the future in favor of a method
            in the graph class returning a tikz picture.

        INPUT:

        - ``graph`` -- graph
        - ``merge_multiedges`` -- boolean (default: ``True``); if the graph
          has multiple edges, whether to merge the multiedges into one
          single edge
        - ``merge_label_function`` -- function (default: ``tuple``); a
          function to apply to each list of labels to be merged. It is
          ignored if ``merge_multiedges`` is not ``True`` or if the graph
          has no multiple edges.

        Other inputs are used for latex drawing with dot2tex and graphviz:

        - ``prog`` -- string (default: ``'dot'``); the program used for the
          layout corresponding to one of the software of the graphviz
          suite: ``'dot'``, ``'neato'``, ``'twopi'``, ``'circo'`` or ``'fdp'``
        - ``edge_labels`` -- boolean (default: ``True``)
        - ``color_by_label`` -- boolean (default: ``False``)
        - ``rankdir`` -- string (default: ``'down'``)
        - ``subgraph_clusters`` -- (default: ``[]``) a list of lists of
          vertices, if supported by the layout engine, nodes belonging to
          the same cluster subgraph are drawn together, with the entire
          drawing of the cluster contained within a bounding rectangle.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: from sage.misc.latex_standalone import TikzPicture
            sage: g = graphs.PetersenGraph()
            sage: tikz = TikzPicture.from_graph(g)                              # optional - dot2tex graphviz
            doctest:...: FutureWarning: This class/method/function is marked as experimental.
            It, its functionality or its interface might change without a formal deprecation.
            See https://github.com/sagemath/sage/issues/20343 for details.
            sage: _ = tikz.pdf()      # not tested

        Using ``prog``::

            sage: # needs sage.graphs
            sage: tikz = TikzPicture.from_graph(g, prog='neato',        # long time (3s), optional - dot2tex graphviz
            ....:                               color_by_label=True)
            sage: _ = tikz.pdf()      # not tested

        Using ``rankdir``::

            sage: tikz = TikzPicture.from_graph(g, rankdir='right')     # long time (3s), optional - dot2tex graphviz, needs sage.graphs
            sage: _ = tikz.pdf()      # not tested

        Using ``merge_multiedges``::

            sage: # needs sage.graphs sage.modules sage.symbolic
            sage: alpha = var('alpha')
            sage: m = matrix(2, range(4)); m.set_immutable()
            sage: G = DiGraph([(0,1,alpha), (0,1,0), (0,2,9), (0,2,m)],
            ....:             multiedges=True)
            sage: tikz = TikzPicture.from_graph(G, merge_multiedges=True)       # optional - dot2tex graphviz
            sage: _ = tikz.pdf()      # not tested

        Using ``merge_multiedges`` with ``merge_label_function``::

            sage: # needs sage.graphs
            sage: fn = lambda L: LatexExpr(','.join(map(str, L)))
            sage: edges = [(0,1,'a'), (0,1,'b'), (0,2,'c'), (0,2,'d')]
            sage: G = DiGraph(edges, multiedges=True)
            sage: tikz = TikzPicture.from_graph(G,                              # optional - dot2tex graphviz
            ....:           merge_multiedges=True, merge_label_function=fn)
            sage: _ = tikz.pdf()      # not tested

        Using subgraphs clusters (broken when using labels, see
        :issue:`22070`)::

            sage: S = FiniteSetMaps(5)
            sage: I = S((0,1,2,3,4))
            sage: a = S((0,1,3,0,0))
            sage: b = S((0,2,4,1,0))
            sage: roots = [I]
            sage: succ = lambda v: [v*a,v*b,a*v,b*v]
            sage: R = RecursivelyEnumeratedSet(roots, succ)
            sage: G = R.to_digraph()                                                    # needs sage.graphs
            sage: G                                                                     # needs sage.graphs
            Looped multi-digraph on 27 vertices
            sage: C = G.strongly_connected_components()                                 # needs sage.graphs
            sage: tikz = TikzPicture.from_graph(G,                              # optional - dot2tex graphviz, needs sage.graphs
            ....:              merge_multiedges=False, subgraph_clusters=C)
            sage: _ = tikz.pdf()      # not tested

        An example coming from ``graphviz_string`` documentation in SageMath::

            sage: # needs sage.graphs sage.symbolic
            sage: f(x) = -1 / x
            sage: g(x) = 1 / (x + 1)
            sage: G = DiGraph()
            sage: G.add_edges((i, f(i), f) for i in (1, 2, 1/2, 1/4))
            sage: G.add_edges((i, g(i), g) for i in (1, 2, 1/2, 1/4))
            sage: tikz = TikzPicture.from_graph(G)      # optional - dot2tex graphviz
            sage: _ = tikz.pdf()      # not tested
            sage: def edge_options(data):
            ....:     u, v, label = data
            ....:     options = {"color": {f: "red", g: "blue"}[label]}
            ....:     if (u,v) == (1/2, -2): options["label"]       = "coucou"; options["label_style"] = "string"
            ....:     if (u,v) == (1/2,2/3): options["dot"]         = "x=1,y=2"
            ....:     if (u,v) == (1,   -1): options["label_style"] = "latex"
            ....:     if (u,v) == (1,  1/2): options["dir"]         = "back"
            ....:     return options
            sage: tikz = TikzPicture.from_graph(G, edge_options=edge_options)   # optional - dot2tex graphviz
            sage: _ = tikz.pdf()      # not tested
        """
        from sage.features.latex import pdflatex
        pdflatex().require()
        from sage.features.graphviz import Graphviz
        Graphviz().require()
        from sage.features import PythonModule
        PythonModule("dot2tex").require()

        if merge_multiedges and graph.has_multiple_edges():
            from collections import defaultdict
            d = defaultdict(list)
            for (u, v, label) in graph.edges(sort=False):
                d[(u, v)].append(label)
            edges = [(u, v, merge_label_function(label_list)) for (u, v), label_list in d.items()]
            loops = graph.has_loops()
            if graph.is_directed():
                from sage.graphs.digraph import DiGraph
                graph = DiGraph(edges, format='list_of_edges', loops=loops)
            else:
                from sage.graphs.graph import Graph
                graph = Graph(edges, format='list_of_edges', loops=loops)

        options = {'format': 'dot2tex', 'edge_labels': True,
                   'color_by_label': False, 'prog': 'dot', 'rankdir': 'down'}
        options.update(kwds)

        graph.latex_options().set_options(**options)
        tikz = graph._latex_()
        return TikzPicture(tikz, standalone_config=["border=4mm"])

    @classmethod
    @experimental(issue_number=20343)
    def from_graph_with_pos(cls, graph, scale=1, merge_multiedges=True,
                            merge_label_function=tuple):
        r"""
        Convert a graph with positions defined for vertices to a tikzpicture.

        .. WARNING::

            This method might be deleted in the future in favor of a method
            in the graph class returning a tikz picture.

        INPUT:

        - ``graph`` -- graph (with predefined positions)
        - ``scale`` -- number (default: ``1``); tikzpicture scale
        - ``merge_multiedges`` -- boolean (default: ``True``); if the graph
          has multiple edges, whether to merge the multiedges into one
          single edge
        - ``merge_label_function`` -- function (default: ``tuple``); a
          function to apply to each list of labels to be merged. It is
          ignored if ``merge_multiedges`` is not ``True`` or if the graph
          has no multiple edges.

        EXAMPLES::

            sage: from sage.misc.latex_standalone import TikzPicture

            sage: # needs sage.graphs
            sage: g = graphs.PetersenGraph()
            sage: tikz = TikzPicture.from_graph_with_pos(g)
            doctest:...: FutureWarning: This class/method/function is marked as experimental.
            It, its functionality or its interface might change without a formal deprecation.
            See https://github.com/sagemath/sage/issues/20343 for details.

        ::

            sage: # needs sage.graphs
            sage: edges = [(0,0,'a'),(0,1,'b'),(0,1,'c')]
            sage: kwds = dict(format='list_of_edges', loops=True, multiedges=True)
            sage: G = DiGraph(edges, **kwds)
            sage: G.set_pos({0:(0,0), 1:(1,0)})
            sage: f = lambda label:','.join(label)
            sage: TikzPicture.from_graph_with_pos(G, merge_label_function=f)
            \documentclass[tikz]{standalone}
            \standaloneconfig{border=4mm}
            \begin{document}
            \begin{tikzpicture}
            [auto,scale=1]
            % vertices
            \node (node_0) at (0, 0) {0};
            \node (node_1) at (1, 0) {1};
            % edges
            \draw[->] (node_0) -- node {b,c} (node_1);
            % loops
            \draw (node_0) edge [loop above] node {a} ();
            \end{tikzpicture}
            \end{document}

        TESTS::

            sage: # needs sage.graphs
            sage: edges = [(0,0,'a'),(0,1,'b'),(0,1,'c')]
            sage: kwds = dict(format='list_of_edges', loops=True, multiedges=True)
            sage: G = DiGraph(edges, **kwds)
            sage: TikzPicture.from_graph_with_pos(G)
            Traceback (most recent call last):
            ...
            ValueError: vertex positions need to be set first
        """
        pos = graph.get_pos()
        if pos is None:
            raise ValueError('vertex positions need to be set first')

        if merge_multiedges and graph.has_multiple_edges():
            from collections import defaultdict
            d = defaultdict(list)
            for (u, v, label) in graph.edges(sort=True):
                d[(u, v)].append(label)
            edges = [(u, v, merge_label_function(label_list)) for (u, v), label_list in d.items()]
            loops = graph.has_loops()
            if graph.is_directed():
                from sage.graphs.digraph import DiGraph
                graph = DiGraph(edges, format='list_of_edges', loops=loops)
            else:
                from sage.graphs.graph import Graph
                graph = Graph(edges, format='list_of_edges', loops=loops)

        keys_for_vertices = graph._keys_for_vertices()

        lines = []
        lines.append(r'\begin{tikzpicture}')
        lines.append(r'[auto,scale={}]'.format(scale))

        # vertices
        lines.append(r'% vertices')
        for u in graph.vertices(sort=False):
            line = r'\node ({}) at {} {{{}}};'.format(keys_for_vertices(u),
                                                      pos[u], u)
            lines.append(line)

        # edges
        lines.append(r'% edges')
        arrow = '->' if graph.is_directed() else ''
        for (u, v, label) in graph.edges(sort=True):
            if u == v:
                # loops are done below
                continue
            if label:
                line = r'\draw[{}] ({}) -- node {{{}}} ({});'.format(arrow,
                                                    keys_for_vertices(u),
                                                    label,
                                                    keys_for_vertices(v))
            else:
                line = r'\draw[{}] ({}) -- ({});'.format(arrow,
                                                    keys_for_vertices(u),
                                                    keys_for_vertices(v))
            lines.append(line)

        # loops
        lines.append(r'% loops')
        for (u, v, label) in graph.loop_edges():
            line = r'\draw ({}) edge [loop above] node {{{}}} ();'.format(
                keys_for_vertices(u), label)
            lines.append(line)

        lines.append(r'\end{tikzpicture}')
        tikz = '\n'.join(lines)
        return TikzPicture(tikz, standalone_config=["border=4mm"])

    @classmethod
    @experimental(issue_number=20343)
    def from_poset(cls, poset, **kwds):
        r"""
        Convert a poset to a tikzpicture using graphviz and dot2tex.

        .. NOTE::

            Prerequisite: dot2tex optional Sage package and graphviz must be
            installed.

        .. WARNING::

            This method might be deleted in the future in favor of a method
            in the graph class returning a tikz picture.

        INPUT:

        - ``poset`` -- poset
        - ``prog`` -- string (default: ``'dot'``); the program used for the
          layout corresponding to one of the software of the graphviz
          suite: ``'dot'``, ``'neato'``, ``'twopi'``, ``'circo'`` or ``'fdp'``
        - ``edge_labels`` -- boolean (default: ``True``)
        - ``color_by_label`` -- boolean (default: ``False``)
        - ``rankdir`` -- string (default: ``'down'``)

        EXAMPLES::

            sage: from sage.misc.latex_standalone import TikzPicture

            sage: # needs sage.graphs sage.modules
            sage: P = posets.PentagonPoset()
            sage: tikz = TikzPicture.from_poset(P)                              # optional - dot2tex graphviz
            doctest:...: FutureWarning: This class/method/function is marked as experimental.
            It, its functionality or its interface might change without a formal deprecation.
            See https://github.com/sagemath/sage/issues/20343 for details.

        ::

            sage: tikz = TikzPicture.from_poset(P, prog='neato',        # long time (3s), optional - dot2tex, needs sage.graphs sage.modules
            ....:                               color_by_label=True)

        ::

            sage: # needs sage.graphs
            sage: P = posets.SymmetricGroupWeakOrderPoset(4)
            sage: tikz = TikzPicture.from_poset(P)                      # long time (4s), optional - dot2tex graphviz
            sage: tikz = TikzPicture.from_poset(P, prog='neato')        # long time (4s), optional - dot2tex graphviz
        """
        graph = poset.hasse_diagram()
        return cls.from_graph(graph, **kwds)

    def tikz_picture_code(self):
        r"""
        EXAMPLES::

            sage: from sage.misc.latex_standalone import TikzPicture
            sage: s = "\\begin{tikzpicture}\n\\draw (0,0) -- (1,1);\n\\end{tikzpicture}"
            sage: t = TikzPicture(s)
            sage: print(t.tikz_picture_code())
            \begin{tikzpicture}
            \draw (0,0) -- (1,1);
            \end{tikzpicture}
        """
        return self.content()
