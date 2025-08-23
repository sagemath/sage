# sage_setup: distribution = sagemath-environment
r"""
Features for testing the presence of ``latex`` and equivalent programs
"""

# ****************************************************************************
#       Copyright (C) 2021 Sebastien Labbe <slabqc@gmail.com>
#                     2021 Matthias Koeppe
#                     2022 Kwankyu Lee
#                     2022 Sebastian Oehms
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from . import StaticFile, Executable, FeatureTestResult, FeatureNotPresentError

latex_url = 'https://www.latex-project.org/'
latex_spkg = 'texlive'


class LaTeX(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of ``latex``.

    EXAMPLES::

        sage: from sage.features.latex import latex
        sage: latex().is_present()             # optional - latex
        FeatureTestResult('latex', True)
    """
    def __init__(self, name):
        r"""
        TESTS::

            sage: from sage.features.latex import latex
            sage: isinstance(latex(), latex)
            True
        """
        super().__init__(name, executable=name, spkg=latex_spkg, url=latex_url)

    def is_functional(self):
        r"""
        Return whether ``latex`` in the path is functional.

        EXAMPLES::

            sage: from sage.features.latex import latex
            sage: latex().is_functional()             # optional - latex
            FeatureTestResult('latex', True)

        When the feature is not functional, more information on the reason
        can be obtained as follows::

            sage: result = latex().is_functional()    # not tested
            sage: print(result.reason)                # not tested
            Running latex on a sample file
            (with command='latex -interaction=nonstopmode tmp_wmpos8ak.tex')
            returned nonzero exit status='1' with stderr=''
            and stdout='This is pdfTeX,
            ...
            Runaway argument?
            {document
            ! File ended while scanning use of \end.
            ...
            No pages of output.
            Transcript written on tmp_wmpos8ak.log.'
        """
        lines = []
        lines.append(r"\documentclass{article}")
        lines.append(r"\begin{document}")
        lines.append(r"$\alpha+2$")
        lines.append(r"\end{document}")
        content = '\n'.join(lines)

        # create a simple tex file with the content
        from sage.misc.temporary_file import tmp_filename
        base_filename_tex = tmp_filename(ext='.tex')
        with open(base_filename_tex, 'w') as f:
            f.write(content)
        import os
        base, filename_tex = os.path.split(base_filename_tex)

        # running latex
        from subprocess import run
        cmd = [self.name, '-interaction=nonstopmode', filename_tex]
        cmd = ' '.join(cmd)
        result = run(cmd, shell=True, cwd=base, capture_output=True, text=True)

        # return
        if result.returncode == 0:
            return FeatureTestResult(self, True)
        else:
            return FeatureTestResult(self, False, reason="Running latex on "
                                     "a sample file (with command='{}') returned nonzero "
                                     "exit status='{}' with stderr='{}' "
                                     "and stdout='{}'".format(result.args,
                                                              result.returncode,
                                                              result.stderr.strip(),
                                                              result.stdout.strip()))


class latex(LaTeX):
    r"""
    A :class:`~sage.features.Feature` describing the presence of ``latex``.

    EXAMPLES::

        sage: from sage.features.latex import latex
        sage: latex().is_present()             # optional - latex
        FeatureTestResult('latex', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.latex import latex
            sage: isinstance(latex(), latex)
            True
        """
        super().__init__("latex")


class pdflatex(LaTeX):
    r"""
    A :class:`~sage.features.Feature` describing the presence of ``pdflatex``.

    EXAMPLES::

        sage: from sage.features.latex import pdflatex
        sage: pdflatex().is_present()             # optional - pdflatex
        FeatureTestResult('pdflatex', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.latex import pdflatex
            sage: isinstance(pdflatex(), pdflatex)
            True
        """
        super().__init__("pdflatex")


class xelatex(LaTeX):
    r"""
    A :class:`~sage.features.Feature` describing the presence of ``xelatex``.

    EXAMPLES::

        sage: from sage.features.latex import xelatex
        sage: xelatex().is_present()             # optional - xelatex
        FeatureTestResult('xelatex', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.latex import xelatex
            sage: isinstance(xelatex(), xelatex)
            True
        """
        super().__init__("xelatex")


class lualatex(LaTeX):
    r"""
    A :class:`~sage.features.Feature` describing the presence of ``lualatex``.

    EXAMPLES::

        sage: from sage.features.latex import lualatex
        sage: lualatex().is_present()             # optional - lualatex
        FeatureTestResult('lualatex', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.latex import lualatex
            sage: isinstance(lualatex(), lualatex)
            True
        """
        super().__init__("lualatex")


class dvips(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of ``dvips``.

    EXAMPLES::

        sage: from sage.features.latex import dvips
        sage: dvips().is_present()             # optional - dvips
        FeatureTestResult('dvips', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.latex import dvips
            sage: isinstance(dvips(), dvips)
            True
        """
        Executable.__init__(self, 'dvips', executable='dvips',
                            url='https://tug.org/texinfohtml/dvips.html')


class TeXFile(StaticFile):
    r"""
    A :class:`sage.features.Feature` describing the presence of a TeX file.

    EXAMPLES::

        sage: from sage.features.latex import TeXFile
        sage: TeXFile('x', 'x.tex').is_present()  # optional - latex
        FeatureTestResult('x', True)
    """
    def __init__(self, name, filename, **kwds):
        r"""
        Initialize.

        TESTS::

            sage: from sage.features.latex import TeXFile
            sage: TeXFile('nonexisting', 'xxxxxx-nonexisting-file.tex').is_present()  # optional - latex
            FeatureTestResult('nonexisting', False)
        """
        StaticFile.__init__(self, name, filename, search_path=[], **kwds)

    def absolute_filename(self) -> str:
        r"""
        The absolute path of the file.

        EXAMPLES::

            sage: from sage.features.latex import TeXFile
            sage: feature = TeXFile('latex_class_article', 'article.cls')
            sage: feature.absolute_filename()  # optional - latex
            '.../latex/base/article.cls'
        """
        from subprocess import run, CalledProcessError
        try:
            proc = run(['kpsewhich', self.filename],
                       capture_output=True, text=True, check=True)
            return proc.stdout.strip()
        except CalledProcessError:
            reason = "{filename!r} not found by kpsewhich".format(filename=self.filename)
            raise FeatureNotPresentError(self, reason)

    def _is_present(self):
        r"""
        Test for the presence of the TeX file.

        EXAMPLES::

            sage: from sage.features.latex import LaTeXPackage, latex
            sage: f = LaTeXPackage("tkz-graph")
            sage: g = latex()
            sage: not f.is_present() or bool(g.is_present())  # indirect doctest
            True
        """
        return latex().is_present() and super()._is_present()


class LaTeXPackage(TeXFile):
    r"""
    A :class:`sage.features.Feature` describing the presence of a LaTeX package
    (``.sty`` file).

    EXAMPLES::

        sage: from sage.features.latex import LaTeXPackage
        sage: LaTeXPackage('graphics').is_present()  # optional - latex
        FeatureTestResult('latex_package_graphics', True)
    """
    @staticmethod
    def __classcall__(cls, package_name, **kwds):
        """
        TESTS::

            sage: from sage.features.latex import LaTeXPackage
            sage: LaTeXPackage('graphics') is LaTeXPackage('graphics')
            True
        """
        return TeXFile.__classcall__(cls,
                                     f'latex_package_{package_name}'.replace('-', '_'),
                                     f'{package_name}.sty',
                                     **kwds)


def all_features():
    return [latex(),
            pdflatex(),
            xelatex(),
            lualatex(),
            dvips(),
            LaTeXPackage("tkz-graph")]
