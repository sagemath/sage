# -*- coding: utf-8 -*-
r"""
Check for pdflatex and equivalent programs
"""
# ****************************************************************************
#       Copyright (C) 2021 Sebastien Labbe <slabqc@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from . import Executable, FeatureTestResult

class latex(Executable):
    r"""
    A :class:`sage.features.Feature` describing the presence of ``latex``

    EXAMPLES::

        sage: from sage.features.latex import latex
        sage: latex().is_present()             # optional: latex
        FeatureTestResult('latex', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.latex import latex
            sage: isinstance(latex(), latex)
            True
        """
        Executable.__init__(self, "latex", executable="latex",
                            url="https://www.latex-project.org/")

    def is_functional(self):
        r"""
        Return whether `latex` in the path is functional.

        EXAMPLES:

            sage: from sage.features.latex import latex
            sage: latex().is_functional()             # optional: latex
            FeatureTestResult('latex', True)
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
        cmd = ['latex', '-interaction=nonstopmode', filename_tex]
        cmd = ' '.join(cmd)
        result = run(cmd, shell=True, cwd=base, capture_output=True, text=True)

        # return
        if result.returncode == 0:
            return FeatureTestResult(self, True)
        else:
            return FeatureTestResult(self, False, reason="Running latex on "
                        "a sample file returned non-zero "
                        "exit status {}".format(result.returncode))

class pdflatex(Executable):
    r"""
    A :class:`sage.features.Feature` describing the presence of ``pdflatex``

    EXAMPLES::

        sage: from sage.features.latex import pdflatex
        sage: pdflatex().is_present()             # optional: pdflatex
        FeatureTestResult('pdflatex', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.latex import pdflatex
            sage: isinstance(pdflatex(), pdflatex)
            True
        """
        Executable.__init__(self, "pdflatex", executable="pdflatex",
                            url="https://www.latex-project.org/")

class xelatex(Executable):
    r"""
    A :class:`sage.features.Feature` describing the presence of ``xelatex``

    EXAMPLES::

        sage: from sage.features.latex import xelatex
        sage: xelatex().is_present()             # optional: xelatex
        FeatureTestResult('xelatex', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.latex import xelatex
            sage: isinstance(xelatex(), xelatex)
            True
        """
        Executable.__init__(self, "xelatex", executable="xelatex",
                            url="https://www.latex-project.org/")

class lualatex(Executable):
    r"""
    A :class:`sage.features.Feature` describing the presence of ``lualatex``

    EXAMPLES::

        sage: from sage.features.latex import lualatex
        sage: lualatex().is_present()             # optional: lualatex
        FeatureTestResult('lualatex', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.latex import lualatex
            sage: isinstance(lualatex(), lualatex)
            True
        """
        Executable.__init__(self, "lualatex", executable="lualatex",
                            url="https://www.latex-project.org/")


def all_features():
    return [latex(),
            pdflatex(),
            xelatex(),
            lualatex()]
