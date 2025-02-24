# sage_setup: distribution = sagemath-environment
r"""
Check for poppler features

poppler-utils is a collection of tools built on Poppler's library API, to
manage PDF and extract contents:

 - ``pdfattach`` -- add a new embedded file (attachment) to an existing PDF
 - ``pdfdetach`` -- extract embedded documents from a PDF
 - ``pdffonts`` -- lists the fonts used in a PDF
 - ``pdfimages`` -- extract all embedded images at native resolution from a PDF
 - ``pdfinfo`` -- list all information of a PDF
 - ``pdfseparate`` -- extract single pages from a PDF
 - ``pdftocairo`` -- convert single pages from a PDF to vector or bitmap formats using cairo
 - ``pdftohtml`` -- convert PDF to HTML format retaining formatting
 - ``pdftoppm`` -- convert a PDF page to a bitmap
 - ``pdftops`` -- convert PDF to printable PS format
 - ``pdftotext`` -- extract all text from PDF
 - ``pdfunite`` -- merges several PDF

Currently we only check for the presence of ``pdftocairo``.
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

from . import Executable

class pdftocairo(Executable):
    r"""
    A :class:`sage.features.Feature` describing the presence of
    ``pdftocairo``

    EXAMPLES::

        sage: from sage.features.poppler import pdftocairo
        sage: pdftocairo().is_present()             # optional: pdftocairo
        FeatureTestResult('pdftocairo', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.poppler import pdftocairo
            sage: isinstance(pdftocairo(), pdftocairo)
            True
        """
        Executable.__init__(self, "pdftocairo", executable='pdftocairo',
                            url='https://poppler.freedesktop.org/')

def all_features():
    return [pdftocairo()]
