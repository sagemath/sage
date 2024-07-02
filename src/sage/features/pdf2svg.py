# sage_setup: distribution = sagemath-environment
r"""
Feature for testing the presence of ``pdf2svg``
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

class pdf2svg(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of :ref:`pdf2svg <spkg_pdf2svg>`.

    EXAMPLES::

        sage: from sage.features.pdf2svg import pdf2svg
        sage: pdf2svg().is_present()             # optional - pdf2svg
        FeatureTestResult('pdf2svg', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.pdf2svg import pdf2svg
            sage: isinstance(pdf2svg(), pdf2svg)
            True
        """
        Executable.__init__(self, "pdf2svg", executable='pdf2svg',
                            spkg='pdf2svg',
                            url='http://www.cityinthesky.co.uk/opensource/pdf2svg/')


def all_features():
    return [pdf2svg()]
