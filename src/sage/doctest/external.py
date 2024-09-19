# sage_setup: distribution = sagemath-repl
"""
Detecting external software

This module makes up a list of external software that Sage interfaces. Availability
of each software is tested only when necessary. This is mainly used for the doctests
which require certain external software installed on the system.

Even though the functions in this module should also work when an external
software is not present, most doctests in this module are only tested if
testing of external software is explicitly enabled in order to avoid invoking
external software otherwise. See :issue:`28819` for details.

AUTHORS:

- Kwankyu Lee (2016-03-09) -- initial version, based on code by Robert Bradshaw and Nathann Cohen
"""

#*****************************************************************************
#       Copyright (C) 2016      Kwankyu Lee <ekwankyu@gmail.com>
#                     2018      Thierry Monteil
#                     2018-2021 Sébastien Labbé
#                     2019      Markus Wageringel
#                     2020      John H. Palmieri
#                     2021      Matthias Koeppe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import multiprocessing
import os

# With OS X, Python 3.8 defaults to use 'spawn' instead of 'fork' in
# multiprocessing, and Sage doctesting doesn't work with 'spawn'. See
# trac #27754.
if os.uname().sysname == 'Darwin':
    multiprocessing.set_start_method('fork', force=True)
Array = multiprocessing.Array

# Functions in this module whose name is of the form 'has_xxx' tests if the
# software xxx is available to Sage.
prefix = 'has_'

def has_internet():
    """
    Test if Internet is available.

    Failure of connecting to the site "https://www.sagemath.org" within a second
    is regarded as internet being not available.

    EXAMPLES::

        sage: from sage.doctest.external import has_internet
        sage: has_internet()  # random, optional -- internet
        FeatureTestResult('internet', True)
    """
    from sage.features.internet import Internet
    return Internet().is_present()

def has_latex():
    """
    Test if Latex is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_latex
        sage: has_latex() # optional - latex
        FeatureTestResult('latex', True)
    """
    from sage.features.latex import latex
    return latex().is_present()

def has_xelatex():
    """
    Test if xelatex is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_xelatex
        sage: has_xelatex()   # optional - xelatex
        FeatureTestResult('xelatex', True)
    """
    from sage.features.latex import xelatex
    return xelatex().is_present()

def has_pdflatex():
    """
    Test if pdflatex is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_pdflatex
        sage: has_pdflatex()   # optional - pdflatex
        FeatureTestResult('pdflatex', True)
    """
    from sage.features.latex import pdflatex
    return pdflatex().is_present()

def has_lualatex():
    """
    Test if lualatex is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_lualatex
        sage: has_lualatex()   # optional - lualatex
        FeatureTestResult('lualatex', True)
    """
    from sage.features.latex import lualatex
    return lualatex().is_present()

def has_magma():
    """
    Test if Magma is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_magma
        sage: has_magma() # random, optional - magma
        True
    """
    from sage.features.interfaces import Magma
    return Magma().is_present()

def has_matlab():
    """
    Test if Matlab is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_matlab
        sage: has_matlab() # random, optional - matlab
        True
    """
    from sage.features.interfaces import Matlab
    return Matlab().is_present()

def has_mathematica():
    """
    Test if Mathematica is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_mathematica
        sage: has_mathematica() # random, optional - mathematica
        True
    """
    from sage.features.interfaces import Mathematica
    return Mathematica().is_present()

def has_maple():
    """
    Test if Maple is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_maple
        sage: has_maple() # random, optional - maple
        True
    """
    from sage.features.interfaces import Maple
    return Maple().is_present()

def has_macaulay2():
    """
    Test if Macaulay2 is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_macaulay2
        sage: has_macaulay2() # random, optional - macaulay2
        True
    """
    from sage.features.interfaces import Macaulay2
    return Macaulay2().is_present()

def has_octave():
    """
    Test if Octave is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_octave
        sage: has_octave() # random, optional - octave
        True
    """
    from sage.features.interfaces import Octave
    return Octave().is_present()

def has_pandoc():
    """
    Test if pandoc is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_pandoc
        sage: has_pandoc()      # optional -- pandoc
        FeatureTestResult('pandoc', True)
    """
    from sage.features.pandoc import Pandoc
    return Pandoc().is_present()

def has_scilab():
    """
    Test if Scilab is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_scilab
        sage: has_scilab() # random, optional - scilab
        True
    """
    from sage.interfaces.scilab import scilab
    try:
        scilab('2+3')
        return True
    except Exception:
        return False

def has_cplex():
    """
    Test if CPLEX is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_cplex
        sage: has_cplex() # random, optional - CPLEX
        FeatureTestResult('cplex', True)
    """
    from sage.features.mip_backends import CPLEX
    return CPLEX().is_present()

def has_gurobi():
    """
    Test if Gurobi is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_gurobi
        sage: has_gurobi() # random, optional - Gurobi
        FeatureTestResult('gurobi', True)
    """
    from sage.features.mip_backends import Gurobi
    return Gurobi().is_present()

def has_graphviz():
    """
    Test if graphviz (dot, twopi, neato) are available.

    EXAMPLES::

        sage: from sage.doctest.external import has_graphviz
        sage: has_graphviz()   # optional -- graphviz
        FeatureTestResult('graphviz', True)
    """
    from sage.features.graphviz import Graphviz
    return Graphviz().is_present()

def has_ffmpeg():
    """
    Test if ffmpeg is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_ffmpeg
        sage: has_ffmpeg()      # optional -- ffmpeg
        FeatureTestResult('ffmpeg', True)
    """
    from sage.features.ffmpeg import FFmpeg
    return FFmpeg().is_present()

def has_imagemagick():
    """
    Test if ImageMagick (command magick or convert) is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_imagemagick
        sage: has_imagemagick() # optional -- imagemagick
        FeatureTestResult('imagemagick', True)
    """
    from sage.features.imagemagick import ImageMagick
    return ImageMagick().is_present()

def has_dvipng():
    """
    Test if dvipng is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_dvipng
        sage: has_dvipng() # optional -- dvipng
        FeatureTestResult('dvipng', True)
    """
    from sage.features.dvipng import dvipng
    return dvipng().is_present()

def has_pdf2svg():
    """
    Test if pdf2svg is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_pdf2svg
        sage: has_pdf2svg() # optional -- pdf2svg
        FeatureTestResult('pdf2svg', True)
    """
    from sage.features.pdf2svg import pdf2svg
    return pdf2svg().is_present()

def has_rubiks():
    """
    Test if the rubiks package (``cu2``, ``cubex``, ``dikcube``,
    ``mcube``, ``optimal``, and ``size222``) is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_rubiks
        sage: has_rubiks()   # optional -- rubiks
        FeatureTestResult('rubiks', True)
    """
    from sage.features.rubiks import Rubiks
    return Rubiks().is_present()

def has_4ti2():
    """
    Test if the 4ti2 package is available.

    EXAMPLES::

        sage: from sage.doctest.external import has_4ti2
        sage: has_4ti2()   # optional -- 4ti2
        FeatureTestResult('4ti2', True)
    """
    from sage.features.four_ti_2 import FourTi2
    return FourTi2().is_present()

def external_features():
    r"""
    Generate the features that are only to be tested if ``--optional=external`` is used.

    EXAMPLES::

        sage: from sage.doctest.external import external_features
        sage: next(external_features())
        Feature('internet')
    """
    from sage.features.internet import Internet
    yield Internet()
    import sage.features.latex
    yield from sage.features.latex.all_features()
    import sage.features.ffmpeg
    yield from sage.features.ffmpeg.all_features()
    import sage.features.interfaces
    for feature in sage.features.interfaces.all_features():
        if feature.name != 'mathics':
            yield feature
    from sage.features.mip_backends import CPLEX, Gurobi
    yield CPLEX()
    yield Gurobi()

def external_software() -> list[str]:
    """
    Return the alphabetical list of external software supported by this module.

    EXAMPLES::

        sage: from sage.doctest.external import external_software
        sage: sorted(external_software) == external_software
        True
    """
    return sorted(f.name for f in external_features())


external_software = external_software()


class AvailableSoftware():
    """
    This class keeps the set of available software whose availability is detected lazily
    from the list of external software.

    EXAMPLES::

        sage: from sage.doctest.external import external_software, available_software
        sage: external_software
        ['asir',
         'cplex',
         'dvips',
         'ffmpeg',
         'gurobi',
         'internet',
         'latex',
         'latex_package_tkz_graph',
         'lualatex',
         'macaulay2',
         'magma',
         'maple',
         'mathematica',
         'matlab',
         'octave',
         'pdflatex',
         'scilab',
         'xelatex']
        sage: 'internet' in available_software # random, optional - internet
        True
        sage: available_software.issuperset(set(['internet','latex'])) # random, optional - internet latex
        True
    """
    def __init__(self):
        """
        Initialization.

        EXAMPLES::

            sage: from sage.doctest.external import AvailableSoftware
            sage: S = AvailableSoftware()
            sage: S.seen() # random
            []
        """
        self._allow_external = True
        # For multiprocessing of doctests, the data self._seen should be
        # shared among subprocesses. Thus we use Array class from the
        # multiprocessing module.
        from sage.features.all import all_features
        self._external_features = set(external_features())
        features = set(self._external_features)
        features.update(all_features())
        self._features = sorted(features, key=lambda feature: feature.name)
        self._indices = {feature.name: idx for idx, feature in enumerate(self._features)}
        self._seen = Array('i', len(self._features)) # initialized to zeroes
        self._hidden = Array('i', len(self._features)) # initialized to zeroes

    def __contains__(self, item):
        """
        Return ``True`` if ``item`` is available on the system.

        EXAMPLES::

            sage: from sage.doctest.external import available_software
            sage: 'internet' in available_software # random, optional - internet
            True
        """
        try:
            idx = self._indices[item]
        except KeyError:
            return False
        feature = self._features[idx]
        if feature.is_hidden():
            if not self._hidden[idx]:
                self._hidden[idx] = 1
            available = False  # a hidden feature is considered to be not available
        else:
            if not self._allow_external and feature in self._external_features:
                # an external feature is considered to be not available
                # if this is not allowed
                available = False
            elif feature.is_present():
                available = True
            else:
                available = False
        if available:
            if not self._seen[idx]:
                self._seen[idx] = 1
            return True
        else:
            return False

    def issuperset(self, other):
        """
        Return ``True`` if ``other`` is a subset of ``self``.

        EXAMPLES::

            sage: from sage.doctest.external import available_software
            sage: available_software.issuperset(set(['internet','latex','magma'])) # random, optional - internet latex magma
            True
        """
        for item in other:
            if item not in self:
                return False
        return True

    def detectable(self):
        """
        Return the list of names of those features for which testing their presence is allowed.
        """
        return [feature.name
                for feature, seen in zip(self._features, self._seen)
                if seen >= 0 and (self._allow_external or feature not in self._external_features)]

    def seen(self):
        """
        Return the list of detected external software.

        EXAMPLES::

            sage: from sage.doctest.external import available_software
            sage: available_software.seen() # random
            ['internet', 'latex', 'magma']
        """
        return [feature.name
                for feature, seen in zip(self._features, self._seen)
                if seen > 0]

    def hidden(self):
        """
        Return the list of detected hidden external software.

        EXAMPLES::

            sage: # needs conway_polynomials database_cremona_mini_ellcurve database_ellcurves database_graphs
            sage: from sage.doctest.external import available_software
            sage: from sage.features.databases import all_features
            sage: for f in all_features():
            ....:    f.hide()
            ....:    if f._spkg_type() == 'standard':
            ....:         test = f.name in available_software
            ....:    f.unhide()
            sage: sorted(available_software.hidden())
            [...'conway_polynomials',...
             'database_cremona_mini_ellcurve',...
             'database_ellcurves',...
             'database_graphs'...]
        """
        return [feature.name
                for feature, hidden in zip(self._features, self._hidden)
                if hidden > 0]


available_software = AvailableSoftware()
