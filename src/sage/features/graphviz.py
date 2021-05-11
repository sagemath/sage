# -*- coding: utf-8 -*-
r"""
Check for graphviz
"""
# ****************************************************************************
#       Copyright (C) 2018 Sebastien Labbe <slabqc@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from . import Feature, Executable, FeatureTestResult


class dot(Executable):
    r"""
    A :class:`sage.features.Executable` describing the presence of
    ``dot``

    EXAMPLES::

        sage: from sage.features.graphviz import dot
        sage: dot().is_present()  # optional: graphviz
        FeatureTestResult('dot', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.graphviz import dot
            sage: isinstance(dot(), dot)
            True
        """
        Executable.__init__(self, "dot", executable="dot",
                            spkg="graphviz",
                            url="https://www.graphviz.org/")


class neato(Executable):
    r"""
    A :class:`sage.features.Executable` describing the presence of
    ``neato``

    EXAMPLES::

        sage: from sage.features.graphviz import neato
        sage: neato().is_present()  # optional: graphviz
        FeatureTestResult('neato', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.graphviz import neato
            sage: isinstance(neato(), neato)
            True
        """
        Executable.__init__(self, "neato", executable="neato",
                            spkg="graphviz",
                            url="https://www.graphviz.org/")


class twopi(Executable):
    r"""
    A :class:`sage.features.Executable` describing the presence of
    ``twopi``

    EXAMPLES::

        sage: from sage.features.graphviz import twopi
        sage: twopi().is_present()  # optional: graphviz
        FeatureTestResult('twopi', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.graphviz import twopi
            sage: isinstance(twopi(), twopi)
            True
        """
        Executable.__init__(self, "twopi", executable="twopi",
                            spkg="graphviz",
                            url="https://www.graphviz.org/")


class Graphviz(Feature):
    r"""
    A :class:`sage.features.Feature` describing the presence of
    ``dot``, ``neato`` and ``twopi``.

    EXAMPLES::

        sage: from sage.features.graphviz import Graphviz
        sage: Graphviz().is_present()  # optional: graphviz
        FeatureTestResult('Graphviz', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.graphviz import Graphviz
            sage: isinstance(Graphviz(), Graphviz)
            True
        """
        Feature.__init__(self, "Graphviz",
                         spkg="graphviz",
                         url="https://www.graphviz.org/")

    def _is_present(self):
        r"""
        EXAMPLES::

            sage: from sage.features.graphviz import Graphviz
            sage: Graphviz()._is_present() # optional: graphviz
            FeatureTestResult('Graphviz', True)
        """
        test = (dot()._is_present() and
                neato()._is_present() and
                twopi()._is_present())
        if not test:
            return test
        else:
            return FeatureTestResult(self, True)
