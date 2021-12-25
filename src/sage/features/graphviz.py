# -*- coding: utf-8 -*-
r"""
Features for testing the presence of ``graphviz``
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
from . import Executable
from .join_feature import JoinFeature


class dot(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of ``dot``

    EXAMPLES::

        sage: from sage.features.graphviz import dot
        sage: dot().is_present()  # optional - graphviz
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
    A :class:`~sage.features.Feature` describing the presence of ``neato``

    EXAMPLES::

        sage: from sage.features.graphviz import neato
        sage: neato().is_present()  # optional - graphviz
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
    A :class:`~sage.features.Feature` describing the presence of ``twopi``

    EXAMPLES::

        sage: from sage.features.graphviz import twopi
        sage: twopi().is_present()  # optional - graphviz
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


class Graphviz(JoinFeature):
    r"""
    A :class:`~sage.features.Feature` describing the presence of
    the ``dot``, ``neato``, and ``twopi`` executables from the
    ``graphviz`` package.

    EXAMPLES::

        sage: from sage.features.graphviz import Graphviz
        sage: Graphviz().is_present()  # optional - graphviz
        FeatureTestResult('graphviz', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.graphviz import Graphviz
            sage: isinstance(Graphviz(), Graphviz)
            True
        """
        JoinFeature.__init__(self, "graphviz",
                             [dot(), neato(), twopi()],
                             spkg="graphviz",
                             url="https://www.graphviz.org/")


def all_features():
    return [Graphviz()]
