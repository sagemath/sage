r"""
Feature for testing the presence of ``graph_genus``.
"""

# ****************************************************************************
#       Copyright (C) 2026 Alexander Metzger
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from . import PythonModule


class GraphGenus(PythonModule):
    r"""
    A :class:`~sage.features.Feature` describing the presence of
    :ref:`graph_genus <spkg_graph_genus>`.

    EXAMPLES::

        sage: from sage.features.graph_genus import GraphGenus
        sage: GraphGenus().is_present()                  # optional - graph_genus
        FeatureTestResult('graph_genus', True)
    """

    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.graph_genus import GraphGenus
            sage: isinstance(GraphGenus(), GraphGenus)
            True
        """
        PythonModule.__init__(self, "graph_genus", spkg="graph_genus", type="optional")


def all_features():
    return [GraphGenus()]
