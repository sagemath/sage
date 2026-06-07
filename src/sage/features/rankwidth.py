r"""
Feature for the rank-width graph decomposition
"""

# ****************************************************************************
#       Copyright (C) 2026 Michael Orlitzky <michael@orlitzky.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.config import rankwidth_enabled
from sage.features import PythonModule
from sage.features.build_feature import BuildFeature


class RankWidth(BuildFeature):
    r"""
    A :class:`~sage.features.Feature` indicating whether or not the
    rank-width graph decomposition is available.

    This is a proxy for the build-time availability of the rankwidth
    (librw) library from the "rw" package; the associated module
    :mod:`sage.graphs.graph_decompositions.rankwidth` is built if and
    only if the library is detected, usable, and enabled.

    EXAMPLES::

        sage: from sage.features.rankwidth import RankWidth
        sage: RankWidth().is_present()  # needs rankwidth
        FeatureTestResult('rankwidth', True)
        sage: RankWidth().is_present()  # needs !rankwidth
        FeatureTestResult('rankwidth', False)

    """
    _enabled_in_build = rankwidth_enabled

    def __init__(self):
        r"""
        EXAMPLES::

            sage: from sage.features.rankwidth import RankWidth
            sage: RankWidth()
            Feature('rankwidth')

        """
        super().__init__("rankwidth", spkg="rw", type="standard")

    def is_present_at_runtime(self):
        r"""
        TESTS::

            sage: from sage.features import FeatureTestResult
            sage: from sage.features.rankwidth import RankWidth
            sage: rw = RankWidth()
            sage: result = rw.is_present_at_runtime()
            sage: isinstance(result, FeatureTestResult)
            True
            sage: result  # needs rankwidth
            FeatureTestResult('rankwidth', True)

        """
        cython_modname = "sage.graphs.graph_decompositions.rankwidth"
        result = PythonModule(cython_modname)._is_present()
        result.feature = self
        return result

def all_features():
    return [RankWidth()]
