r"""
Forest Posets

AUTHORS:

- Stefan Grosser (06-2020): initial implementation
"""

# ****************************************************************************
#       Copyright (C) 2020 Stefan Grosser <stefan.grosser1@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.combinat.posets.posets import FinitePoset
from sage.combinat.posets.linear_extensions import LinearExtensionsOfForest


class ForestPoset(FinitePoset):
    r"""
    A forest poset is a poset where the underlying Hasse diagram and is
    directed acyclic graph.
    """
    _lin_ext_type = LinearExtensionsOfForest
    _desc = 'Finite forest poset'
