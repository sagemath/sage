r"""
Matching covered graphs

This module implements functions and operations pertaing to matching covered graphs.

AUTHORS:

- Janmenjaya Panda (2024-06-14): initial version

"""
# ****************************************************************************
#         Copyright (C) 2024 Janmenjaya Panda <janmenjaya.panda.22@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from .graph import Graph


class MatchingCoveredGraph(Graph):
    def __init__(self):
        pass

    def _upgrade_from_graph(self):
        pass

    def __repr__(self):
        pass

    def plot(self):
        pass

    def maximal_barrier(self):
        pass

    def canonical_partition(self):
        pass

    def tight_cut_decomposition(self):
        pass

    def bricks_and_braces(self):
        pass

    def number_of_bricks(self):
        pass

    def number_of_braces(self):
        pass

    def number_of_petersen_bricks(self):
        pass

    def is_brick(self):
        pass

    def is_brace(self):
        pass

    def is_removable_edge(self):
        pass

    def removable_edges(self):
        pass

    def is_removable_doubleton(self):
        pass

    def removable_doubletons(self):
        pass

    def retract(self):
        pass

    def ear_decomposition(self):
        pass

    def is_b_invariant_edge(self):
        pass

    def is_thin_edge(self):
        pass

    def is_strictly_thin_edge(self):
        pass

    def is_mccuaig_brace(self):
        pass

    def brace_generation_sequence(self):
        pass

    def is_norine_thomas_brick(self):
        pass

    def brick_generation_sequence(self):
        pass