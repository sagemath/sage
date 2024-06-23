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
    def __init__(self, *_, **__):
        raise NotImplementedError()

    def _upgrade_from_graph(self, *_, **__):
        raise NotImplementedError()

    def __repr__(self, *_, **__):
        raise NotImplementedError()

    def plot(self, *_, **__):
        raise NotImplementedError()

    def maximal_barrier(self, *_, **__):
        raise NotImplementedError()

    def canonical_partition(self, *_, **__):
        raise NotImplementedError()

    def tight_cut_decomposition(self, *_, **__):
        raise NotImplementedError()

    def bricks_and_braces(self, *_, **__):
        raise NotImplementedError()

    def number_of_bricks(self, *_, **__):
        raise NotImplementedError()

    def number_of_braces(self, *_, **__):
        raise NotImplementedError()

    def number_of_petersen_bricks(self, *_, **__):
        raise NotImplementedError()

    def is_brick(self, *_, **__):
        raise NotImplementedError()

    def is_brace(self, *_, **__):
        raise NotImplementedError()

    def is_removable_edge(self, *_, **__):
        raise NotImplementedError()

    def removable_edges(self, *_, **__):
        raise NotImplementedError()

    def is_removable_doubleton(self, *_, **__):
        raise NotImplementedError()

    def removable_doubletons(self, *_, **__):
        raise NotImplementedError()

    def retract(self, *_, **__):
        raise NotImplementedError()

    def ear_decomposition(self, *_, **__):
        raise NotImplementedError()

    def is_b_invariant_edge(self, *_, **__):
        raise NotImplementedError()

    def is_thin_edge(self, *_, **__):
        raise NotImplementedError()

    def is_strictly_thin_edge(self, *_, **__):
        raise NotImplementedError()

    def is_mccuaig_brace(self, *_, **__):
        raise NotImplementedError()

    def brace_generation_sequence(self, *_, **__):
        raise NotImplementedError()

    def is_norine_thomas_brick(self, *_, **__):
        raise NotImplementedError()

    def brick_generation_sequence(self, *_, **__):
        raise NotImplementedError()