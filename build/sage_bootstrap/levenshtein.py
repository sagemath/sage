# -*- coding: utf-8 -*-
u"""
Levenshtein Distance

The Levenshtein distance between two words is the minimal number of
edits that turn one word into the other. Here, "edit" means a
single-letter addition, single-letter deletion, or exchange of a
letter with another letter.

https://en.wikipedia.org/wiki/Levenshtein_distance

EXAMPLES::

    >>> from sage_bootstrap.levenshtein import Levenshtein
    >>> Levenshtein(5)(u'Queensryche', u'Queensrÿche')
    1
"""

# ****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


class DistanceExceeded(Exception):
    pass


class Levenshtein(object):

    def __init__(self, limit):
        """
        Levenshtein Distance with Maximum Distance Cutoff

        Args:
            limit (int): if the distance exceeds the limit, a
                :class:`DistanceExceeded` is raised and the
                computation is aborted.

        EXAMPLES::

            >>> from sage_bootstrap.levenshtein import Levenshtein
            >>> lev3 = Levenshtein(3)
            >>> lev3(u'saturday', u'sunday')
            3
            >>> lev3(u'kitten', u'sitting')
            3
            >>> lev2 = Levenshtein(2)
            >>> lev2(u'kitten', u'sitting')
            Traceback (most recent call last):
            ...
            DistanceExceeded
        """
        self._limit = limit

    def __call__(self, a, b):
        """
        calculate the levenshtein distance

        args:
            a,b (str): the two strings to compare

        returns:
            int: the Levenshtein distance if it is less or equal to
                the distance limit.

        Example::

            >>> from app.scoring.levenshtein import Levenshtein
            >>> lev3 = Levenshtein(3)
            >>> lev3(u'Saturday', u'Sunday')
            3
        """
        n, m = len(a), len(b)
        if n > m:
            # Optimization to use O(min(n,m)) space
            a, b, n, m = b, a, m, n
        curr = range(n + 1)
        for i in range(1, m + 1):
            prev, curr = curr, [i] + [0] * n
            for j in range(1, n + 1):
                cost_add, cost_del = prev[j] + 1, curr[j - 1] + 1
                cost_change = prev[j - 1]
                if a[j - 1] != b[i - 1]:
                    cost_change += 1
                curr[j] = min(cost_add, cost_del, cost_change)
            if min(curr) > self._limit:
                raise DistanceExceeded
        if curr[n] > self._limit:
            raise DistanceExceeded
        return curr[n]
