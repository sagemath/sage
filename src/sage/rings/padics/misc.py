"""
Miscellaneous Functions

This file contains some miscellaneous functions used by p-adics.

- ``min`` -- a version of ``min`` that returns `\infty` on empty input.
- ``max`` -- a version of ``max`` that returns `-\infty` on empty input.

AUTHORS:

- David Roe
"""

#*****************************************************************************
#       Copyright (C) 2007-2013 David Roe <roed.math@gmail.com>
#                               William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __builtin__ import min as python_min
from __builtin__ import max as python_max
from sage.rings.infinity import infinity

def min(*L):
    """
    Returns the minimum of the inputs, where the minimum of the empty
    list is ``infinity``.

    EXAMPLES::

        sage: from sage.rings.padics.misc import min
        sage: min()
        +Infinity
        sage: min(2,3)
        2
    """
    if len(L) == 1 and isinstance(L[0], (list, tuple)):
        L = L[0]
    try:
        return python_min(L)
    except ValueError:
        return infinity

def max(*L):
    """
    Returns the maximum of the inputs, where the maximum of the empty
    list is ``-infinity``.

    EXAMPLES::

        sage: from sage.rings.padics.misc import max
        sage: max()
        -Infinity
        sage: max(2,3)
        3
    """
    if len(L) == 1 and isinstance(L[0], (list, tuple)):
        L = L[0]
    try:
        return python_max(L)
    except ValueError:
        return -infinity
