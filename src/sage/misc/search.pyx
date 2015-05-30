"""
Searching a sorted list

This is like the bisect library module, but also returns whether
or not the element is in the list, which saves having to do an
extra comparison.   Also, the function names make more sense.
"""

#*****************************************************************************
#  Sage: System for Algebra and Geometry Computation (c) William Stein, 2004
#
#       Copyright (C) 2004 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import bisect

def search(object v, object x):
    """
    Return (True,i) where i is such that v[i] == x if there is such an i,
    or (False,j) otherwise, where j is the position that a should be inserted
    so that v remains sorted.

    INPUT:
        v -- a list, which is assumed sorted
        x -- Python object
    OUTPUT:
        bool, int
    """
    # This is implemented using the built-in bisect module.  I tried to
    # implement this directly (using binary search) in Pyrex, but it
    # was about half as fast.  (Maybe I did things in an inefficient
    # manner?)   In any case, simply using the bisect module works well,
    # and for some strange reason in tests this is faster than using
    # the bisect module directly from Python (!?).
    cdef int i, n
    i = bisect.bisect_left(v, x)
    n = len(v)
    if i >= n:
        return False, n
    if v[i] != x:
        return False, i
    return True, i

