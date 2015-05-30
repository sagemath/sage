"""
Derecated splits

Authors:

- Mike Hansen (2007): original version

- Vincent Delecroix (2014): deprecation
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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

def SplitNK(n, k):
    """
    Returns the combinatorial class of splits of a the set range(n)
    into a set of size k and a set of size n-k.

    This was deprecated in :trac:`10534`. Use instead
    :class:`OrderedSetPartitions`.

    EXAMPLES::

        sage: from sage.combinat.split_nk import SplitNK
        sage: S = SplitNK(5,2)
        doctest:...: DeprecationWarning:  SplitNk is deprecated and will be
        removed. Use OrderedSetPartitions instead.
        See http://trac.sagemath.org/10534 for details.
        sage: S
        Ordered set partitions of {0, 1, 2, 3, 4} into parts of size [2, 3]
        sage: S.first()
        [{0, 1}, {2, 3, 4}]
        sage: S.last()
        [{3, 4}, {0, 1, 2}]
        sage: S.list()
        [[{0, 1}, {2, 3, 4}],
         [{0, 2}, {1, 3, 4}],
         [{0, 3}, {1, 2, 4}],
         [{0, 4}, {1, 2, 3}],
         [{1, 2}, {0, 3, 4}],
         [{1, 3}, {0, 2, 4}],
         [{1, 4}, {0, 2, 3}],
         [{2, 3}, {0, 1, 4}],
         [{2, 4}, {0, 1, 3}],
         [{3, 4}, {0, 1, 2}]]
    """
    from sage.misc.superseded import deprecation
    from sage.combinat.set_partition_ordered import OrderedSetPartitions
    deprecation(10534, "SplitNk is deprecated and will be removed. Use OrderedSetPartitions instead.")
    return OrderedSetPartitions(range(n), [k, n-k])
