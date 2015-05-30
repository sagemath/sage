"""
Deprecated low-level permutations
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

def PermutationsNK(n, k):
    r"""
    This is deprecated in :trac:`16472`. Use :class:`Permutations` instead
    (or ``itertools.permutations`` for iteration).

    EXAMPLES::

        sage: from sage.combinat.permutation_nk import PermutationsNK
        sage: P = PermutationsNK(10,4)
        doctest:...: DeprecationWarning: PermutationsNK is deprecated. Please
        use Permutations instead (or itertools.permutations for iteration).
        See http://trac.sagemath.org/16472 for details.
        sage: [ p for p in PermutationsNK(3,2)]
        [[0, 1], [0, 2], [1, 0], [1, 2], [2, 0], [2, 1]]
        sage: PermutationsNK(3,2).cardinality()
        6
        sage: PermutationsNK(5,4).cardinality()
        120
    """
    from sage.misc.superseded import deprecation
    deprecation(16472, "PermutationsNK is deprecated. Please use Permutations instead (or itertools.permutations for iteration).")
    from sage.combinat.permutation import Permutations
    return Permutations(range(n),k)

