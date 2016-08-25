r"""
Finite combinatorial classes
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
from __future__ import absolute_import

from .combinat import CombinatorialClass

class FiniteCombinatorialClass(CombinatorialClass):
    """
    INPUT:
     - l a list or iterable

    Returns l, wrapped as a combinatorial class

    EXAMPLES::

        sage: F = FiniteCombinatorialClass([1,2,3])
        sage: F.list()
        [1, 2, 3]
        sage: F.cardinality()
        3
        sage: F.random_element()
        1
        sage: F.first()
        1
        sage: F.last()
        3
    """
    def __init__(self, l):
        """
        TESTS::

            sage: F = FiniteCombinatorialClass([1,2,3])
            sage: F == loads(dumps(F))
            True
        """
        self.l = list(l) # Probably would be better to use a tuple

    def _element_constructor_(self, x):
        """
        EXAMPLES::

            sage: F = FiniteCombinatorialClass([1,2,3])
            sage: F._element_constructor_(1)
            1
            sage: F(1)
            1
        """
        return x

    def __repr__(self):
        """
        TESTS::

            sage: F = FiniteCombinatorialClass([1,2,3])
            sage: repr(F)
            'Combinatorial class with elements in [1, 2, 3]'
        """
        return "Combinatorial class with elements in %s"%self.l

    def __contains__(self, x):
        """
        EXAMPLES::

            sage: F = FiniteCombinatorialClass([1,2,3])
            sage: 1 in F
            True
            sage: 2 in F
            True
            sage: 4 in F
            False
            sage: ZZ in F
            False
        """
        return x in self.l

    def list(self):
        """
        TESTS::

            sage: F = FiniteCombinatorialClass([1,2,3])
            sage: F.list()
            [1, 2, 3]
        """
        return self.l

    def cardinality(self):
        """
        EXAMPLES::

            sage: F = FiniteCombinatorialClass([1,2,3])
            sage: F.cardinality()
            3
        """
        return len(self.l)


    def __getitem__(self, i): # TODO: optimize
        """
        EXAMPLES::

            sage: F = FiniteCombinatorialClass(["a", "b", "c"])
            sage: F[2]
            'c'
        """
        return self.l[i]

    def keys(self):
        """
        EXAMPLES::

            sage: F = FiniteCombinatorialClass([1,2,3])
            sage: F.keys()
            [0, 1, 2]
        """
        return range(len(self.l))

# Backward compatibility pointer
# Needed for unpickling.
FiniteCombinatorialClass_l = FiniteCombinatorialClass
