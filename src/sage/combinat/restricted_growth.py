# sage_setup: distribution = sagemath-combinat
"""
Restricted growth arrays

These combinatorial objects are in bijection with set partitions.
"""
#*****************************************************************************
#       Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>,
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
from sage.combinat.combinat import bell_number
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
import copy


class RestrictedGrowthArrays(UniqueRepresentation, Parent):
    def __init__(self, n):
        """
        EXAMPLES::

            sage: from sage.combinat.restricted_growth import RestrictedGrowthArrays
            sage: R = RestrictedGrowthArrays(3)
            sage: R == loads(dumps(R))
            True
            sage: TestSuite(R).run(skip=['_test_an_element',                            # needs sage.libs.flint
            ....:   '_test_enumerated_set_contains', '_test_some_elements'])
        """
        self._n = n
        self._name = "Restricted growth arrays of size %s" % n
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def __iter__(self):
        """
        EXAMPLES::

            sage: from sage.combinat.restricted_growth import RestrictedGrowthArrays
            sage: R = RestrictedGrowthArrays(3)
            sage: R.list()
            [[1, 0, 0], [2, 0, 1], [2, 1, 0], [2, 1, 1], [3, 1, 2]]
        """
        n = self._n
        a = [1] + [0] * (n - 1)
        m = [0] + [1] * (n - 1)
        while True:
            yield copy.copy(a)
            # Search for maximum i with a[i] != m[i]
            i = n - 1
            while a[i] == m[i] and i >= 0:
                i -= 1
            if i == 0:
                break
            # Update arrays a and m
            a[i] += 1
            a[0] = a[i] + 1 if a[i] == m[i] else m[i]
            mi = a[0]
            for j in range(i + 1, n):
                a[j] = 0
                m[j] = mi

    def cardinality(self):
        """
        EXAMPLES::

            sage: from sage.combinat.restricted_growth import RestrictedGrowthArrays
            sage: R = RestrictedGrowthArrays(6)
            sage: R.cardinality()                                                       # needs sage.libs.flint
            203
        """
        return bell_number(self._n)
