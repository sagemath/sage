"""
Low-level multichoose
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>
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
from sage.arith.all import binomial
import sage.misc.prandom as rnd

class MultichooseNK(CombinatorialClass):
    def __init__(self, n, k):
        """
        TESTS::

            sage: a = MultichooseNK(3,2)
            sage: a == loads(dumps(a))
            True
        """
        self._n = n
        self._k = k

    def cardinality(self):
        """
        Returns the number of multichoices of k things from a list of n
        things.

        EXAMPLES::

            sage: MultichooseNK(3,2).cardinality()
            6
        """
        n,k = self._n, self._k
        return binomial(n+k-1,k)

    def __iter__(self):
        """
        An iterator for all multichoices of k things from range(n).

        EXAMPLES::

            sage: [c for c in MultichooseNK(3,2)]
            [[0, 0], [0, 1], [0, 2], [1, 1], [1, 2], [2, 2]]
        """
        n,k = self._n, self._k
        dif = 0
        if k == 0:
            yield []
            return

        if n < 1+(k-1)*dif:
            return
        else:
            subword = [ i*dif for i in range(k) ]

        yield subword[:]
        finished = False

        while not finished:
            #Find the biggest element that can be increased
            if subword[-1] < n-1:
                subword[-1] += 1
                yield subword[:]
                continue

            finished = True
            for i in reversed(range(k-1)):
                if subword[i]+dif < subword[i+1]:
                    subword[i] += 1
                    #Reset the bigger elements
                    for j in range(1,k-i):
                        subword[i+j] = subword[i]+j*dif
                    yield subword[:]
                    finished = False
                    break

        return

    def random_element(self):
        """
        Returns a random multichoice of k things from range(n).

        EXAMPLES::

            sage: MultichooseNK(5,2).random_element()
            [0, 2]
            sage: MultichooseNK(5,2).random_element()
            [0, 1]
        """
        n,k = self._n, self._k
        rng = range(n)
        r = []
        for i in range(k):
            r.append( rnd.choice(rng))

        r.sort()
        return r
