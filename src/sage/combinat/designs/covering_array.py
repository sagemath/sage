r"""
Covering Arrays (CA)

A Covering Array, denoted CA(N,k,v,t), is an n by k array with entries from a
set of v elements with theproperty that in every selection of t columns, each
row contains every sequence of t-elements at least once.

An Orthogonal Array, denoted OA(N,k,v,t) is a covering array with the
property that each row contains every sequence of t-elements exactly once

REFERENCES:

.. [Col2004] \C.J. Colbourn. “Combinatorial aspects of covering arrays”.
            Matematiche (Catania) 59 (2004), pp. 125–172.

.. [Sher2006] \G.B. Sherwood, S.S Martirosyan, and C.J. Colbourn, "Covering
              arrays of higher strength from permutation vectors". J. Combin.
              Designs, 14 (2006) pp. 202-213.

.. [Wal2007] \R.A. Walker II, and C.J. Colbourn, "Perfect Hash Families:
             Constructions and Existence". J. Math. Crypt. 1 (2007),
             pp.125-150

AUTHORS:

- Aaron Dwyer and brett stevens (2022): initial version

Classes and methods
-------------------
"""

# ****************************************************************************
#       Copyright (C) 2022 Aaron Dwyer <aarondwyer@cmail.carleton.ca>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.finite_rings.finite_field_constructor import GF
import itertools
import copy

class CoveringArray():
    r"""
    Covering Array (CA)

    INPUT:

    - ``Array`` -- The N by k array itself stored as a tuple of tuples.
      The N and k parameters are derived from this inputted array

    - ``SymbolSet`` -- The collection of symbols that is used in ``Array``.
      If left blank, then a symbol set will be assumed by checking for each
      unique entry in the given ``Array``. In such a case it will be stored
      as a list of symbols but any appropriate object may be used as long as
      it has `len()` as a method

    EXAMPLES::

        sage: from sage.combinat.designs.covering_array import CoveringArray
        sage: C = (('a','a','a','b'),\
                   ('a','a','b','a'),\
                   ('a','b','a','a'),\
                   ('b','a','a','a'),\
                   ('b','b','b','b'))
        sage: CoveringArray(C)
        A 5 by 4 Covering Array with entries from ['a', 'b']

        sage: C = ((0,0,0,0,0,0,0,0,0,0),\
                  (1,1,1,1,1,1,1,1,1,1),\
                  (1,1,1,0,1,0,0,0,0,1),\
                  (1,0,1,1,0,1,0,1,0,0),\
                  (1,0,0,0,1,1,1,0,0,0),\
                  (0,1,1,0,0,1,0,0,1,0),\
                  (0,0,1,0,1,0,1,1,1,0),\
                  (1,1,0,1,0,0,1,0,1,0),\
                  (0,0,0,1,1,1,0,0,1,1),\
                  (0,0,1,1,0,0,1,0,0,1),\
                  (0,1,0,1,1,0,0,1,0,0),\
                  (1,0,0,0,0,0,0,1,1,1),\
                  (0,1,0,0,0,1,1,1,0,1))
        sage: CoveringArray(C,[0, 1])
        A 13 by 10 Covering Array with entries from [0, 1]

        sage: C = ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0),\
                   (0, 0, 0, 0, 1, 1, 1, 1, 1, 1),\
                   (0, 1, 1, 1, 0, 0, 0, 1, 1, 1),\
                   (1, 0, 1, 1, 0, 1, 1, 0, 0, 1),\
                   (1, 1, 0, 1, 1, 0, 1, 0, 1, 0),\
                   (1, 1, 1, 0, 1, 1, 0, 1, 0, 0))
        sage: CoveringArray(C,GF(2))
        A 6 by 10 Covering Array with entries from Finite Field of size 2

        sage: C = ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),\
                   (0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),\
                   (0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1),\
                   (1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1),\
                   (1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1),\
                   (1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0),\
                   (1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0))
        sage: CoveringArray(C)
        A 7 by 15 Covering Array with entries from [0, 1]

    """
    def __init__(self, Array, SymbolSet=None):
        r"""
        Constructor function

    EXAMPLES::

        sage: from sage.combinat.designs.covering_array import CoveringArray
        sage: C = (('a','a','a','b'),\
                   ('a','a','b','a'),\
                   ('a','b','a','a'),\
                   ('b','a','a','a'),\
                   ('b','b','b','b'))
        sage: CoveringArray(C)
        A 5 by 4 Covering Array with entries from ['a', 'b']

        sage: C = ((0,0,0,0,0,0,0,0,0,0),\
                  (1,1,1,1,1,1,1,1,1,1),\
                  (1,1,1,0,1,0,0,0,0,1),\
                  (1,0,1,1,0,1,0,1,0,0),\
                  (1,0,0,0,1,1,1,0,0,0),\
                  (0,1,1,0,0,1,0,0,1,0),\
                  (0,0,1,0,1,0,1,1,1,0),\
                  (1,1,0,1,0,0,1,0,1,0),\
                  (0,0,0,1,1,1,0,0,1,1),\
                  (0,0,1,1,0,0,1,0,0,1),\
                  (0,1,0,1,1,0,0,1,0,0),\
                  (1,0,0,0,0,0,0,1,1,1),\
                  (0,1,0,0,0,1,1,1,0,1))
        sage: CoveringArray(C,[0, 1])
        A 13 by 10 Covering Array with entries from [0, 1]

        sage: C = ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0),\
                   (0, 0, 0, 0, 1, 1, 1, 1, 1, 1),\
                   (0, 1, 1, 1, 0, 0, 0, 1, 1, 1),\
                   (1, 0, 1, 1, 0, 1, 1, 0, 0, 1),\
                   (1, 1, 0, 1, 1, 0, 1, 0, 1, 0),\
                   (1, 1, 1, 0, 1, 1, 0, 1, 0, 0))
        sage: CoveringArray(C,GF(2))
        A 6 by 10 Covering Array with entries from Finite Field of size 2

        sage: C = ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),\
                   (0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),\
                   (0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1),\
                   (1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1),\
                   (1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1),\
                   (1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0),\
                   (1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0))
        sage: CoveringArray(C)
        A 7 by 15 Covering Array with entries from [0, 1]

    """
        #From the array input, grab the dimensions of the array
        N=len(Array)
        self.__n=N
        k=len(Array[0])
        self.__k=k

        #Array input is a tuple of tuples, the first thing to do is to sort
        #the tuples lexicographically increasing
        L=list(Array)
        L.sort()
        SortedArray=tuple(L)
        self.__array=SortedArray

        for row in Array:
            assert len(row)==len(Array[0]), "Not all rows have same length"

        #If no symbol set is given, then it may be assumed from what
        #symbols are in the array by flattening the array and counting the
        #number of unique entries.
        if SymbolSet is None:
            SymbolSet = list({x for l in Array for x in l})
            SymbolSet.sort()
        self.__sset=SymbolSet

    def numrows(self):
        r"""
        Return the number of rows, N, of the covering array

        EXAMPLES::

            sage: from sage.combinat.designs.covering_array import CoveringArray
            sage: C = ((1,1,1,0),\
                       (1,1,0,1),\
                       (1,0,1,1),\
                       (0,1,1,1),\
                       (0,0,0,0))
            sage: CA = CoveringArray(C,GF(2))
            sage: CA.numrows()
            5
        """
        return self.__n

    def numcols(self):
        r"""
        Returns the number of columns, k, of the covering array

        EXAMPLES::

            sage: from sage.combinat.designs.covering_array import CoveringArray
            sage: C = ((1,1,1,0),\
                       (1,1,0,1),\
                       (1,0,1,1),\
                       (0,1,1,1),\
                       (0,0,0,0))
            sage: CA = CoveringArray(C,GF(2))
            sage: CA.numcols()
            4
        """
        return self.__k

    def symbolset(self):
        r"""
        Return the symbol set of the array.

        EXAMPLES::

            sage: from sage.combinat.designs.covering_array import CoveringArray
            sage: C = ((1,1,1,0),\
                       (1,1,0,1),\
                       (1,0,1,1),\
                       (0,1,1,1),\
                       (0,0,0,0))
            sage: CA = CoveringArray(C,GF(2))
            sage: CA.symbolset()
            Finite Field of size 2
            sage: CA = CoveringArray(C)
            sage: CA.symbolset()
            [0, 1]
        """
        return self.__sset

    def is_covering_array(self,strength):
        r"""
        Check whether the tuple of tuples in ``Array`` forms a covering array
        with the given strength.

        EXAMPLES::

            sage: from sage.combinat.designs.covering_array import CoveringArray
            sage: C1 = ((1,1,1,0),\
                    (1,1,0,1),\
                    (1,0,1,1),\
                    (0,1,1,1),\
                    (0,0,0,0))
            sage: CA1 = CoveringArray(C1,GF(2))
            sage: C2 = ((1, 0, 0, 2, 0, 2, 1, 2, 2, 1, 0, 2, 2),\
                    (1, 1, 0, 0, 2, 0, 2, 1, 2, 2, 1, 0, 2),\
                    (1, 1, 1, 0, 0, 2, 0, 2, 1, 2, 2, 1, 0),\
                    (0, 1, 1, 1, 0, 0, 2, 0, 2, 1, 2, 2, 1),\
                    (2, 0, 1, 1, 1, 0, 0, 2, 0, 2, 1, 2, 2),\
                    (1, 2, 0, 1, 1, 1, 0, 0, 2, 0, 2, 1, 2),\
                    (1, 1, 2, 0, 1, 1, 1, 0, 0, 2, 0, 2, 1),\
                    (2, 1, 1, 2, 0, 1, 1, 1, 0, 0, 2, 0, 2),\
                    (1, 2, 1, 1, 2, 0, 1, 1, 1, 0, 0, 2, 0),\
                    (0, 1, 2, 1, 1, 2, 0, 1, 1, 1, 0, 0, 2),\
                    (1, 0, 1, 2, 1, 1, 2, 0, 1, 1, 1, 0, 0),\
                    (0, 1, 0, 1, 2, 1, 1, 2, 0, 1, 1, 1, 0),\
                    (0, 0, 1, 0, 1, 2, 1, 1, 2, 0, 1, 1, 1),\
                    (2, 0, 0, 1, 0, 1, 2, 1, 1, 2, 0, 1, 1),\
                    (2, 2, 0, 0, 1, 0, 1, 2, 1, 1, 2, 0, 1),\
                    (2, 2, 2, 0, 0, 1, 0, 1, 2, 1, 1, 2, 0),\
                    (0, 2, 2, 2, 0, 0, 1, 0, 1, 2, 1, 1, 2),\
                    (1, 0, 2, 2, 2, 0, 0, 1, 0, 1, 2, 1, 1),\
                    (2, 1, 0, 2, 2, 2, 0, 0, 1, 0, 1, 2, 1),\
                    (2, 2, 1, 0, 2, 2, 2, 0, 0, 1, 0, 1, 2),\
                    (1, 2, 2, 1, 0, 2, 2, 2, 0, 0, 1, 0, 1),\
                    (2, 1, 2, 2, 1, 0, 2, 2, 2, 0, 0, 1, 0),\
                    (0, 2, 1, 2, 2, 1, 0, 2, 2, 2, 0, 0, 1),\
                    (2, 0, 2, 1, 2, 2, 1, 0, 2, 2, 2, 0, 0),\
                    (0, 2, 0, 2, 1, 2, 2, 1, 0, 2, 2, 2, 0),\
                    (0, 0, 2, 0, 2, 1, 2, 2, 1, 0, 2, 2, 2),\
                    (1, 1, 0, 2, 1, 1, 2, 1, 0, 1, 0, 0, 2),\
                    (1, 1, 1, 0, 2, 1, 1, 2, 1, 0, 1, 0, 0),\
                    (0, 1, 1, 1, 0, 2, 1, 1, 2, 1, 0, 1, 0),\
                    (0, 0, 1, 1, 1, 0, 2, 1, 1, 2, 1, 0, 1),\
                    (2, 0, 0, 1, 1, 1, 0, 2, 1, 1, 2, 1, 0),\
                    (0, 2, 0, 0, 1, 1, 1, 0, 2, 1, 1, 2, 1),\
                    (2, 0, 2, 0, 0, 1, 1, 1, 0, 2, 1, 1, 2),\
                    (1, 2, 0, 2, 0, 0, 1, 1, 1, 0, 2, 1, 1),\
                    (2, 1, 2, 0, 2, 0, 0, 1, 1, 1, 0, 2, 1),\
                    (2, 2, 1, 2, 0, 2, 0, 0, 1, 1, 1, 0, 2),\
                    (1, 2, 2, 1, 2, 0, 2, 0, 0, 1, 1, 1, 0),\
                    (0, 1, 2, 2, 1, 2, 0, 2, 0, 0, 1, 1, 1),\
                    (2, 0, 1, 2, 2, 1, 2, 0, 2, 0, 0, 1, 1),\
                    (2, 2, 0, 1, 2, 2, 1, 2, 0, 2, 0, 0, 1),\
                    (2, 2, 2, 0, 1, 2, 2, 1, 2, 0, 2, 0, 0),\
                    (0, 2, 2, 2, 0, 1, 2, 2, 1, 2, 0, 2, 0),\
                    (0, 0, 2, 2, 2, 0, 1, 2, 2, 1, 2, 0, 2),\
                    (1, 0, 0, 2, 2, 2, 0, 1, 2, 2, 1, 2, 0),\
                    (0, 1, 0, 0, 2, 2, 2, 0, 1, 2, 2, 1, 2),\
                    (1, 0, 1, 0, 0, 2, 2, 2, 0, 1, 2, 2, 1),\
                    (2, 1, 0, 1, 0, 0, 2, 2, 2, 0, 1, 2, 2),\
                    (1, 2, 1, 0, 1, 0, 0, 2, 2, 2, 0, 1, 2),\
                    (1, 1, 2, 1, 0, 1, 0, 0, 2, 2, 2, 0, 1),\
                    (2, 1, 1, 2, 1, 0, 1, 0, 0, 2, 2, 2, 0),\
                    (0, 2, 1, 1, 2, 1, 0, 1, 0, 0, 2, 2, 2),\
                    (1, 0, 2, 1, 1, 2, 1, 0, 1, 0, 0, 2, 2),\
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            sage: CA2 = CoveringArray(C2,GF(3))
            sage: CA1.is_covering_array(2)
            True
            sage: CA1.is_covering_array(3)
            False
            sage: CA2.is_covering_array(3)
            True
        """
        tupledict={}
        a=[ttuple for ttuple in itertools.product(self.symbolset(),repeat=strength)]
        for item in a:
            tupledict.update({item:0})
        for comb in itertools.combinations(range(self.numcols()),strength):
            wdict=copy.deepcopy(tupledict)
            for row in self.__array:
                wdict[tuple([row[ti] for ti in comb])]+=1
            if 0 in wdict.values():
                return False
        return True

    def strength(self):
        r"""
        Return the strength of the covering array, which is the parameter
        t, such that in any selection of t columns of the array, every
        t tuple appears at least once.

        EXAMPLES::

            sage: from sage.combinat.designs.covering_array import CoveringArray
            sage: C = ((1,1,1,0),\
                       (1,1,0,1),\
                       (1,0,1,1),\
                       (0,1,1,1),\
                       (0,0,0,0))
            sage: CA = CoveringArray(C,GF(2))
            sage: CA.strength()
            2


            sage: C2 = ((1, 0, 0, 2, 0, 2, 1, 2, 2, 1, 0, 2, 2),\
                    (1, 1, 0, 0, 2, 0, 2, 1, 2, 2, 1, 0, 2),\
                    (1, 1, 1, 0, 0, 2, 0, 2, 1, 2, 2, 1, 0),\
                    (0, 1, 1, 1, 0, 0, 2, 0, 2, 1, 2, 2, 1),\
                    (2, 0, 1, 1, 1, 0, 0, 2, 0, 2, 1, 2, 2),\
                    (1, 2, 0, 1, 1, 1, 0, 0, 2, 0, 2, 1, 2),\
                    (1, 1, 2, 0, 1, 1, 1, 0, 0, 2, 0, 2, 1),\
                    (2, 1, 1, 2, 0, 1, 1, 1, 0, 0, 2, 0, 2),\
                    (1, 2, 1, 1, 2, 0, 1, 1, 1, 0, 0, 2, 0),\
                    (0, 1, 2, 1, 1, 2, 0, 1, 1, 1, 0, 0, 2),\
                    (1, 0, 1, 2, 1, 1, 2, 0, 1, 1, 1, 0, 0),\
                    (0, 1, 0, 1, 2, 1, 1, 2, 0, 1, 1, 1, 0),\
                    (0, 0, 1, 0, 1, 2, 1, 1, 2, 0, 1, 1, 1),\
                    (2, 0, 0, 1, 0, 1, 2, 1, 1, 2, 0, 1, 1),\
                    (2, 2, 0, 0, 1, 0, 1, 2, 1, 1, 2, 0, 1),\
                    (2, 2, 2, 0, 0, 1, 0, 1, 2, 1, 1, 2, 0),\
                    (0, 2, 2, 2, 0, 0, 1, 0, 1, 2, 1, 1, 2),\
                    (1, 0, 2, 2, 2, 0, 0, 1, 0, 1, 2, 1, 1),\
                    (2, 1, 0, 2, 2, 2, 0, 0, 1, 0, 1, 2, 1),\
                    (2, 2, 1, 0, 2, 2, 2, 0, 0, 1, 0, 1, 2),\
                    (1, 2, 2, 1, 0, 2, 2, 2, 0, 0, 1, 0, 1),\
                    (2, 1, 2, 2, 1, 0, 2, 2, 2, 0, 0, 1, 0),\
                    (0, 2, 1, 2, 2, 1, 0, 2, 2, 2, 0, 0, 1),\
                    (2, 0, 2, 1, 2, 2, 1, 0, 2, 2, 2, 0, 0),\
                    (0, 2, 0, 2, 1, 2, 2, 1, 0, 2, 2, 2, 0),\
                    (0, 0, 2, 0, 2, 1, 2, 2, 1, 0, 2, 2, 2),\
                    (1, 1, 0, 2, 1, 1, 2, 1, 0, 1, 0, 0, 2),\
                    (1, 1, 1, 0, 2, 1, 1, 2, 1, 0, 1, 0, 0),\
                    (0, 1, 1, 1, 0, 2, 1, 1, 2, 1, 0, 1, 0),\
                    (0, 0, 1, 1, 1, 0, 2, 1, 1, 2, 1, 0, 1),\
                    (2, 0, 0, 1, 1, 1, 0, 2, 1, 1, 2, 1, 0),\
                    (0, 2, 0, 0, 1, 1, 1, 0, 2, 1, 1, 2, 1),\
                    (2, 0, 2, 0, 0, 1, 1, 1, 0, 2, 1, 1, 2),\
                    (1, 2, 0, 2, 0, 0, 1, 1, 1, 0, 2, 1, 1),\
                    (2, 1, 2, 0, 2, 0, 0, 1, 1, 1, 0, 2, 1),\
                    (2, 2, 1, 2, 0, 2, 0, 0, 1, 1, 1, 0, 2),\
                    (1, 2, 2, 1, 2, 0, 2, 0, 0, 1, 1, 1, 0),\
                    (0, 1, 2, 2, 1, 2, 0, 2, 0, 0, 1, 1, 1),\
                    (2, 0, 1, 2, 2, 1, 2, 0, 2, 0, 0, 1, 1),\
                    (2, 2, 0, 1, 2, 2, 1, 2, 0, 2, 0, 0, 1),\
                    (2, 2, 2, 0, 1, 2, 2, 1, 2, 0, 2, 0, 0),\
                    (0, 2, 2, 2, 0, 1, 2, 2, 1, 2, 0, 2, 0),\
                    (0, 0, 2, 2, 2, 0, 1, 2, 2, 1, 2, 0, 2),\
                    (1, 0, 0, 2, 2, 2, 0, 1, 2, 2, 1, 2, 0),\
                    (0, 1, 0, 0, 2, 2, 2, 0, 1, 2, 2, 1, 2),\
                    (1, 0, 1, 0, 0, 2, 2, 2, 0, 1, 2, 2, 1),\
                    (2, 1, 0, 1, 0, 0, 2, 2, 2, 0, 1, 2, 2),\
                    (1, 2, 1, 0, 1, 0, 0, 2, 2, 2, 0, 1, 2),\
                    (1, 1, 2, 1, 0, 1, 0, 0, 2, 2, 2, 0, 1),\
                    (2, 1, 1, 2, 1, 0, 1, 0, 0, 2, 2, 2, 0),\
                    (0, 2, 1, 1, 2, 1, 0, 1, 0, 0, 2, 2, 2),\
                    (1, 0, 2, 1, 1, 2, 1, 0, 1, 0, 0, 2, 2),\
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
            sage: CA2 = CoveringArray(C2,GF(3))
            sage: CA2.strength()
            3
        """
        finished=False
        strength=1
        while not finished:
            if self.is_covering_array(strength):
                strength+=1
            else:
                strength-=1
                finished=True
        return strength

    def levels(self):
        r"""
        Return the number of levels for the covering array, which is
        the parameter v, such that v is the size of the symbol set of the
        array.

        EXAMPLES::

            sage: from sage.combinat.designs.covering_array import CoveringArray
            sage: C = ((1,1,1,0),\
                       (1,1,0,1),\
                       (1,0,1,1),\
                       (0,1,1,1),\
                       (0,0,0,0))
            sage: CA = CoveringArray(C,GF(2))
            sage: CA.levels()
            2
        """
        return len(self.__sset)

    def array_representation(self):
        r"""
        Return the covering array as a tuple of tuples, but where
        the output is such that each row of the array is sorted in
        lexicographic order

        EXAMPLES::

            sage: from sage.combinat.designs.covering_array import CoveringArray
            sage: C = ((1,1,1,0),\
                       (0,0,0,0),\
                       (1,0,1,1),\
                       (1,1,0,1),\
                       (0,1,1,1),)
            sage: CA = CoveringArray(C,GF(2))
            sage: CA.array_representation()
            ((0, 0, 0, 0), (0, 1, 1, 1), (1, 0, 1, 1), (1, 1, 0, 1), (1, 1, 1, 0))

        """
        return self.__array


    def __repr__(self):
        r"""
        Returns a string that describes self

        EXAMPLES::

            sage: from sage.combinat.designs.covering_array import CoveringArray
            sage: C = (('a','a','a','b'),\
                    ('a','a','b','a'),\
                    ('a','b','a','a'),\
                    ('b','a','a','a'),\
                    ('b','b','b','b'))
            sage: CoveringArray(C)
            A 5 by 4 Covering Array with entries from ['a', 'b']

            sage: C = ((0,0,0,0,0,0,0,0,0,0),\
                    (1,1,1,1,1,1,1,1,1,1),\
                    (1,1,1,0,1,0,0,0,0,1),\
                    (1,0,1,1,0,1,0,1,0,0),\
                    (1,0,0,0,1,1,1,0,0,0),\
                    (0,1,1,0,0,1,0,0,1,0),\
                    (0,0,1,0,1,0,1,1,1,0),\
                    (1,1,0,1,0,0,1,0,1,0),\
                    (0,0,0,1,1,1,0,0,1,1),\
                    (0,0,1,1,0,0,1,0,0,1),\
                    (0,1,0,1,1,0,0,1,0,0),\
                    (1,0,0,0,0,0,0,1,1,1),\
                    (0,1,0,0,0,1,1,1,0,1))
            sage: CoveringArray(C,[0, 1])
            A 13 by 10 Covering Array with entries from [0, 1]

            sage: C = ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0),\
                    (0, 0, 0, 0, 1, 1, 1, 1, 1, 1),\
                    (0, 1, 1, 1, 0, 0, 0, 1, 1, 1),\
                    (1, 0, 1, 1, 0, 1, 1, 0, 0, 1),\
                    (1, 1, 0, 1, 1, 0, 1, 0, 1, 0),\
                    (1, 1, 1, 0, 1, 1, 0, 1, 0, 0))
            sage: CoveringArray(C,GF(2))
            A 6 by 10 Covering Array with entries from Finite Field of size 2

            sage: C = ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),\
                    (0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),\
                    (0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1),\
                    (1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1),\
                    (1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1),\
                    (1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0),\
                    (1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0))
            sage: CoveringArray(C)
            A 7 by 15 Covering Array with entries from [0, 1]

    """
        return 'A {} by {} Covering Array with entries from {}'.format(
            self.numrows(), self.numcols(), self.symbolset())

    __str__=__repr__

    def pprint(self):
        r"""
        Prints the covering array in a format easy for users to read

        EXAMPLES::

            sage: from sage.combinat.designs.covering_array import CoveringArray
            sage: C = ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),\
                       (0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),\
                       (0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1),\
                       (1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1),\
                       (1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1),\
                       (1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0),\
                       (1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0))
            sage: CA = CoveringArray(C)
            sage: CA.pprint()
            (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
            (0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
            (0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
            (1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1)
            (1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1)
            (1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0)
            (1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0)
        """
        for i in self.__array:
            print(str(i))

    def __hash__(self):
        r"""
        Hashs the tuple of tuples and all tuples inside

        EXAMPLES::

            sage: from sage.combinat.designs.covering_array import CoveringArray
            sage: C = ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),\
                       (0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),\
                       (0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1),\
                       (1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1),\
                       (1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1),\
                       (1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0),\
                       (1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0))
            sage: CA = CoveringArray(C)
            sage: hash(CA)
            4367534393624660384
        """
        return hash((self.array_representation(),tuple(self.symbolset())))

    def __eq__(self, other):
        r"""
        Return whether two covering arrays are equal by considering the
        array with rows sorted in lexicographic order

        EXAMPLES::

            sage: from sage.combinat.designs.covering_array import CoveringArray
            sage: C1 = ((1,1,1,0),\
                       (0,0,0,0),\
                       (1,0,1,1),\
                       (1,1,0,1),\
                       (0,1,1,1))
            sage: C2 = ((1,1,1,0),\
                       (1,1,0,1),\
                       (1,0,1,1),\
                       (0,1,1,1),\
                       (0,0,0,0))
            sage: C3 = ((0,0,0,0,0,0,0,0,0,0),\
                       (1,1,1,1,1,1,1,1,1,1),\
                       (1,1,1,0,1,0,0,0,0,1),\
                       (1,0,1,1,0,1,0,1,0,0),\
                       (1,0,0,0,1,1,1,0,0,0),\
                       (0,1,1,0,0,1,0,0,1,0),\
                       (0,0,1,0,1,0,1,1,1,0),\
                       (1,1,0,1,0,0,1,0,1,0),\
                       (0,0,0,1,1,1,0,0,1,1),\
                       (0,0,1,1,0,0,1,0,0,1),\
                       (0,1,0,1,1,0,0,1,0,0),\
                       (1,0,0,0,0,0,0,1,1,1),\
                       (0,1,0,0,0,1,1,1,0,1))
            sage: CA1 = CoveringArray(C1)
            sage: CA2 = CoveringArray(C2)
            sage: CA3 = CoveringArray(C3)
            sage: CA1==CA2
            True
            sage: CA1==CA3
            False
        """
        if self.array_representation() == other.array_representation():
            return True
        else:
            return False

    def __neq__(self, other):
        r"""
        Return whether two covering arrays are not equal by considering
        the array with rows sorted in lexicographic order

        EXAMPLES::

            sage: from sage.combinat.designs.covering_array import CoveringArray
            sage: C1 = ((1,1,1,0),\
                       (0,0,0,0),\
                       (1,0,1,1),\
                       (1,1,0,1),\
                       (0,1,1,1))
            sage: C2 = ((1,1,1,0),\
                       (1,1,0,1),\
                       (1,0,1,1),\
                       (0,1,1,1),\
                       (0,0,0,0))
            sage: C3 = ((0,0,0,0,0,0,0,0,0,0),\
                       (1,1,1,1,1,1,1,1,1,1),\
                       (1,1,1,0,1,0,0,0,0,1),\
                       (1,0,1,1,0,1,0,1,0,0),\
                       (1,0,0,0,1,1,1,0,0,0),\
                       (0,1,1,0,0,1,0,0,1,0),\
                       (0,0,1,0,1,0,1,1,1,0),\
                       (1,1,0,1,0,0,1,0,1,0),\
                       (0,0,0,1,1,1,0,0,1,1),\
                       (0,0,1,1,0,0,1,0,0,1),\
                       (0,1,0,1,1,0,0,1,0,0),\
                       (1,0,0,0,0,0,0,1,1,1),\
                       (0,1,0,0,0,1,1,1,0,1))
            sage: CA1 = CoveringArray(C1)
            sage: CA2 = CoveringArray(C2)
            sage: CA3 = CoveringArray(C3)
            sage: CA1!=CA2
            False
            sage: CA1!=CA3
            True
        """
        if self.array_representation() != other.array_representation():
            return True
        else:
            return False

    def __lt__(self, other):
        r"""
        Return whether one covering array is less than another
        based on the lexicographic order on the rows

        EXAMPLES::

            sage: from sage.combinat.designs.covering_array import CoveringArray
            sage: C1 = ((1,1,1,0),\
                       (0,0,0,0),\
                       (1,0,1,1),\
                       (1,1,0,1),\
                       (0,1,1,1))
            sage: C2 = ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0),\
                        (0, 0, 0, 0, 1, 1, 1, 1, 1, 1),\
                        (0, 1, 1, 1, 0, 0, 0, 1, 1, 1),\
                        (1, 0, 1, 1, 0, 1, 1, 0, 0, 1),\
                        (1, 1, 0, 1, 1, 0, 1, 0, 1, 0),\
                        (1, 1, 1, 0, 1, 1, 0, 1, 0, 0))
            sage: CA1 = CoveringArray(C1)
            sage: CA2 = CoveringArray(C2)
            sage: CA1<CA2
            True
        """
        if self.array_representation() < other.array_representation():
            return True
        else:
            return False

    def __le__(self, other):
        r"""
        Return whether one covering array is less than or
        equal to another based on the lexicographic order on the rows

        EXAMPLES::

            sage: from sage.combinat.designs.covering_array import CoveringArray
            sage: C1 = ((1,1,1,0),\
                       (0,0,0,0),\
                       (1,0,1,1),\
                       (1,1,0,1),\
                       (0,1,1,1))
            sage: C2 = ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0),\
                        (0, 0, 0, 0, 1, 1, 1, 1, 1, 1),\
                        (0, 1, 1, 1, 0, 0, 0, 1, 1, 1),\
                        (1, 0, 1, 1, 0, 1, 1, 0, 0, 1),\
                        (1, 1, 0, 1, 1, 0, 1, 0, 1, 0),\
                        (1, 1, 1, 0, 1, 1, 0, 1, 0, 0))
            sage: CA1 = CoveringArray(C1)
            sage: CA2 = CoveringArray(C2)
            sage: CA1<=CA2
            True
        """
        if self.array_representation() <= other.array_representation():
            return True
        else:
            return False

    def __gt__(self, other):
        r"""
        Return whether one covering array is greater than
        another based on the lexicographic order on the rows

        EXAMPLES::

            sage: from sage.combinat.designs.covering_array import CoveringArray
            sage: C1 = ((1,1,1,0),\
                       (0,0,0,0),\
                       (1,0,1,1),\
                       (1,1,0,1),\
                       (0,1,1,1))
            sage: C2 = ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0),\
                        (0, 0, 0, 0, 1, 1, 1, 1, 1, 1),\
                        (0, 1, 1, 1, 0, 0, 0, 1, 1, 1),\
                        (1, 0, 1, 1, 0, 1, 1, 0, 0, 1),\
                        (1, 1, 0, 1, 1, 0, 1, 0, 1, 0),\
                        (1, 1, 1, 0, 1, 1, 0, 1, 0, 0))
            sage: CA1 = CoveringArray(C1)
            sage: CA2 = CoveringArray(C2)
            sage: CA1>CA2
            False
        """
        if self.array_representation() > other.array_representation():
            return True
        else:
            return False

    def __ge__(self, other):
        r"""
        Return whether one covering array is greater than or
        equal to another based on the lexicographic order on the rows

        EXAMPLES::

            sage: from sage.combinat.designs.covering_array import CoveringArray
            sage: C1 = ((1,1,1,0),\
                       (0,0,0,0),\
                       (1,0,1,1),\
                       (1,1,0,1),\
                       (0,1,1,1))
            sage: C2 = ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0),\
                        (0, 0, 0, 0, 1, 1, 1, 1, 1, 1),\
                        (0, 1, 1, 1, 0, 0, 0, 1, 1, 1),\
                        (1, 0, 1, 1, 0, 1, 1, 0, 0, 1),\
                        (1, 1, 0, 1, 1, 0, 1, 0, 1, 0),\
                        (1, 1, 1, 0, 1, 1, 0, 1, 0, 0))
            sage: CA1 = CoveringArray(C1)
            sage: CA2 = CoveringArray(C2)
            sage: CA1>=CA2
            False
        """
        if self.array_representation() >= other.array_representation():
            return True
        else:
            return False
