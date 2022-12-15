<<<<<<< HEAD
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

.. NOTES::
    
This is a work in progress, it will be an implementation of a Covering Array 
(CA) class for sagemath. The initial commit will include the definition of 
the Covering Array class and some basic methods to check and return the 
parameters N,k,v,t of the CA, as well as an ordering based on the 
lexicographic ordering of each row.

Later commits will include methods to create CAs from Linear Feedback Shift 
Register (LFSR), Perfect Hash Families and Covering Perfect Hash Families 
(CPHF) as well as recursive methods.
=======
'''
An implementation of a Covering Array (CA) class for sagemath. 
The initial commit will include the definition of the Covering Array class 
and some basic methods to check and return the parameters n,k,v,t of the CA, 
as well as an ordering based on the lexicographic ordering of each row.

Later commits will include methods to create CAs from Perfect Hash Families 
and Covering Perfect Hash Families (CPHF) as well as recursive methods.
>>>>>>> 5a4d10877a5c6414a20d23fcc217872d6287ba27

The Covering Array class may be used as a basis for an Orthogonal Array class
which will be implemented afterwards

<<<<<<< HEAD
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

def is_covering_array(Array,levels,strength): 
    r"""
    Check whether the tuple of tuples in ``Array`` forms a covering array.
    
    INPUT:
    
    - ``Array`` - a tuple of tuples that represents a N x k array
    - ``levels`` - positive integer representing the v value of the CA
    - ``strength`` - positive integer representing the t value of the CA
    
    OUTPUT:
     
    A boolean representing if the input is a covering array or not
    
    EXAMPLES::
        
        sage: from sage.combinat.designs.covering_arrays import is_covering_array
        sage: C = ((1,1,1,0),\
                   (1,1,0,1),\
                   (1,0,1,1),\
                   (0,1,1,1),\
                   (0,0,0,0))
        sage: D = ((1, 0, 0, 2, 0, 2, 1, 2, 2, 1, 0, 2, 2),\
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
        sage: is_covering_array(C,2,2)
        True
        sage: is_covering_array(C,2,3)
        False
        sage: is_covering_array(D,3,3)
        True
    """
    tupledict={}
    a=[ttuple for ttuple in itertools.product(GF(levels),repeat=strength)]
    for item in a:
        tupledict.update({item:0})
    for comb in itertools.combinations(range(len(Array[0])),strength):
        wdict=copy.deepcopy(tupledict)
        for row in Array:
            wdict[tuple([row[ti] for ti in comb])]+=1
        if 0 in wdict.values():
            return False
    return True

class CoveringArray():
    r"""
    Covering Array (CA)
    
    INPUT:
        
    - ``Array`` -- The N by k array itself stored as a tuple of tuples.
      The N and k parameters are derived from this inputted array
    
    - ``strength`` -- The parameter t, such that in any selection of t columns
      of the array, every t tuple appears at least once. If ``None`` then the
      maxiumum t is found by iterating from t=0 until array is not a CA with 
      strength t+1.
    
    - ``levels`` -- The paramter v, such that v is the size of the symbol set 
      of the array. If ``None`` then a v will be assumed by counting number of
      unique elements appear in the array.
    
    EXAMPLES::
        
        sage: from sage.combinat.designs.covering_arrays import CoveringArray
        sage: C = ((1,1,1,0),\
                   (1,1,0,1),\
                   (1,0,1,1),\
                   (0,1,1,1),\
                   (0,0,0,0))          
        sage: CoveringArray(C,strength=2)
        A 5 by 4 Covering Array of strength 2 with 2 levels
        
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
        sage: CoveringArray(C,3,2)
        A 13 by 10 Covering Array of strength 3 with 2 levels
        
        sage: C = ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0),\
                   (0, 0, 0, 0, 1, 1, 1, 1, 1, 1),\
                   (0, 1, 1, 1, 0, 0, 0, 1, 1, 1),\
                   (1, 0, 1, 1, 0, 1, 1, 0, 0, 1),\
                   (1, 1, 0, 1, 1, 0, 1, 0, 1, 0),\
                   (1, 1, 1, 0, 1, 1, 0, 1, 0, 0))
        sage: CoveringArray(C)
        A 6 by 10 Covering Array of strength 2 with 2 levels

        sage: C = ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),\
                   (0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),\
                   (0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1),\
                   (1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1),\
                   (1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1),\
                   (1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0),\
                   (1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0))
        sage: CoveringArray(C)
        A 7 by 15 Covering Array of strength 2 with 2 levels

    """
    def __init__(self, Array, strength=None, levels=None):
        r"""
        Constructor function    
        
        EXAMPLES::
        
            sage: from sage.combinat.designs.covering_arrays import CoveringArray
            sage: C = ((1,1,1,0),\
                       (1,1,0,1),\
                       (1,0,1,1),\
                       (0,1,1,1),\
                       (0,0,0,0))           
            sage: CoveringArray(C,strength=2)
            A 5 by 4 Covering Array of strength 2 with 2 levels

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
            sage: CoveringArray(C,3,2)
            A 13 by 10 Covering Array of strength 3 with 2 levels

            sage: C = ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0),\
                       (0, 0, 0, 0, 1, 1, 1, 1, 1, 1),\
                       (0, 1, 1, 1, 0, 0, 0, 1, 1, 1),\
                       (1, 0, 1, 1, 0, 1, 1, 0, 0, 1),\
                       (1, 1, 0, 1, 1, 0, 1, 0, 1, 0),\
                       (1, 1, 1, 0, 1, 1, 0, 1, 0, 0))
            sage: CoveringArray(C)
            A 6 by 10 Covering Array of strength 2 with 2 levels

            sage: C = ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),\
                       (0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),\
                       (0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1),\
                       (1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1),\
                       (1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1),\
                       (1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0),\
                       (1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0))
            sage: CoveringArray(C)
            A 7 by 15 Covering Array of strength 2 with 2 levels

    """
        #From the array input, grab the dimensions of the array
        N=len(Array)
        self.__n=N
=======
'''

import itertools
import copy

def is_covering_array(Array,v ,t): 
        '''
        Check if the tuple of tuples 'Array' is a Covering Array CA(n,k,v,t).
        
        A CA(n,k,v,t) is an n by k array with entries from a set of v elements
        with the property that in every selection of t columns, each row 
        contains every sequence of t-elements at least once.
        
        
        '''
        
        tupledict={}
        a=[ttuple for ttuple in itertools.product(GF(v),repeat=t)]
        for item in a:
            tupledict.update({item:0})
        for comb in itertools.combinations(range(len(Array[0])),t):
            wdict=copy.deepcopy(tupledict)
            for row in Array:
                wdict[tuple([row[ti] for ti in comb])]+=1
            if 0 in wdict.values():
                return False
        return True

class CoveringArray():
    def __init__(self, Array, t=None, v=None):
        #From the array input, grab the dimensions of the array
        n=len(Array)
        self.__n=n
>>>>>>> 5a4d10877a5c6414a20d23fcc217872d6287ba27
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
        

        #If v is not inputted then the symbol set may be assumed from what 
        #symbols are in the array by flattening the Array and counting the 
        #number of unique entries.
<<<<<<< HEAD
        if levels==None:
            levels = len({x for l in Array for x in l})
        self.__v=levels
        
        #If t is not inputted then try all t in {1,2,3...} until the array is
        #not  a covering array for that t
        if strength==None:
            finished=False
            strength=1
            while finished==False:
                if is_covering_array(Array, levels ,strength):
                    strength+=1
                else:
                    strength-=1
=======
        if v==None:
            v = len({x for l in Array for x in l})
        self.__v=v
        
        #If t is not inputted then try all t in {1,2,3...} until the array is
        #not  a covering array for that t
        if t==None:
            finished=False
            t=1
            while finished==False:
                if is_covering_array(Array, v ,t):
                    t+=1
                else:
                    t-=1
>>>>>>> 5a4d10877a5c6414a20d23fcc217872d6287ba27
                    finished=True
        
        #If t is given, make sure that t is the highest strength of the given
        #Covering Array, that is the Covering Array is not strength t+1
        else:
            finished=False
            while finished==False:
<<<<<<< HEAD
                if is_covering_array(Array, levels ,strength):
                    strength+=1
                else:
                    strength-=1
                    finished=True
        self.__t=strength
        
    def numrows(self):
        r"""
        Return the number of rows, N, of the covering array
        
        EXAMPLES::
               
            sage: from sage.combinat.designs.covering_arrays import CoveringArray
            sage: C = ((1,1,1,0),\
                       (1,1,0,1),\
                       (1,0,1,1),\
                       (0,1,1,1),\
                       (0,0,0,0))         
            sage: CA = CoveringArray(C,2,2)
            sage: CA.numrows()
            5
        """
        return self.__n
        
    def numcols(self):
        r"""
        Returns the number of columns, k, of the covering array
        
        EXAMPLES::
            
            sage: from sage.combinat.designs.covering_arrays import CoveringArray
            sage: C = ((1,1,1,0),\
                       (1,1,0,1),\
                       (1,0,1,1),\
                       (0,1,1,1),\
                       (0,0,0,0))          
            sage: CA = CoveringArray(C,2,2)
            sage: CA.numcols()
            4
        """
        return self.__k
    
    def levels(self):
        r"""
        Return the number of levels for the covering array, which is 
        the paramter v, such that v is the size of the symbol set of the 
        array.
    
        EXAMPLES::
                
            sage: from sage.combinat.designs.covering_arrays import CoveringArray
            sage: C = ((1,1,1,0),\
                       (1,1,0,1),\
                       (1,0,1,1),\
                       (0,1,1,1),\
                       (0,0,0,0))         
            sage: CA = CoveringArray(C,2,2)
            sage: CA.levels()
            2
        """
        return self.__v
        
    def strength(self):
        r"""
        Return the strength of the covering array, which is the paramter 
        t, such that in any selection of t columns of the array, every 
        t tuple appears at least once.
    
        EXAMPLES::
                
            sage: from sage.combinat.designs.covering_arrays import CoveringArray
            sage: C = ((1,1,1,0),\
                       (1,1,0,1),\
                       (1,0,1,1),\
                       (0,1,1,1),\
                       (0,0,0,0))        
            sage: CA = CoveringArray(C,2,2)
            sage: CA.strength()
            2
        """
        return self.__t
    
    def array_representation(self):
        r"""
        Return the covering array as a tuple of tuples, but where
        the output is such that each row of the array is sorted in 
        lexicographic order
        
        EXAMPLES::
                
            sage: from sage.combinat.designs.covering_arrays import CoveringArray
            sage: C = ((1,1,1,0),\
                       (0,0,0,0),\
                       (1,0,1,1),\
                       (1,1,0,1),\
                       (0,1,1,1),)         
            sage: CA = CoveringArray(C,2,2)
            sage: CA.array_representation()
            ((0, 0, 0, 0), (0, 1, 1, 1), (1, 0, 1, 1), (1, 1, 0, 1), (1, 1, 1, 0))
            
        """
=======
                if is_covering_array(Array, v ,t):
                    t+=1
                else:
                    t-=1
                    finished=True
        self.__t=t
        
    def n(self):
        return self.__n
        
    def k(self):
        return self.__k
    
    def v(self):
        return self.__v
        
    def t(self):
        return self.__t
    
    def ArrayRepresentation(self):
>>>>>>> 5a4d10877a5c6414a20d23fcc217872d6287ba27
        return self.__array
    
    
    def __repr__(self):
<<<<<<< HEAD
        r"""
        Returns a string that describes self
        
        EXAMPLES::
        
            sage: from sage.combinat.designs.covering_arrays import CoveringArray
            sage: C = ((1,1,1,0),\
                       (1,1,0,1),\
                       (1,0,1,1),\
                       (0,1,1,1),\
                       (0,0,0,0))           
            sage: CoveringArray(C,strength=2)
            A 5 by 4 Covering Array of strength 2 with 2 levels

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
            sage: CoveringArray(C,3,2)
            A 13 by 10 Covering Array of strength 3 with 2 levels

            sage: C = ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0),\
                       (0, 0, 0, 0, 1, 1, 1, 1, 1, 1),\
                       (0, 1, 1, 1, 0, 0, 0, 1, 1, 1),\
                       (1, 0, 1, 1, 0, 1, 1, 0, 0, 1),\
                       (1, 1, 0, 1, 1, 0, 1, 0, 1, 0),\
                       (1, 1, 1, 0, 1, 1, 0, 1, 0, 0))
            sage: CoveringArray(C)
            A 6 by 10 Covering Array of strength 2 with 2 levels

            sage: C = ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),\
                       (0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),\
                       (0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1),\
                       (1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1),\
                       (1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1),\
                       (1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0),\
                       (1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0))
            sage: CoveringArray(C)
            A 7 by 15 Covering Array of strength 2 with 2 levels


    """
        return 'A {} by {} Covering Array of strength {} with {} levels'.format(
            self.numrows(), self.numcols(), self.strength(), self.levels())
        
    __str__=__repr__
    
    def pprint(self):
        r"""
        Prints the covering array in a format easy for users to read
        
        EXAMPLES::
                
            sage: from sage.combinat.designs.covering_arrays import CoveringArray
            sage: C = ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),\
                       (0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),\
                       (0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1),\
                       (1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1),\
                       (1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1),\
                       (1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0),\
                       (1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0))
            sage: CA = CoveringArray(C,2,2)
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
    
            sage: from sage.combinat.designs.covering_arrays import CoveringArray
            sage: C = ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),\
                       (0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),\
                       (0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1),\
                       (1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1),\
                       (1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1),\
                       (1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0),\
                       (1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0))
            sage: CA = CoveringArray(C,2,2)
            sage: hash(CA)
            -2522140066511050633
        """
        return hash((self.array_representation(),self.strength(),self.levels()))
        
    def __eq__(self, other):
        r"""
        Return whether two covering arrays are equal by considering the 
        array with rows sorted in lexicographic order
        
        EXAMPLES::
        
            sage: from sage.combinat.designs.covering_arrays import CoveringArray
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
=======
        """
        A print method
        
        EXAMPLES::
        
        sage: CA = CoveringArray(5,4,2,2,[[1 , 1 , 1 , 0],[1 , 1 , 0 , 1],[1 , 0 , 1 , 1],[0 , 1 , 1 , 1],[0 , 0 , 0 , 0]])
        sage: CA
        A 5 by 4 Covering Array of strength 2 with 2 levels
        """
        return 'A {} by {} Covering Array of strength {} with {} levels'.format(
            self.nrows, self.ncols, self.strength, self.levels)
        
    __str__=__repr__
    
    def __eq__(self, other):
        if self.ArrayRepresentation() == other.ArrayRepresentation():
>>>>>>> 5a4d10877a5c6414a20d23fcc217872d6287ba27
            return True
        else:
            return False
        
    def __neq__(self, other):
<<<<<<< HEAD
        r"""
        Return whether two covering arrays are not equal by considering 
        the array with rows sorted in lexicographic order
        
        EXAMPLES::
        
            sage: from sage.combinat.designs.covering_arrays import CoveringArray
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
=======
        if self.ArrayRepresentation() != other.ArrayRepresentation():
>>>>>>> 5a4d10877a5c6414a20d23fcc217872d6287ba27
            return True
        else:
            return False
        
    def __lt__(self, other):
<<<<<<< HEAD
        r"""
        Return whether one covering array is less than another
        based on the lexicographic order on the rows
        
        EXAMPLES::
        
            sage: from sage.combinat.designs.covering_arrays import CoveringArray
            sage: C1 = ((1,1,1,0),\
                       (0,0,0,0),\
                       (1,0,1,1),\
                       (1,1,0,1),\
                       (0,1,1,1),)
            sage: C2 = ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0),\
                        (0, 0, 0, 0, 1, 1, 1, 1, 1, 1),\
                        (0, 1, 1, 1, 0, 0, 0, 1, 1, 1),\
                        (1, 0, 1, 1, 0, 1, 1, 0, 0, 1),\
                        (1, 1, 0, 1, 1, 0, 1, 0, 1, 0),\
                        (1, 1, 1, 0, 1, 1, 0, 1, 0, 0))
            sage: CA1 = CoveringArray(C1,2,2)
            sage: CA2 = CoveringArray(C2,2,2)
            sage: CA1<CA2
            True
        """
        if self.array_representation() < other.array_representation():
=======
        if self.ArrayRepresentation() < other.ArrayRepresentation():
>>>>>>> 5a4d10877a5c6414a20d23fcc217872d6287ba27
            return True
        else:
            return False
        
    def __le__(self, other):
<<<<<<< HEAD
        r"""
        Return whether one covering array is less than or 
        equal to another based on the lexicographic order on the rows
        
        EXAMPLES::
        
            sage: from sage.combinat.designs.covering_arrays import CoveringArray
            sage: C1 = ((1,1,1,0),\
                       (0,0,0,0),\
                       (1,0,1,1),\
                       (1,1,0,1),\
                       (0,1,1,1),)
            sage: C2 = ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0),\
                        (0, 0, 0, 0, 1, 1, 1, 1, 1, 1),\
                        (0, 1, 1, 1, 0, 0, 0, 1, 1, 1),\
                        (1, 0, 1, 1, 0, 1, 1, 0, 0, 1),\
                        (1, 1, 0, 1, 1, 0, 1, 0, 1, 0),\
                        (1, 1, 1, 0, 1, 1, 0, 1, 0, 0))
            sage: CA1 = CoveringArray(C1,2,2)
            sage: CA2 = CoveringArray(C2,2,2)
            sage: CA1<=CA2
            True
        """
        if self.array_representation() <= other.array_representation():
=======
        if self.ArrayRepresentation() <= other.ArrayRepresentation():
>>>>>>> 5a4d10877a5c6414a20d23fcc217872d6287ba27
            return True
        else:
            return False
        
    def __gt__(self, other):
<<<<<<< HEAD
        r"""
        Return whether one covering array is greater than 
        another based on the lexicographic order on the rows
        
        EXAMPLES::
        
            sage: from sage.combinat.designs.covering_arrays import CoveringArray
            sage: C1 = ((1,1,1,0),\
                       (0,0,0,0),\
                       (1,0,1,1),\
                       (1,1,0,1),\
                       (0,1,1,1),)
            sage: C2 = ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0),\
                        (0, 0, 0, 0, 1, 1, 1, 1, 1, 1),\
                        (0, 1, 1, 1, 0, 0, 0, 1, 1, 1),\
                        (1, 0, 1, 1, 0, 1, 1, 0, 0, 1),\
                        (1, 1, 0, 1, 1, 0, 1, 0, 1, 0),\
                        (1, 1, 1, 0, 1, 1, 0, 1, 0, 0))
            sage: CA1 = CoveringArray(C1,2,2)
            sage: CA2 = CoveringArray(C2,2,2)
            sage: CA1>CA2
            False
        """
        if self.array_representation() > other.array_representation():
=======
        if self.ArrayRepresentation() > other.ArrayRepresentation():
>>>>>>> 5a4d10877a5c6414a20d23fcc217872d6287ba27
            return True
        else:
            return False
        
    def __ge__(self, other):
<<<<<<< HEAD
        r"""
        Return whether one covering array is greater than or 
        equal to another based on the lexicographic order on the rows
        
        EXAMPLES::
        
            sage: from sage.combinat.designs.covering_arrays import CoveringArray
            sage: C1 = ((1,1,1,0),\
                       (0,0,0,0),\
                       (1,0,1,1),\
                       (1,1,0,1),\
                       (0,1,1,1),)
            sage: C2 = ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0),\
                        (0, 0, 0, 0, 1, 1, 1, 1, 1, 1),\
                        (0, 1, 1, 1, 0, 0, 0, 1, 1, 1),\
                        (1, 0, 1, 1, 0, 1, 1, 0, 0, 1),\
                        (1, 1, 0, 1, 1, 0, 1, 0, 1, 0),\
                        (1, 1, 1, 0, 1, 1, 0, 1, 0, 0))
            sage: CA1 = CoveringArray(C1,2,2)
            sage: CA2 = CoveringArray(C2,2,2)
            sage: CA1>=CA2
            False
        """
        if self.array_representation() >= other.array_representation():
            return True
        else:
            return False
=======
        if self.ArrayRepresentation() >= other.ArrayRepresentation():
            return True
        else:
            return False
    
>>>>>>> 5a4d10877a5c6414a20d23fcc217872d6287ba27
