
"""
Covering Arrays

A Covering Array, denoted CA(N,k,v,t), is an n by k array with entries from a 
set of v elements with theproperty that in every selection of t columns, each 
row contains every sequence of t-elements at least once.

An Orthogonal Array, denoted OA(N,k,v,t) is a covering array with the property
that each row contains every sequence of t-elements exactly once

REFERENCES:
    
- reference 1

- reference 2

AUTHORS:
    
- Aaron Dwyer and brett stevens

.. NOTES::
    
This is a work in progress, it will be an implementation of a Covering Array (CA) 
class for sagemath. The initial commit will include the definition of the 
Covering Array class and some basic methods to check and return the parameters 
N,k,v,t of the CA, as well as an ordering based on the lexicographic ordering 
of each row.

Later commits will include methods to create CAs from Perfect Hash Families 
and Covering Perfect Hash Families (CPHF) as well as recursive methods.

The Covering Array class may be used as a basis for an Orthogonal Array class
which will be implemented afterwards

"""

import itertools
import copy

def is_covering_array(Array,levels,strength): 
    """
    Input is a tuple of tuples called 'Array' that represents a N x k array, 
    and integers ``levels`` and ``strength`` which represent the desired 
    covering arrays values for v and t respectively.
    
    Checks if 'Array' with given parameters is a Covering Array CA(n,k,v,t) 
    and outputs True or False.
    
    EXAMPLES::
        
        sage: C = ((1,1,1,0),
                   (1,1,0,1),
                   (1,0,1,1),
                   (0,1,1,1),
                   (0,0,0,0))        
        sage: is_covering_array(C,2,2)
        True
        sage: is_covering_array(C,2,3)
        False
        
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
    """
    A class for covering array objects. The objects will contain the array 
    itself and its paramters v, and t
    
    INPUT:
        
    -``Array`` -- The N x k array itself. This is stored as a tuple of tuples.
    the N and k parameters are derived from this inputted array
    
    -``strength`` -- The parameter t, such that in any selection of t columns
    of the array, every t tuple appears at least once. 
    
    If a covering array has strength t, it also has strength t-1,t-2...1. So 
    we wish to only store the maximum strength of the given covering array, so 
    the class automatically checks if the CA(N,K,v,t) is a CA(N,k,v,t+1) until
    it finds the maximum t
    
    -----------------all arrays are CA t=1??------
    
    -``levels`` -- The paramter v, such that v is the size of the symbol set 
    of the array.
    
    If no such symbol set or size is given, then a v will be assumed by 
    counting how many unique elements appear in the array.
    
    EXAMPLES::
        
        sage: C = ((1,1,1,0),
                   (1,1,0,1),
                   (1,0,1,1),
                   (0,1,1,1),
                   (0,0,0,0))           
        sage: CoveringArray(C,strength=2)
        A 5 by 4 Covering Array of strength 2 with 2 levels
        
        sage: C= ((0,0,0,0,0,0,0,0,0,0),
                 (1,1,1,1,1,1,1,1,1,1),
                 (1,1,1,0,1,0,0,0,0,1),
                 (1,0,1,1,0,1,0,1,0,0),
                 (1,0,0,0,1,1,1,0,0,0),
                 (0,1,1,0,0,1,0,0,1,0),
                 (0,0,1,0,1,0,1,1,1,0),
                 (1,1,0,1,0,0,1,0,1,0),
                 (0,0,0,1,1,1,0,0,1,1),
                 (0,0,1,1,0,0,1,0,0,1),
                 (0,1,0,1,1,0,0,1,0,0),
                 (1,0,0,0,0,0,0,1,1,1),
                 (0,1,0,0,0,1,1,1,0,1))
        sage: CoveringArray(C,2,3)
        A 5 by 4 Covering Array of strength 3 with 2 levels
    """
    def __init__(self, Array, strength=None, levels=None):
        """
        Constructor function    
        
        EXAMPLES::
            
            sage: C = ((1,1,1,0),
                       (1,1,0,1),
                       (1,0,1,1),
                       (0,1,1,1),
                       (0,0,0,0))           
            sage: CoveringArray(C,strength=2)
            A 5 by 4 Covering Array of strength 2 with 2 levels
            
            sage: C= ((0,0,0,0,0,0,0,0,0,0),
                     (1,1,1,1,1,1,1,1,1,1),
                     (1,1,1,0,1,0,0,0,0,1),
                     (1,0,1,1,0,1,0,1,0,0),
                     (1,0,0,0,1,1,1,0,0,0),
                     (0,1,1,0,0,1,0,0,1,0),
                     (0,0,1,0,1,0,1,1,1,0),
                     (1,1,0,1,0,0,1,0,1,0),
                     (0,0,0,1,1,1,0,0,1,1),
                     (0,0,1,1,0,0,1,0,0,1),
                     (0,1,0,1,1,0,0,1,0,0),
                     (1,0,0,0,0,0,0,1,1,1),
                     (0,1,0,0,0,1,1,1,0,1))
            sage: CoveringArray(C,2,3)
            A 5 by 4 Covering Array of strength 3 with 2 levels
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
        

        #If v is not inputted then the symbol set may be assumed from what 
        #symbols are in the array by flattening the Array and counting the 
        #number of unique entries.
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
                    finished=True
        
        #If t is given, make sure that t is the highest strength of the given
        #Covering Array, that is the Covering Array is not strength t+1
        else:
            finished=False
            while finished==False:
                if is_covering_array(Array, levels ,strength):
                    strength+=1
                else:
                    strength-=1
                    finished=True
        self.__t=strength
        
    def numrows(self):
        """
        Method that returns the number of rows, N, of the covering array
        
        EXAMPLES::
                
            sage: C = ((1,1,1,0),
                       (1,1,0,1),
                       (1,0,1,1),
                       (0,1,1,1),
                       (0,0,0,0))         
            sage: CA=CoveringArray(C,2,2)
            sage: CA.numrows()
            5
        """
        return self.__n
        
    def numcols(self):
        """
        Method that returns the number of columns, k, of the covering array
        
        EXAMPLES::
                
            sage: C = ((1,1,1,0),
                       (1,1,0,1),
                       (1,0,1,1),
                       (0,1,1,1),
                       (0,0,0,0))          
            sage: CA=CoveringArray(C,2,2)
            sage: CA.numcols()
            4
        """
        return self.__k
    
    def levels(self):
        """
        Method that returns the number of levels for the covering array, which 
        is the paramter v, such that v is the size of the symbol set of the 
        array.
    
        EXAMPLES::
                
            sage: C = ((1,1,1,0),
                       (1,1,0,1),
                       (1,0,1,1),
                       (0,1,1,1),
                       (0,0,0,0))         
            sage: CA=CoveringArray(C,2,2)
            sage: CA.levels()
            2
        """
        return self.__v
        
    def strength(self):
        """
        Method that returns the strength of the covering array, which is the 
        paramter t, such that in any selection of t columns of the array, 
        every t tuple appears at least once.
    
        EXAMPLES::
                
            sage: C = ((1,1,1,0),
                       (1,1,0,1),
                       (1,0,1,1),
                       (0,1,1,1),
                       (0,0,0,0))        
            sage: CA=CoveringArray(C,2,2)
            sage: CA.strength()
            2
        """
        return self.__t
    
    def array_representation(self):
        """
        Method that returns the covering array as a tuple of tuples, but where
        the output is such that each row of the array is sorted in 
        lexicographic order
        
        EXAMPLES::
                
            sage: C = ((1,1,1,0),
                       (0,0,0,0),
                       (1,0,1,1),
                       (1,1,0,1),
                       (0,1,1,1),)         
            sage: CA=CoveringArray(C,2,2)
            sage: CA.array_representation()
            ((0,0,0,0),(0,1,1,1),(1,0,1,1),(1,1,0,1),(1,1,1,0))
            
        """
        return self.__array
    
    
    def __repr__(self):
        """
        EXAMPLES::
            
            sage: C = ((1,1,1,0),
                       (1,1,0,1),
                       (1,0,1,1),
                       (0,1,1,1),
                       (0,0,0,0))           
            sage: CoveringArray(C,strength=2)
            A 5 by 4 Covering Array of strength 2 with 2 levels
            
            sage: C= ((0,0,0,0,0,0,0,0,0,0),
                     (1,1,1,1,1,1,1,1,1,1),
                     (1,1,1,0,1,0,0,0,0,1),
                     (1,0,1,1,0,1,0,1,0,0),
                     (1,0,0,0,1,1,1,0,0,0),
                     (0,1,1,0,0,1,0,0,1,0),
                     (0,0,1,0,1,0,1,1,1,0),
                     (1,1,0,1,0,0,1,0,1,0),
                     (0,0,0,1,1,1,0,0,1,1),
                     (0,0,1,1,0,0,1,0,0,1),
                     (0,1,0,1,1,0,0,1,0,0),
                     (1,0,0,0,0,0,0,1,1,1),
                     (0,1,0,0,0,1,1,1,0,1))
            sage: CoveringArray(C,2,3)
            A 5 by 4 Covering Array of strength 3 with 2 levels
        """
        return 'A {} by {} Covering Array of strength {} with {} levels'.format(
            self.nrows, self.ncols, self.strength, self.levels)
        
    __str__=__repr__
    
    
    def __eq__(self, other):
        """
        A method that desrcibes whether two covering arrays are equal by 
        considering the array with rows sorted in lexicographic order
        
        EXMAPLES::
            sage: C1 = ((1,1,1,0),
                       (0,0,0,0),
                       (1,0,1,1),
                       (1,1,0,1),
                       (0,1,1,1),)
            sage: C2 = ((1,1,1,0),
                       (1,1,0,1),
                       (1,0,1,1),
                       (0,1,1,1),
                       (0,0,0,0))
            sage: C3 = ((0,0,0,0,0,0,0,0,0,0),
                       (1,1,1,1,1,1,1,1,1,1),
                       (1,1,1,0,1,0,0,0,0,1),
                       (1,0,1,1,0,1,0,1,0,0),
                       (1,0,0,0,1,1,1,0,0,0),
                       (0,1,1,0,0,1,0,0,1,0),
                       (0,0,1,0,1,0,1,1,1,0),
                       (1,1,0,1,0,0,1,0,1,0),
                       (0,0,0,1,1,1,0,0,1,1),
                       (0,0,1,1,0,0,1,0,0,1),
                       (0,1,0,1,1,0,0,1,0,0),
                       (1,0,0,0,0,0,0,1,1,1),
                       (0,1,0,0,0,1,1,1,0,1))
            sage: C1=C2
            True
            sage: C1=C3
            False
        """
        if self.ArrayRepresentation() == other.ArrayRepresentation():
            return True
        else:
            return False
        
    def __neq__(self, other):
        """
        A method that desrcibes whether two covering arrays are not equal by 
        considering the array with rows sorted in lexicographic order
        
        EXMAPLES::
            sage: C1 = ((1,1,1,0),
                       (0,0,0,0),
                       (1,0,1,1),
                       (1,1,0,1),
                       (0,1,1,1),)
            sage: C2 = ((1,1,1,0),
                       (1,1,0,1),
                       (1,0,1,1),
                       (0,1,1,1),
                       (0,0,0,0))
            sage: C3 = ((0,0,0,0,0,0,0,0,0,0),
                       (1,1,1,1,1,1,1,1,1,1),
                       (1,1,1,0,1,0,0,0,0,1),
                       (1,0,1,1,0,1,0,1,0,0),
                       (1,0,0,0,1,1,1,0,0,0),
                       (0,1,1,0,0,1,0,0,1,0),
                       (0,0,1,0,1,0,1,1,1,0),
                       (1,1,0,1,0,0,1,0,1,0),
                       (0,0,0,1,1,1,0,0,1,1),
                       (0,0,1,1,0,0,1,0,0,1),
                       (0,1,0,1,1,0,0,1,0,0),
                       (1,0,0,0,0,0,0,1,1,1),
                       (0,1,0,0,0,1,1,1,0,1))
            sage: C1!=C2
            False
            sage: C1!=C3
            True
        """
        if self.ArrayRepresentation() != other.ArrayRepresentation():
            return True
        else:
            return False
        
    def __lt__(self, other):
        """
        A method that desrcibes whether one covering array is less than another
        based on the lexicographic order on the rows
        
        EXMAPLES::
            sage: C1 = ((1,1,1,0),
                       (0,0,0,0),
                       (1,0,1,1),
                       (1,1,0,1),
                       (0,1,1,1),)
            
        """
        if self.ArrayRepresentation() < other.ArrayRepresentation():
            return True
        else:
            return False
        
    def __le__(self, other):
        """
        A method that desrcibes whether one covering array is less than or 
        equal to another based on the lexicographic order on the rows
        
        EXMAPLES::
            sage: C1 = ((1,1,1,0),
                       (0,0,0,0),
                       (1,0,1,1),
                       (1,1,0,1),
                       (0,1,1,1),)
            
        """
        if self.ArrayRepresentation() <= other.ArrayRepresentation():
            return True
        else:
            return False
        
    def __gt__(self, other):
        """
        A method that desrcibes whether one covering array is greater than 
        another based on the lexicographic order on the rows
        
        EXMAPLES::
            sage: C1 = ((1,1,1,0),
                       (0,0,0,0),
                       (1,0,1,1),
                       (1,1,0,1),
                       (0,1,1,1),)
            
        """
        if self.ArrayRepresentation() > other.ArrayRepresentation():
            return True
        else:
            return False
        
        
        (0,1,0,0)
        
        (1,0)
        
    def __ge__(self, other):
        """
        A method that desrcibes whether one covering array is greater than or 
        equal to another based on the lexicographic order on the rows
        
        EXMAPLES::
            sage: C1 = ((1,1,1,0),
                       (0,0,0,0),
                       (1,0,1,1),
                       (1,1,0,1),
                       (0,1,1,1),)
            
        """
        if self.ArrayRepresentation() >= other.ArrayRepresentation():
            return True
        else:
            return False
    