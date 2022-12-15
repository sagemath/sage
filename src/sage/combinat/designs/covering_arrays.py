'''
An implementation of a Covering Array (CA) class for sagemath. 
The initial commit will include the definition of the Covering Array class 
and some basic methods to check and return the parameters n,k,v,t of the CA, 
as well as an ordering based on the lexicographic ordering of each row.

Later commits will include methods to create CAs from Perfect Hash Families 
and Covering Perfect Hash Families (CPHF) as well as recursive methods.

The Covering Array class may be used as a basis for an Orthogonal Array class
which will be implemented afterwards

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
                    finished=True
        
        #If t is given, make sure that t is the highest strength of the given
        #Covering Array, that is the Covering Array is not strength t+1
        else:
            finished=False
            while finished==False:
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
        return self.__array
    
    
    def __repr__(self):
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
            return True
        else:
            return False
        
    def __neq__(self, other):
        if self.ArrayRepresentation() != other.ArrayRepresentation():
            return True
        else:
            return False
        
    def __lt__(self, other):
        if self.ArrayRepresentation() < other.ArrayRepresentation():
            return True
        else:
            return False
        
    def __le__(self, other):
        if self.ArrayRepresentation() <= other.ArrayRepresentation():
            return True
        else:
            return False
        
    def __gt__(self, other):
        if self.ArrayRepresentation() > other.ArrayRepresentation():
            return True
        else:
            return False
        
    def __ge__(self, other):
        if self.ArrayRepresentation() >= other.ArrayRepresentation():
            return True
        else:
            return False
    