import itertools
import math

def tuple_weight(t):
    '''
    given a tuple t returns the integer sum of all of its values
    '''
    ctr=0
    for i in t:
        ctr+=i
    return ctr

def Strength2(n):
    '''
    Given an integer n, returns a binary covering array with strength 2
    To do this we form a matrix represented as a tuple of tuples where the 
    columns are all binary n-tuples of weight ceiling(n/2) that have 0 
    in 1st position
    
    Reference is Combinatorial aspects...
    '''
    resulttupleT=()
    x = [0,1]
    for p in itertools.product(x, repeat=n-1):
        if tuple_weight(p)==math.ceil(n/2):
            column=(0,)+p
            resulttupleT+=(column,)
    M=Matrix(resulttupleT)
    M1=M.transpose()
    return(tuple(M1))

