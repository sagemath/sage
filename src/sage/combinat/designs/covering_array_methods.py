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

def Kleitman_Specer_Katona(n): #Kleitman and Spencer [76] and Katona
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

def LFSR_for_strength3(n,p):
    r"""
    Returns a matrix that represents a CA()
    
    This method forms a sequence of entries of the ring GF(p) by using LFSR, 
    where the expression is given by the coeffients of the irreducible 
    polynomial of degree n with coefficients in GF(p). We form the circulant 
    matrix of the sequence and the reverse sequence, and remove the last 
    p^2+p+1 columns. The result is a covering array
    """
    
    R=PolynomialRing(GF(p),x)
    #find the irreducible polynomial degree n
    irr=R.irreducible_element(n)
    #we want to use coefficients of irr poly and solve for the coeff for the 
    #highest order term, this is done by taking the inverse of everything and
    #then dropping the first term
    a=[j:=-i for i in irr]
    a.pop()
    #for use of sages LFSR we set the fill to be all 0s except the first term
    fill=[0]*(n-1)
    fill.insert(0,1)
    #and the key is the expression from the coefficients found above
    key=a
    k=(p**n)-1 #the number of columns to be created
    LSFRseq = lfsr_sequence(key,fill,k)
    #We create a deep copy of this sequence so we can reverse it  
    LSFRrevseq = copy.deepcopy(L1) 
    LSFRrevseq.reverse() 
    #The array will be formed by the circ matrix of each sequence, plus the 
    #all zero row at the end
    M=block_matrix(3,[matrix.circulant(LSFRseq),matrix.circulant(LSFRrevseq)
                      ,matrix(GF(3),[0]*k)])
    #We then delete the excess columns to obtain the CA
    M=M.delete_columns([i for i in range (k-(p*p+p+1),k)])
    M2=[]
    for i in M:
        M2.append(tuple(i))
    return tuple(M2)
