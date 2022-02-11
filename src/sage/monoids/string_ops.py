"Utility functions on strings"

#*****************************************************************************
#       Copyright (C) 2007 David Kohel <kohel@maths.usyd.edu.au>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import RealField
from .string_monoid_element import StringMonoidElement

def strip_encoding(S):
    """
    The upper case string of S stripped of all non-alphabetic characters.

    EXAMPLES::

        sage: S = "The cat in the hat."
        sage: strip_encoding(S)
        'THECATINTHEHAT'
    """
    if not isinstance(S,str):
        raise TypeError("Argument S (= %s) must be a string.")
    X = ''
    for i in range(len(S)):
        C = S[i]
        if C.isalpha():
            X += S[i].upper()
    return X

def frequency_distribution(S, n=1, field=None):
    """
    The probability space of frequencies of n-character substrings of S.
    """
    if isinstance(S,tuple):
        S = list(S)
    elif isinstance(S,(str,StringMonoidElement)):
        S = [ S[i:i+n] for i in range(len(S)-n+1) ]
    if field is None:
        field = RealField()
    if isinstance(S,list):
        P = {}
        N = len(S)
        eps = field(1)/N
        for i in range(N):
            c = S[i]
            if c in P:
                P[c] += eps
            else:
                P[c] = eps
        from sage.probability.random_variable import DiscreteProbabilitySpace
        return DiscreteProbabilitySpace(S,P,field)
    raise TypeError("Argument S (= %s) must be a string, list, or tuple.")

def coincidence_index(S,n=1):
    """
    The coincidence index of the string S.

    EXAMPLES::

        sage: S = strip_encoding("The cat in the hat.")
        sage: coincidence_index(S)
        0.120879120879121
    """
    if not isinstance(S,str):
        try:
            S.coincidence_index(n)
        except AttributeError:
            raise TypeError("Argument S (= %s) must be a string.")
    S = strip_encoding(S)
    N = len(S)-n+1
    X = {}
    for i in range(N):
        c = S[i:i+n]
        if c in X:
            X[c] += 1
        else:
            X[c] = 1
    RR = RealField()
    return RR(sum([ m*(m-1) for m in X.values() ]))/RR(N*(N-1))

def coincidence_discriminant(S,n=2):
    """
    Input: A tuple of strings, e.g. produced as decimation of transposition
    ciphertext, or a sample plaintext.
    Output: A measure of the difference of probability of association of
    character pairs, relative to their independent one-character probabilities.

    EXAMPLES::

        sage: S = strip_encoding("The cat in the hat.")
        sage: coincidence_discriminant([ S[i:i+2] for i in range(len(S)-1) ])
        0.0827001855677322
    """
    if not isinstance(S,(list,tuple)):
        raise TypeError("Argument S (= %s) must be a list or tuple" % S)
    if n != 2:
        raise ValueError("Argument n (= %s) is only implemented for n = 2" % n)
    truth = True
    for bool in ( isinstance(c,(str,StringMonoidElement)) for c in S ):
        truth = truth and bool
    if not truth:
        raise TypeError("Argument S (= %s) must be a list of strings.")
    for bool in ( len(c) == n for c in S ):
        truth = truth and bool
    if not truth:
        raise ValueError("Argument S (= %s) must be a list of strings of length 2" % S)
    X1 = [ frequency_distribution([ s[i] for s in S]) for i in range(2) ]
    XX = frequency_distribution(S)
    if isinstance(S[0],StringMonoidElement):
        M = S[0].parent()
        n = M.ngens()
        return sum([ (XX(M([i,j]))-X1[0](M([i]))*X1[1](M([j])))**2 for i in range(n) for j in range(n) ])
    AZ = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    return sum([ (XX(AZ[i]+AZ[j])-X1[0](AZ[i])*X1[1](AZ[j]))**2 for i in range(26) for j in range(26) ])
