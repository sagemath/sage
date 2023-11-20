r"""
Conditional Oriented Matroids

Conditional oriented matroids are abstractions
for diverse mathematical objects like apartments of hyperplane arrangements also
called realizable COMs, and partial cubes with gated antipodal subgraphs. They are
common generalizations of oriented matroids and lopsided systems in particular. The
following programs allow to manipulate of these objects. Among other functions, the
checks if a set corresponds to a conditional oriented matroid, an oriented matroid,
or a lopsided system are implemented. Functions generating of a conditional oriented
matroid from its tope set, computing the Varchenko determinant of a conditional
oriented matroid, and resolving the Aguiar-Mahajan linear system of an oriented
matroid are also programmed.

REFERENCES: For more information on conditional oriented matroids, see [BCK2018]_.

AUTHOR:

- Hery Randriamaro (2023-11-20): initial version

"""

# ****************************************************************************
#       Copyright (C) 2023 Hery Randriamaro <hery.randriamaro@mathematik.uni-kassel.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************




# 1 Conditional Oriented Matroids


# 1.1 Sign Systems


def neg(X):
    """
    Return the negative of a sign vector.
    
    INPUT:
    
    - ``X`` -- a sign vector
    
    EXAMPLES::
    
        sage: neg((0, 1, -1, -1, 0, 0, 1))
        (0, -1, 1, 1, 0, 0, -1)
    """
    return tuple(-X[i] for i in range(len(X)))
    

def Zero(X):
    """
    Return the zero set of a sign vector.
    
    INPUT:
    
    - ``X`` -- a sign vector
    
    EXAMPLES::
    
        sage: Zero((0, 1, -1, -1, 0, 0, 1))
        {0, 4, 5}
    """
    z = Set()
    for i in range(len(X)):
        if X[i] == 0:
            z = z.union(Set([i]))
    return z


def Support(X):
    """
    Return the support of a sign vector.
    
    INPUT:
    
    - ``X`` -- a sign vector
    
    EXAMPLES::
    
        sage: Support((0, 1, -1, -1, 0, 0, 1))
        {1, 2, 3, 6}
    """
    s = Set()
    for i in range(len(X)):
        if X[i] in Set([-1, 1]):
            s = s.union(Set([i]))
    return s


def Separation(X, Y):
    """
    Return the separation set of two sign vectors.
    
    INPUT:
    
    - ``X`` -- a sign vector
    - ``Y`` -- a sign vector
    
    EXAMPLES::
    
        sage: Separation((1, -1, 0, 0, 1), (0, 1, 0, 1, -1))
        {1, 4}
    """
    s = Set()
    for i in range(len(X)):
        if X[i] in Set([-1, 1]):
            if X[i] == -Y[i]:
                s = s.union(Set([i]))
    return s


def Composition(X, Y):
    """
    Return the composition of the first sign vector with the second.
    
    INPUT:
    
    - ``X`` -- a sign vector
    - ``Y`` -- a sign vector
    
    EXAMPLES::
    
        sage: Composition((1, -1, 0, 0, 1), (0, 1, 0, 1, -1))
        (1, -1, 0, 1, 1)
    """
    def sigma(a, b):
        if a == 0:
            return b
        else:
            return a
    return tuple(sigma(X[i], Y[i]) for i in range(len(X)))
    
    
def prec(X, Y):
    """
    Return whether the first sign vector is smaller that the second.
    
    INPUT:
    
    - ``X`` -- a sign vector
    - ``Y`` -- a sign vector
    
    EXAMPLES::
    
        sage: prec((1, -1, 0, 0, 1), (0, 1, 0, 1, -1))
        False
    """
    a = True
    for i in range(len(X)):
        a = a & (X[i] in Set([0, Y[i]]))
    return a
    
    
def IsMax(L, X):
    """
    Return whether this sign vector is maximal in this sign system.
    
    INPUT:
    
    - ``L`` -- a sign system
    - ``X`` -- a sign vector
    
    EXAMPLES::
    
        sage: L = Set([(1, -1, 0, 0, 1), (0, 1, 1, -1, 0), (1, 1, 1, -1, 0),
        ....: (0, 1, 0, 0, 1), (0, 0, 1, -1, 0), (0, 1, 0, 0, 0), (1, 1, 0,
        ....: -1, -1), (0, 1, 0, 1, -1)])
        sage: X = (1, -1, 0, 0, 1)
        sage: IsMax(L, X)
        True
    """
    a = True
    if X not in L:
        a = False
    for Y in L.difference(Set([X])):
        a = a & (not prec(X, Y))
    return a


def Face(L, X):
    """
    Return the face of a sign vector in a sign system.
    
    INPUT:
    
    - ``L`` -- a sign system
    - ``X`` -- a sign vector
    
    EXAMPLES::
    
        sage: L = Set([(1, -1, 0, 0, 1), (0, 1, 1, -1, 0), (1, 1, 1, -1, 0),
        ....: (0, 1, 0, 0, 1), (0, 0, 1, -1, 0), (0, 1, 0, 0, 0), (1, 1, 0,
        ....: -1, -1), (0, 1, 0, 1, -1)])
        sage: X = (0, 1, 0, 0, 0)
        sage: Face(L, X)
        {(1, 1, 1, -1, 0), (0, 1, 0, 0, 0), (1, 1, 0, -1, -1), (0, 1, 0, 1, -1),
        (0, 1, 0, 0, 1), (0, 1, 1, -1, 0)}
    """
    S = Set()
    for Y in L:
        if prec(X, Y):
            S = S.union(Set([Y]))
    return S


def Restriction(X, A):
    """
    Return the restriction of a sign vector relative to a set.
    
    INPUT:
    
    - ``X`` -- a sign vector
    - ``A`` -- a set
    
    EXAMPLES::
    
        sage: Restriction((1, -1, 0, 0, 1), Set([1, 3, 4]))
        (1, 0)
    """
    return tuple(X[i] for i in Set(range(len(X))).difference(A))


def Fiber(L, X, A):
    """
    Return the fiber in a sign system relative to a sign vector and a set.
    
    INPUT:
    
    - ``L`` -- a sign system
    - ``X`` -- a sign vector
    - ``A`` -- a set
    
    EXAMPLES::
    
        sage: L = Set([(1, -1, 0, 0, 1), (0, 1, 1, -1, 0), (1, 1, 1, -1, 0),
        ....: (0, 1, 0, 0, 1), (0, 0, 1, -1, 0), (0, 1, 0, 0, 0), (1, 1, 0,
        ....: -1, -1), (0, 1, 0, 1, -1)])
        sage: X = (1, -1, 0, 0, 1)
        sage: A = Set([1, 3, 4])
        sage: Fiber(L, X, A)
        {(1, 1, 0, -1, -1), (1, -1, 0, 0, 1)}
    """
    F = Set()
    for Y in L:
        if Restriction(Y, A) == Restriction(X, A):
            F = F.union(Set([Y]))
    return F


# 1.2 Conditional Oriented Matroids


def COM(L):
    """
    Return whether this sign system is a conditional oriented matroid.
    
    INPUT:
    
    - ``L`` -- a sign system
    
    EXAMPLES::
    
        sage: L1 = Set([(0, 1, 1, -1, 0), (1, 1, 1, -1, 0), (0, 1, 0, 0, 1), (0, 0,
        ....: 1, -1, 0), (0, 1, 0, 0, 0), (1, 1, 0, -1, -1), (0, 1, 0, 1, -1)])
        sage: L2 = Set([(1, 1, 1, 1, 0), (1, 0, 1, 1, 0), (1, -1, 1, 1, 0), (1, -1,
        ....: 0, 1, 0), (1, -1, 0, 0, 0), (1, -1, 1, 0, 0), (1, -1, -1, 1, 0), (1,
        ....: -1, -1, -1, 0), (1, -1, 1, -1, 0), (1, -1, -1, 0, 0), (1, -1, 0, -1,
        ....: 0)])
        sage: COM(L1)
        False
        sage: COM(L2)
        True
    """
    def FS(M):
        a = True
        for X in M:
            for Y in M:
                a = a & (Composition(X, neg(Y)) in M)
        return a
    def SE(N):
        a = True
        for X in N:
            for Y in N:
                S = Set()
                for Z in Fiber(N, Composition(X, Y), Separation(X, Y)):
                    S = S.union(Zero(Z))
                a = a & (S.intersection(Separation(X, Y)) == Separation(X, Y))
        return a
    return FS(L) & SE(L)


# 1.3 Oriented Matroids


def OM(L):
    """
    Return whether this sign system is an oriented matroid.
    
    INPUT:
    
    - ``L`` -- a sign system
    
    EXAMPLES::
    
        sage: L2 = Set([(1, 1, 1, 1, 0), (1, 0, 1, 1, 0), (1, -1, 1, 1, 0), (1, -1,
        ....: 0, 1, 0), (1, -1, 0, 0, 0), (1, -1, 1, 0, 0), (1, -1, -1, 1, 0), (1,
        ....: -1, -1, -1, 0), (1, -1, 1, -1, 0), (1, -1, -1, 0, 0), (1, -1, 0, -1,
        ....: 0)])
        sage: L3 = Set([(0, 0, 0, 0, 0), (0, 0, -1, -1, 0), (0, 1, -1, -1, 0), (0,
        ....: 1, 0, -1, 0), (0, 1, 1, -1, 0), (0, 1, 1, 0, 0), (0, 1, 1, 1, 0), (0,
        ....: 0, 1, 1, 0), (0, -1, 1, 1, 0), (0, -1, 0, 1, 0), (0, -1, -1, 1, 0),
        ....: (0, -1, -1, 0, 0), (0, -1, -1, -1, 0)])
        sage: OM(L2)
        False
        sage: OM(L3)
        True
    """
    def Z(M):
        return (tuple(0 for i in range(len(M.random_element()))) in M)
    return COM(L) & Z(L)


# 1.4 Lopsided Systems


def LS(L):
    """
    Return whether this sign system is a lopsided system.
    
    INPUT:
    
    - ``L`` -- a sign system
    
    EXAMPLES::
    
        sage: L3 = Set([(0, 0, 0, 0, 0), (0, 0, -1, -1, 0), (0, 1, -1, -1, 0), (0,
        ....: 1, 0, -1, 0), (0, 1, 1, -1, 0), (0, 1, 1, 0, 0), (0, 1, 1, 1, 0), (0,
        ....: 0, 1, 1, 0), (0, -1, 1, 1, 0), (0, -1, 0, 1, 0), (0, -1, -1, 1, 0),
        ....: (0, -1, -1, 0, 0), (0, -1, -1, -1, 0)])
        sage: L4 = Set([(1, -1, 1, 0, 0), (1, -1, 1, -1, 0), (1, -1, 1, 1, 0), (1,
        ....: -1, 1, 0, -1), (1, -1, 1, 0, 1), (1, -1, 1, -1, -1), (1, -1, 1, -1,
        ....: 1), (1, -1, 1, 1, -1), (1, -1, 1, 1, 1)])
        sage: LS(L3)
        False
        sage: LS(L4)
        True
    """
    def TC(M):
        a = True
        for X in M:
            for Y in Tuples([-1, 1], len(M.random_element())):
                a = a & (Composition(X, Y) in M)
        return a
    return COM(L) & TC(L)
    
    
def Tope(L):
    """
    Return the tope set of a conditional oriented matroid.
    
    INPUT:
    
    - ``L`` -- a conditional oriented matroid
    
    EXAMPLES::
    
        sage: L2 = Set([(1, 1, 1, 1, 0), (1, 0, 1, 1, 0), (1, -1, 1, 1, 0), (1, -1,
        ....: 0, 1, 0), (1, -1, 0, 0, 0), (1, -1, 1, 0, 0), (1, -1, -1, 1, 0), (1,
        ....: -1, -1, -1, 0), (1, -1, 1, -1, 0), (1, -1, -1, 0, 0), (1, -1, 0, -1,
        ....: 0)])
        sage: Tope(L2)
        {(1, -1, -1, -1, 0), (1, -1, -1, 1, 0), (1, -1, 1, 1, 0), (1, 1, 1, 1, 0),
        (1, -1, 1, -1, 0)}
    """
    S = Set()
    for X in L:
        if IsMax(L, X):
            S = S.union(Set([X]))
    return S


# 2 Construction of Conditional Oriented Matroids


# 2.1 Deletion and Contraction


def Deletion(L, A):
    """
    Return the deletion of a conditional oriented matroid relative to a set.
    
    INPUT:
    
    - ``L`` -- a conditional oriented matroid
    - ``A`` -- a set
        
    EXAMPLES::
    
        sage: L2 = Set([(1, 1, 1, 1, 0), (1, 0, 1, 1, 0), (1, -1, 1, 1, 0), (1, -1,
        ....: 0, 1, 0), (1, -1, 0, 0, 0), (1, -1, 1, 0, 0), (1, -1, -1, 1, 0), (1,
        ....: -1, -1, -1, 0), (1, -1, 1, -1, 0), (1, -1, -1, 0, 0), (1, -1, 0, -1,
        ....: 0)])
        sage: Deletion(L2, Set([2, 4]))
        {(1, 0, 1), (1, -1, 0), (1, -1, 1), (1, 1, 1), (1, -1, -1)}
    """
    return Set([Restriction(X, A) for X in L])


def Contraction(L, A):
    """
    Return the deletion of a conditional oriented matroid relative to a set.
    
    INPUT:
    
    - ``L`` -- a conditional oriented matroid
    - ``A`` -- a set
        
    EXAMPLES::
    
        sage: L2 = Set([(1, 1, 1, 1, 0), (1, 0, 1, 1, 0), (1, -1, 1, 1, 0), (1, -1,
        ....: 0, 1, 0), (1, -1, 0, 0, 0), (1, -1, 1, 0, 0), (1, -1, -1, 1, 0), (1,
        ....: -1, -1, -1, 0), (1, -1, 1, -1, 0), (1, -1, -1, 0, 0), (1, -1, 0, -1,
        ....: 0)])
        sage: Contraction(L2, Set([2, 4]))
        {(1, -1, 0), (1, -1, 1), (1, -1, -1)}
    """
    C = Set()
    for X in L:
        if Support(X).intersection(A) == Set():
            C = C.union(Set([Restriction(X, A)]))
    return C


# 2.2 Simplification


def Coloop(L):
    """
    Return the coloop set for a conditional oriented matroid.
    
    INPUT:
    
    - ``L`` -- a conditional oriented matroid
        
    EXAMPLES::
    
        sage: L2 = Set([(1, 1, 1, 1, 0), (1, 0, 1, 1, 0), (1, -1, 1, 1, 0), (1, -1,
        ....: 0, 1, 0), (1, -1, 0, 0, 0), (1, -1, 1, 0, 0), (1, -1, -1, 1, 0), (1,
        ....: -1, -1, -1, 0), (1, -1, 1, -1, 0), (1, -1, -1, 0, 0), (1, -1, 0, -1,
        ....: 0)])
        sage: Coloop(L2)
        {0, 4}
    """
    C = Set()
    for e in range(len(L.random_element())):
        Xe = Set([X[e] for X in L])
        if len(Xe) == 1:
            C = C.union(Set([e]))
    return C


def Parallel_Element(L):
    """
    Return the list of parallel elements for a conditional oriented matroid.
    
    INPUT:
    
    - ``L`` -- a conditional oriented matroid
        
    EXAMPLES::
    
        sage: L3 = Set([(0, 0, 0, 0, 0), (0, 0, -1, -1, 0), (0, 1, -1, -1, 0), (0,
        ....: 1, 0, -1, 0), (0, 1, 1, -1, 0), (0, 1, 1, 0, 0), (0, 1, 1, 1, 0), (0,
        ....: 0, 1, 1, 0), (0, -1, 1, 1, 0), (0, -1, 0, 1, 0), (0, -1, -1, 1, 0),
        ....: (0, -1, -1, 0, 0), (0, -1, -1, -1, 0)])
        sage: Parallel_Element(L3)
        [{0, 4}, {1}, {2}, {3}]
    """
    def parallel(M, e, f):
        Le = tuple(X[e] for X in M)
        Lf = tuple(X[f] for X in M)
        return (Le==Lf) | (Le==neg(Lf))
    E = Set(range(len(L.random_element())))
    P = []
    while E.cardinality() > 0:
        for e in E:
            Pe = Set([])
            for f in E:
                if parallel(L, e, f):
                    Pe = Pe.union(Set([f]))
            if Pe.cardinality() > 0:
                P.append(Pe)
            E = E.difference(Pe)
    return P


def Simplification(L):
    """
    Return the simplification of a conditional oriented matroid.
    
    INPUT:
    
    - ``L`` -- a conditional oriented matroid
    
    OUTPUT: the simplification of L
        
    EXAMPLES::
    
        sage: L4 = Set([(1, -1, 1, 0, 0), (1, -1, 1, -1, 0), (1, -1, 1, 1, 0), (1,
        ....: -1, 1, 0, -1), (1, -1, 1, 0, 1), (1, -1, 1, -1, -1), (1, -1, 1, -1,
        ....: 1), (1, -1, 1, 1, -1), (1, -1, 1, 1, 1)])
        sage: Simplification(L4)
        {(0, 1), (-1, -1), (0, 0), (-1, 1), (1, 1), (1, -1), (-1, 0), (1, 0),
        (0, -1)}
    """
    P = Parallel_Element(L)
    F = Set()
    for Q in P:
        R = Q.difference(Set([Q.random_element()]))
        F = F.union(R)
    E = Set(range(len(L.random_element())))
    A = F.union(Coloop(L))
    return Deletion(L, A)


# 2.3 Generating from Topes


def Generalized_Mandel(T):
    """
    Return the conditional oriented matroid generated by a tope set.
    
    It implements Algorithm 1 of the following article: H. Randriamaro,
    A Generalization of a Theorem of Mandel, Eng. Math. Lett. (2023), Article ID 1.
    
    INPUT:
    
    - ``T`` -- a tope set
    
    OUTPUT: the conditional oriented matroid generated by T
        
    EXAMPLES::
    
        sage: T = Set([(-1, -1, -1, -1), (1, -1, -1, -1), (1, 1, -1, -1), (1, -1, 1
        ....: , -1), (1, 1, 1, -1), (1, 1, 1, 1)])
        sage: Generalized_Mandel(T)
        {(1, -1, -1, -1), (1, 0, 1, -1), (0, -1, -1, -1), (1, 0, -1, -1), (-1, -1,
        -1, -1), (1, -1, 1, -1), (1, 1, 1, 0), (1, 0, 0, -1), (1, 1, 0, -1), (1, -1
        , 0, -1), (1, 1, 1, -1), (1, 1, 1, 1), (1, 1, -1, -1)}
        
    REFERENCES:
    
    For more information on its algorithm, see [Ran2023]_.
    """
    M = Set()
    for X in Tuples([-1, 0, 1], len(T.random_element())):
        a = True
        for Y in T:
            a = a & (Composition(X, neg(Y)) in T)
        if a == True:
            M = M.union(Set([tuple(X)]))
    return M


# 3 Special Functions


# 3.1 The Varchenko Matrix


def v(X, Y):
    """
    Return the Aguiar-Mahajan distance between two topes.
    
    INPUT:
    
    - ``X`` -- a tope
    - ``Y`` -- a tope
    
    EXAMPLES::
    
        sage: v((1, -1, -1, 1, 1), (-1, -1, 1, 1, -1))
        a2*b0*b4
    """
    class VariableGenerator(object):
        def __init__(self, prefix):
            self.__prefix = prefix
        @cached_method
        def __getitem__(self, key):
            return SR.var("%s%s"%(self.__prefix,key))
    a = VariableGenerator('a')
    b = VariableGenerator('b')
    def q(k, l):
        if l == -1:
            return(a[k])
        if l == 1:
            return(b[k])
    x=1
    for i in Separation(X, Y):
        x = x*q(i, X[i])
    return x


def V(L):
    """
    Return the Varchenko matrix of a conditional oriented matroid.
    
    INPUT:
    
    - ``L`` -- a conditional oriented matroid
    
    OUTPUT: the Varchenko matrix of the simplification of L
        
    EXAMPLES::
    
        sage: L2 = Set([(1, 1, 1, 1, 0), (1, 0, 1, 1, 0), (1, -1, 1, 1, 0), (1, -1,
        ....: 0, 1, 0), (1, -1, 0, 0, 0), (1, -1, 1, 0, 0), (1, -1, -1, 1, 0), (1,
        ....: -1, -1, -1, 0), (1, -1, 1, -1, 0), (1, -1, -1, 0, 0), (1, -1, 0, -1,
        ....: 0)])
        sage: V(L2)
        [       1       a2       b1    a0*a2    a2*b1]
        [      b2        1    b1*b2       a0       b1]
        [      a1    a1*a2        1 a0*a1*a2       a2]
        [   b0*b2       b0 b0*b1*b2        1    b0*b1]
        [   a1*b2       a1       b2    a0*a1        1]
        
    .. SEEALSO::
        
        :mod:`sage.geometry.hyperplane_arrangement.hyperplane`.
    """
    M = Tope(Simplification(L))
    return matrix([[v(X, Y) for Y in M] for X in M])


def Varchenko(L):
    """
    Return the Varchenko determinant of a conditional oriented matroid.
    
    INPUT:
    
    - ``L`` -- a conditional oriented matroid
    
    OUTPUT: the Varchenko determinant of the simplification of L
        
    EXAMPLES::
    
        sage: L2 = Set([(1, 1, 1, 1, 0), (1, 0, 1, 1, 0), (1, -1, 1, 1, 0), (1, -1,
        ....: 0, 1, 0), (1, -1, 0, 0, 0), (1, -1, 1, 0, 0), (1, -1, -1, 1, 0), (1,
        ....: -1, -1, -1, 0), (1, -1, 1, -1, 0), (1, -1, -1, 0, 0), (1, -1, 0, -1,
        ....: 0)])
        sage: Varchenko(L2)
        -(a0*b0 - 1)*(a1*b1 - 1)^2*(a2*b2 - 1)^2
        
    REFERENCES:
    
    For more information on this function, see Theorem 4.36 of [Ran2024]_.
    
    For more information on Varchenko determinants, see the following references:
    
    - [AM2017]_
    
    - [Ran2022]_
    
    - [Var1993]_
    """
    def Weight(X):
        x=1
        for i in Zero(X):
            x = x*a[i]*b[i]
        return x
    def iBoundary(L, i, X):
        S = Set()
        for Y in L.difference(Set([X])):
            if prec(Y, X) and (Y[i] == 0):
                S = S.union(Set([Y]))
        return S
    def Theta(L, X):
        M = []
        i = Zero(X).random_element()
        for Y in Tope(L):
            if IsMax(iBoundary(L, i, Y), X):
                M.append(Y)
        return len(M)/2
    M = Simplification(L)
    return prod([(1-Weight(X))^(Theta(M, X)) for X in M.difference(Tope(M))])


# 3.2 The Aguiar-Mahajan Equation System


def AguiarMahajan(L, o):
    """
    Return the solution of an Aguiar-Mahajan linear system.
    
    INPUT:
    
    - ``L`` -- an oriented matroid
    - ``o`` -- an initial value associated to the zero vector
    
    OUTPUT: it gives X --> solution corresponding to X, for each covector X of L
    
    EXAMPLES::
    
        sage: L3 = Set([(0, 0, 0, 0, 0), (0, 0, -1, -1, 0), (0, 1, -1, -1, 0), (0,
        ....: 1, 0, -1, 0), (0, 1, 1, -1, 0), (0, 1, 1, 0, 0), (0, 1, 1, 1, 0), (0,
        ....: 0, 1, 1, 0), (0, -1, 1, 1, 0), (0, -1, 0, 1, 0), (0, -1, -1, 1, 0),
        ....: (0, -1, -1, 0, 0), (0, -1, -1, -1, 0)])
        sage: AguiarMahajan(L3, 1)
        (0, -1, -1)  -->  -(b1*b2 - 1)/(a1*a2*b1*b2 - 1)
        (-1, 0, 1)  -->  -(a2*b0 - 1)/(a0*a2*b0*b2 - 1)
        (1, 0, -1)  -->  -(a0*b2 - 1)/(a0*a2*b0*b2 - 1)
        (1, 1, 0)  -->  -(a0*a1 - 1)/(a0*a1*b0*b1 - 1)
        (-1, 1, 1)  -->  -((a0*b2 - 1)*a1*a2*b0/(a0*a2*b0*b2 - 1) + (b1*b2 -
        1)*a1*a2*b0/(a1*a2*b1*b2 - 1) - a1*a2*b0 + (a1*a2 - 1)/(a1*a2*b1*b2 - 1) +
        (a2*b0 - 1)/(a0*a2*b0*b2 - 1) - 1)/(a0*a1*a2*b0*b1*b2 - 1)
        (-1, -1, -1)  -->  -((a0*a1 - 1)*b0*b1*b2/(a0*a1*b0*b1 - 1) + (a1*a2 -
        1)*b0*b1*b2/(a1*a2*b1*b2 - 1) - b0*b1*b2 + (b0*b1 - 1)/(a0*a1*b0*b1 - 1) +
        (b1*b2 - 1)/(a1*a2*b1*b2 - 1) - 1)/(a0*a1*a2*b0*b1*b2 - 1)
        (0, 0, 0)  -->  1
        (-1, -1, 1)  -->  -((a0*a1 - 1)*a2*b0*b1/(a0*a1*b0*b1 - 1) + (a0*b2 -
        1)*a2*b0*b1/(a0*a2*b0*b2 - 1) - a2*b0*b1 + (a2*b0 - 1)/(a0*a2*b0*b2 - 1) +
        (b0*b1 - 1)/(a0*a1*b0*b1 - 1) - 1)/(a0*a1*a2*b0*b1*b2 - 1)
        (1, -1, -1)  -->  -((a1*a2 - 1)*a0*b1*b2/(a1*a2*b1*b2 - 1) + (a2*b0 -
        1)*a0*b1*b2/(a0*a2*b0*b2 - 1) - a0*b1*b2 + (a0*b2 - 1)/(a0*a2*b0*b2 - 1) +
        (b1*b2 - 1)/(a1*a2*b1*b2 - 1) - 1)/(a0*a1*a2*b0*b1*b2 - 1)
        (1, 1, -1)  -->  -((a2*b0 - 1)*a0*a1*b2/(a0*a2*b0*b2 - 1) + (b0*b1 -
        1)*a0*a1*b2/(a0*a1*b0*b1 - 1) - a0*a1*b2 + (a0*a1 - 1)/(a0*a1*b0*b1 - 1) +
        (a0*b2 - 1)/(a0*a2*b0*b2 - 1) - 1)/(a0*a1*a2*b0*b1*b2 - 1)
        (-1, -1, 0)  -->  -(b0*b1 - 1)/(a0*a1*b0*b1 - 1)
        (1, 1, 1)  -->  -((b0*b1 - 1)*a0*a1*a2/(a0*a1*b0*b1 - 1) + (b1*b2 -
        1)*a0*a1*a2/(a1*a2*b1*b2 - 1) - a0*a1*a2 + (a0*a1 - 1)/(a0*a1*b0*b1 - 1) +
        (a1*a2 - 1)/(a1*a2*b1*b2 - 1) - 1)/(a0*a1*a2*b0*b1*b2 - 1)
        (0, 1, 1)  -->  -(a1*a2 - 1)/(a1*a2*b1*b2 - 1)
        
    REFERENCES:
    
    For more information on this function, see Theorem 4.44 of [Ran2024]_.
    
    For more information on Aguiar-Mahajan systems, see the following references:
    
    - [AM2017]_
    
    - [Ran2022]_
    """
    def Min(L):
        X = L.random_element()
        for Y in L:
            if prec(Y, X):
                X=Y
        return X
    def crk(L, X):
        k=0
        F = Face(L, X).difference(Set([X]))
        while not (F == Set()):
            k = k+1
            Y = Min(F)
            F = Face(L, Y).difference(Set([Y]))
        return k
    def rk(L):
        k=0
        for X in L:
            k = max(k, crk(L, X))
        return k
    def drk(L, X):
        return rk(L) - crk(L, X)
    def Inf(L, X):
        S = Set()
        for Y in L.difference(Set([X])):
            if prec(Y, X):
                S = S.union(Set([Y]))
        return S
    def AM(L, X, o):
        S = Inf(L, X)
        if X == Min(L):
            return o
        else:
            return (-1/(1-v(X, neg(X))*v(neg(X), X))) * sum(AM(L, Y, o)+
            (-1)^(drk(L, X))*v(neg(X), X)*AM(L, neg(Y), o) for Y in S)
    M = Simplification(L)
    for X in M:
        print (X, " --> ", AM(M, X, o))
