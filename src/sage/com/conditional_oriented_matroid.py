r"""
Conditional Oriented Matroids

Conditional oriented matroids are abstractions for diverse mathematical objects
like apartments of hyperplane arrangements, and partial cubes with gated antipodal
subgraphs. They are common generalizations of oriented matroids and lopsided
systems in particular. The following programs allow to manipulate of these objects.
Among other functions, the checks if a set corresponds to a conditional oriented
matroid, an oriented matroid, or a lopsided system are implemented. Functions
generating of a conditional oriented matroid from its tope set, computing the
Varchenko determinant of a conditional oriented matroid, and resolving the
Aguiar-Mahajan linear system of an oriented matroid are also programmed.

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


from sage.sets.set import Set, union, intersection, difference, random_element


# 1 Conditional Oriented Matroids


# 1.1 Sign Systems


r"""
A sign system is a pair `(E,\,\mathcal{L})` containing a finite set `E` and a
subset `\mathcal{L}` of `\{-1,\,0,\,1\}^E`. The elements of `\mathcal{L}` are
called sign vectors.
"""


def neg(X):
    r"""
    Return the negative `-X` of a sign vector `X`.
    
    INPUT:
    
    - ``X`` -- a sign vector
    
    .. MATH::
    
        -X := (-X_e\ |\ e \in E)
    
    EXAMPLES::
    
        sage: neg((0, 1, -1, -1, 0, 0, 1))
        (0, -1, 1, 1, 0, 0, -1)
    """
    return tuple(-X[i] for i in range(len(X)))
    

def zero(X):
    r"""
    Return the zero set `X^0` of a sign vector `X`.
    
    INPUT:
    
    - ``X`` -- a sign vector
    
    .. MATH::
    
        X^0 := \{e \in E\ |\ X_e = 0\}
    
    EXAMPLES::
    
        sage: zero((0, 1, -1, -1, 0, 0, 1))
        {0, 4, 5}
    """
    z = Set()
    for i in range(len(X)):
        if X[i] == 0:
            z = z.union(Set([i]))
    return z


def support(X):
    r"""
    Return the support `\underline{X}` of a sign vector `X`.
    
    INPUT:
    
    - ``X`` -- a sign vector
    
    .. MATH::
    
        \underline{X} := E \setminus X^0
        
    EXAMPLES::
    
        sage: support((0, 1, -1, -1, 0, 0, 1))
        {1, 2, 3, 6}
    """
    s = Set()
    for i in range(len(X)):
        if X[i] in Set([-1, 1]):
            s = s.union(Set([i]))
    return s


def separation(X, Y):
    r"""
    Return the separation set `\mathrm{S}(X,\,Y)` of two sign vectors `X` and `Y`.
    
    INPUT:
    
    - ``X`` -- a sign vector
    - ``Y`` -- a sign vector
    
    .. MATH::
    
        \mathrm{S}(X,\,Y) := \big\{e \in E\ \big|\ X_e = -Y_e \neq 0\big\}
        
    EXAMPLES::
    
        sage: separation((1, -1, 0, 0, 1), (0, 1, 0, 1, -1))
        {1, 4}
    """
    s = Set()
    for i in range(len(X)):
        if X[i] in Set([-1, 1]):
            if X[i] == -Y[i]:
                s = s.union(Set([i]))
    return s


def composition(X, Y):
    r"""
    Return the composition `X \circ Y` of the sign vectors `X` and `Y`.
    
    INPUT:
    
    - ``X`` -- a sign vector
    - ``Y`` -- a sign vector
    
    .. MATH::
    
        \forall e \in E,\,
        (X \circ Y)_e := \begin{cases} X_e & \text{if}\ X_e \neq 0, \\
        Y_e & \text{otherwise}.
        
    EXAMPLES::
    
        sage: composition((1, -1, 0, 0, 1), (0, 1, 0, 1, -1))
        (1, -1, 0, 1, 1)
    """
    def sigma(a, b):
        if a == 0:
            return b
        else:
            return a
    return tuple(sigma(X[i], Y[i]) for i in range(len(X)))
    
    
def prec(X, Y):
    r"""
    Return whether `X \preceq Y` for the partial order `\preceq`.
    
    INPUT:
    
    - ``X`` -- a sign vector
    - ``Y`` -- a sign vector
    
    .. MATH::
    
        \forall X, Y \in \mathcal{L}:\ X \preceq Y \ \Longleftrightarrow \ \forall e \in E,\, X_e \in \{0,\, Y_e\}
        
    EXAMPLES::
    
        sage: prec((1, -1, 0, 0, 1), (0, 1, 0, 1, -1))
        False
    """
    a = True
    for i in range(len(X)):
        a = a & (X[i] in Set([0, Y[i]]))
    return a
    
    
def is_max(L, X):
    r"""
    Return whether the sign vector `X` is maximal in this sign system `\mathcal{L}`.
    
    INPUT:
    
    - ``L`` -- a sign system
    - ``X`` -- a sign vector
    
    EXAMPLES::
    
        sage: L = Set([(1, -1, 0, 0, 1), (0, 1, 1, -1, 0), (1, 1, 1, -1, 0),
        ....: (0, 1, 0, 0, 1), (0, 0, 1, -1, 0), (0, 1, 0, 0, 0), (1, 1, 0,
        ....: -1, -1), (0, 1, 0, 1, -1)])
        sage: X = (1, -1, 0, 0, 1)
        sage: is_max(L, X)
        True
    """
    a = True
    if X not in L:
        a = False
    for Y in L.difference(Set([X])):
        a = a & (not prec(X, Y))
    return a


def face(L, X):
    r"""
    Return the face `\mathrm{F}(X)` of a sign vector `X` in a sign system `\mathcal{L}`.
    
    INPUT:
    
    - ``L`` -- a sign system
    - ``X`` -- a sign vector
    
    .. MATH::
    
        \mathrm{F}(X) := \{Y \in \mathcal{L}\ |\ X \preceq Y\}
        
    EXAMPLES::
    
        sage: L = Set([(1, -1, 0, 0, 1), (0, 1, 1, -1, 0), (1, 1, 1, -1, 0),
        ....: (0, 1, 0, 0, 1), (0, 0, 1, -1, 0), (0, 1, 0, 0, 0), (1, 1, 0,
        ....: -1, -1), (0, 1, 0, 1, -1)])
        sage: X = (0, 1, 0, 0, 0)
        sage: face(L, X)
        {(1, 1, 1, -1, 0), (0, 1, 0, 0, 0), (1, 1, 0, -1, -1), (0, 1, 0, 1, -1),
        (0, 1, 0, 0, 1), (0, 1, 1, -1, 0)}
    """
    S = Set()
    for Y in L:
        if prec(X, Y):
            S = S.union(Set([Y]))
    return S


def restriction(X, A):
    r"""
    Return the restriction `X \setminus A` of a sign vector `X` relative to a set `A`.
    
    INPUT:
    
    - ``X`` -- a sign vector
    - ``A`` -- a set
    
    .. MATH::
    
        X \setminus A \in \{-1,\,0,\,1\}^{E \setminus A},\, (X \setminus A)_e = X_e, \, \forall e \in E \setminus A
        
    EXAMPLES::
    
        sage: restriction((1, -1, 0, 0, 1), Set([1, 3, 4]))
        (1, 0)
    """
    return tuple(X[i] for i in Set(range(len(X))).difference(A))


def fiber(L, X, A):
    r"""
    Return the fiber `\mathrm{F}(X,\,A)` in a sign system `\mathcal{L}` relative to a sign vector `X` and a set `A`.
    
    INPUT:
    
    - ``L`` -- a sign system
    - ``X`` -- a sign vector
    - ``A`` -- a set
    
    .. MATH::
    
        \mathrm{F}(X,\,A) := \{Y \in \mathcal{L}\ |\ Y \setminus A = X \setminus A\}
        
    EXAMPLES::
    
        sage: L = Set([(1, -1, 0, 0, 1), (0, 1, 1, -1, 0), (1, 1, 1, -1, 0),
        ....: (0, 1, 0, 0, 1), (0, 0, 1, -1, 0), (0, 1, 0, 0, 0), (1, 1, 0,
        ....: -1, -1), (0, 1, 0, 1, -1)])
        sage: X = (1, -1, 0, 0, 1)
        sage: A = Set([1, 3, 4])
        sage: fiber(L, X, A)
        {(1, 1, 0, -1, -1), (1, -1, 0, 0, 1)}
    """
    F = Set()
    for Y in L:
        if restriction(Y, A) == restriction(X, A):
            F = F.union(Set([Y]))
    return F


# 1.2 Conditional Oriented Matroids


def com(L):
    r"""
    Return whether `\mathcal{L}` corresponds to a conditional oriented matroid.
    A conditional oriented matroid is a sign system `(E,\,\mathcal{L})` such that `\mathcal{L}` satisfies the following conditions:
    (FS) if `X,Y \in \mathcal{L}`, then `X \circ -Y \in \mathcal{L}`,
    (SE) for each pair `X,Y \in \mathcal{L}`, and every `e \in \mathrm{S}(X,\,Y)`, there exists `Z \in \mathcal{L}` such that
        `Z_e = 0` and `\forall f \in E \setminus \mathrm{S}(X,\,Y)`, `Z_f = (X \circ Y)_f = (Y \circ X)_f`.
    
    INPUT:
    
    - ``L`` -- a sign system
    
    EXAMPLES::
    
        sage: L1 = Set([(0, 1, 1, -1, 0), (1, 1, 1, -1, 0), (0, 1, 0, 0, 1), (0, 0,
        ....: 1, -1, 0), (0, 1, 0, 0, 0), (1, 1, 0, -1, -1), (0, 1, 0, 1, -1)])
        sage: L2 = Set([(1, 1, 1, 1, 0), (1, 0, 1, 1, 0), (1, -1, 1, 1, 0), (1, -1,
        ....: 0, 1, 0), (1, -1, 0, 0, 0), (1, -1, 1, 0, 0), (1, -1, -1, 1, 0), (1,
        ....: -1, -1, -1, 0), (1, -1, 1, -1, 0), (1, -1, -1, 0, 0), (1, -1, 0, -1,
        ....: 0)])
        sage: com(L1)
        False
        sage: com(L2)
        True
    """
    def fs(M):
        a = True
        for X in M:
            for Y in M:
                a = a & (composition(X, neg(Y)) in M)
        return a
    def se(N):
        a = True
        for X in N:
            for Y in N:
                S = Set()
                for Z in fiber(N, composition(X, Y), separation(X, Y)):
                    S = S.union(zero(Z))
                a = a & (S.intersection(separation(X, Y)) == separation(X, Y))
        return a
    return fs(L) & se(L)


# 1.3 Oriented Matroids


def om(L):
    r"""
    Return whether `\mathcal{L}` corresponds to an oriented matroid.
    An oriented matroid is a conditional oriented matroid `(E,\,\mathcal{L})` such that `\mathcal{L}` satisfies the following condition:
    (Z) the zero element `(0,\, \dots,\, 0)` belongs to `\mathcal{L}`.
    
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
        sage: om(L2)
        False
        sage: om(L3)
        True
    """
    def z(M):
        return (tuple(0 for i in range(len(M.random_element()))) in M)
    return com(L) & z(L)


# 1.4 Lopsided Systems


def ls(L):
    r"""
    Return whether `\mathcal{L}` corresponds to a lopsided system.
    A lopsided system is a conditional oriented matroid `(E,\,\mathcal{L})` such that `\mathcal{L}` satisfies the following condition:
    (TC) `\forall X \in \mathcal{L},\,\forall Y \in \{-1,\,1\}^E,\ X \circ Y \in \mathcal{L}`
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
        sage: ls(L3)
        False
        sage: ls(L4)
        True
    """
    def tc(M):
        a = True
        for X in M:
            for Y in Tuples([-1, 1], len(M.random_element())):
                a = a & (composition(X, Y) in M)
        return a
    return com(L) & tc(L)
    
    
def tope(L):
    r"""
    Return the tope set `\mathcal{T}` of a conditional oriented matroid `(E,\,\mathcal{L})`.
    
    INPUT:
    
    - ``L`` -- a conditional oriented matroid
    
    .. MATH::
    
        \mathcal{T} := \{X \in \mathcal{L}\ |\ \nexists Y \in \mathcal{L},\, X \prec Y\}
        
    EXAMPLES::
    
        sage: L2 = Set([(1, 1, 1, 1, 0), (1, 0, 1, 1, 0), (1, -1, 1, 1, 0), (1, -1,
        ....: 0, 1, 0), (1, -1, 0, 0, 0), (1, -1, 1, 0, 0), (1, -1, -1, 1, 0), (1,
        ....: -1, -1, -1, 0), (1, -1, 1, -1, 0), (1, -1, -1, 0, 0), (1, -1, 0, -1,
        ....: 0)])
        sage: tope(L2)
        {(1, -1, -1, -1, 0), (1, -1, -1, 1, 0), (1, -1, 1, 1, 0), (1, 1, 1, 1, 0),
        (1, -1, 1, -1, 0)}
    """
    S = Set()
    for X in L:
        if is_max(L, X):
            S = S.union(Set([X]))
    return S


# 2 Construction of Conditional Oriented Matroids


# 2.1 Deletion and Contraction


def deletion(L, A):
    r"""
    Return the set `\mathcal{L}\A` of the deletion `(E \setminus A,\, \mathcal{L} \setminus A)` of a conditional oriented matroid `(E,\, \mathcal{L})` relative to a set `A`.
    
    INPUT:
    
    - ``L`` -- a conditional oriented matroid
    - ``A`` -- a set
    
    .. MATH::
    
        `\mathcal{L} \setminus A = \{X \setminus A\ |\ X \in \mathcal{L}\}`
        
    EXAMPLES::
    
        sage: L2 = Set([(1, 1, 1, 1, 0), (1, 0, 1, 1, 0), (1, -1, 1, 1, 0), (1, -1,
        ....: 0, 1, 0), (1, -1, 0, 0, 0), (1, -1, 1, 0, 0), (1, -1, -1, 1, 0), (1,
        ....: -1, -1, -1, 0), (1, -1, 1, -1, 0), (1, -1, -1, 0, 0), (1, -1, 0, -1,
        ....: 0)])
        sage: deletion(L2, Set([2, 4]))
        {(1, 0, 1), (1, -1, 0), (1, -1, 1), (1, 1, 1), (1, -1, -1)}
    """
    return Set([restriction(X, A) for X in L])


def contraction(L, A):
    r"""
    Return the set `\mathcal{L}/A` of the contraction `(E \setminus A,\, \mathcal{L}/A)` of a conditional oriented matroid `(E,\, \mathcal{L})` relative to a set `A`.
    
    INPUT:
    
    - ``L`` -- a conditional oriented matroid
    - ``A`` -- a set
    
    .. MATH::
    
        `\mathcal{L}/A := \{X \setminus A\ |\ X \in \mathcal{L},\, \underline{X} \cap A = \varnothing\}`
        
    EXAMPLES::
    
        sage: L2 = Set([(1, 1, 1, 1, 0), (1, 0, 1, 1, 0), (1, -1, 1, 1, 0), (1, -1,
        ....: 0, 1, 0), (1, -1, 0, 0, 0), (1, -1, 1, 0, 0), (1, -1, -1, 1, 0), (1,
        ....: -1, -1, -1, 0), (1, -1, 1, -1, 0), (1, -1, -1, 0, 0), (1, -1, 0, -1,
        ....: 0)])
        sage: contraction(L2, Set([2, 4]))
        {(1, -1, 0), (1, -1, 1), (1, -1, -1)}
    """
    C = Set()
    for X in L:
        if support(X).intersection(A) == Set():
            C = C.union(Set([restriction(X, A)]))
    return C


# 2.2 Simplification


def coloop(L):
    r"""
    Return the coloop set for a conditional oriented matroid `(E,\, \mathcal{L})`.
    An element `e \in E` is a coloop if `\{X_e\ |\ X \in \mathcal{L}\} \subseteq \big\{\{-1\},\, \{0\},\, \{1\}\big\}`.
    
    INPUT:
    
    - ``L`` -- a conditional oriented matroid
        
    EXAMPLES::
    
        sage: L2 = Set([(1, 1, 1, 1, 0), (1, 0, 1, 1, 0), (1, -1, 1, 1, 0), (1, -1,
        ....: 0, 1, 0), (1, -1, 0, 0, 0), (1, -1, 1, 0, 0), (1, -1, -1, 1, 0), (1,
        ....: -1, -1, -1, 0), (1, -1, 1, -1, 0), (1, -1, -1, 0, 0), (1, -1, 0, -1,
        ....: 0)])
        sage: coloop(L2)
        {0, 4}
    """
    C = Set()
    for e in range(len(L.random_element())):
        Xe = Set([X[e] for X in L])
        if len(Xe) == 1:
            C = C.union(Set([e]))
    return C


def parallel_element(L):
    r"""
    Return the list of parallel elements for a conditional oriented matroid `(E,\, \mathcal{L})`.
    Two elements `e,f \in E` are parallel if either `X_e = X_f` for all `X \in \mathcal{L}`, or `X_e = -X_f` for all `X \in \mathcal{L}.
    
    INPUT:
    
    - ``L`` -- a conditional oriented matroid
        
    EXAMPLES::
    
        sage: L3 = Set([(0, 0, 0, 0, 0), (0, 0, -1, -1, 0), (0, 1, -1, -1, 0), (0,
        ....: 1, 0, -1, 0), (0, 1, 1, -1, 0), (0, 1, 1, 0, 0), (0, 1, 1, 1, 0), (0,
        ....: 0, 1, 1, 0), (0, -1, 1, 1, 0), (0, -1, 0, 1, 0), (0, -1, -1, 1, 0),
        ....: (0, -1, -1, 0, 0), (0, -1, -1, -1, 0)])
        sage: parallel_element(L3)
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


def simplification(L):
    r"""
    Return a simplification of a conditional oriented matroid `(E,\, \mathcal{L})`.
    A conditional oriented matroid is simple if it has neither coloops nor distinct parallel elements.
    A homomorphism between two conditional oriented matroids `(E,\, \mathcal{L})` and `(F,\, \mathcal{M})` is a function `h: \mathcal{L} \rightarrow \mathcal{M}` such that `\forall X, Y \in \mathcal{L},\ X \preceq Y \Longrightarrow h(X) \preceq h(Y)`.
    One says that the homomorphism `h` is an isomorphism if it is additionally bijective.
    A simplification of `(E,\, \mathcal{L})` is a simple conditional oriented matroid which is isomorphic to `(E,\, \mathcal{L})`.
    
    INPUT:
    
    - ``L`` -- a conditional oriented matroid
    
    OUTPUT: the simplification of L
        
    EXAMPLES::
    
        sage: L4 = Set([(1, -1, 1, 0, 0), (1, -1, 1, -1, 0), (1, -1, 1, 1, 0), (1,
        ....: -1, 1, 0, -1), (1, -1, 1, 0, 1), (1, -1, 1, -1, -1), (1, -1, 1, -1,
        ....: 1), (1, -1, 1, 1, -1), (1, -1, 1, 1, 1)])
        sage: simplification(L4)
        {(0, 1), (-1, -1), (0, 0), (-1, 1), (1, 1), (1, -1), (-1, 0), (1, 0),
        (0, -1)}
        
    REFERENCES:
    
    For more information, see Proposition 2.3 of [Ran2024]_.
    """
    P = parallel_element(L)
    F = Set()
    for Q in P:
        R = Q.difference(Set([Q.random_element()]))
        F = F.union(R)
    E = Set(range(len(L.random_element())))
    A = F.union(coloop(L))
    return deletion(L, A)


# 2.3 Generating from Topes


def generalized_mandel(T):
    r"""
    Return the set `\mathcal{L}` associated to a conditional oriented matroid `(E,\, \mathcal{L})` generated by its tope set `\mathcal{T}`.
    The algorithm is based on the formula `\mathcal{L} = \big\{X \in \{-1,\,0,\,1\}^E\ \big|\ \forall T \in \mathcal{T},\, X \circ -T \in \mathcal{T}\big\}`.
    
    INPUT:
    
    - ``T`` -- a tope set
    
    OUTPUT: the conditional oriented matroid generated by T
        
    EXAMPLES::
    
        sage: T = Set([(-1, -1, -1, -1), (1, -1, -1, -1), (1, 1, -1, -1), (1, -1, 1
        ....: , -1), (1, 1, 1, -1), (1, 1, 1, 1)])
        sage: generalized_mandel(T)
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
            a = a & (composition(X, neg(Y)) in T)
        if a == True:
            M = M.union(Set([tuple(X)]))
    return M


# 3 Special Functions


# 3.1 The Varchenko Matrix


def v(X, Y):
    r"""
    Return the Aguiar-Mahajan distance between two topes.
    
    INPUT:
    
    - ``X`` -- a tope
    - ``Y`` -- a tope
    
    .. MATH::
    
        `\mathrm{v}(X,\,Y) := \prod_{e \in \mathrm{S}(X,\,Y)} q_{e,X_e}` with `q_{e,X_e} = ae` if `X_e = -1` and `q_{e,X_e} = be` otherwise
        
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
    for i in separation(X, Y):
        x = x*q(i, X[i])
    return x


def v_matrix(L):
    r"""
    Return the Varchenko matrix of a conditional oriented matroid.
    
    INPUT:
    
    - ``L`` -- a conditional oriented matroid
    
    OUTPUT: the Varchenko matrix of the simplification of L
    
    .. MATH::
    
        `\mathrm{V}(\mathcal{L}) := \big(\mathrm{v}(U,\,T)\big)_{T,U \in \mathcal{T}}`
        
    EXAMPLES::
    
        sage: L2 = Set([(1, 1, 1, 1, 0), (1, 0, 1, 1, 0), (1, -1, 1, 1, 0), (1, -1,
        ....: 0, 1, 0), (1, -1, 0, 0, 0), (1, -1, 1, 0, 0), (1, -1, -1, 1, 0), (1,
        ....: -1, -1, -1, 0), (1, -1, 1, -1, 0), (1, -1, -1, 0, 0), (1, -1, 0, -1,
        ....: 0)])
        sage: v_matrix(L2)
        [       1       a2       b1    a0*a2    a2*b1]
        [      b2        1    b1*b2       a0       b1]
        [      a1    a1*a2        1 a0*a1*a2       a2]
        [   b0*b2       b0 b0*b1*b2        1    b0*b1]
        [   a1*b2       a1       b2    a0*a1        1]
        
    .. SEEALSO::
        
        :mod:`sage.geometry.hyperplane_arrangement.hyperplane`.
    """
    M = tope(simplification(L))
    return matrix([[v(X, Y) for Y in M] for X in M])


def varchenko(L):
    r"""
    Return the Varchenko determinant of a conditional oriented matroid.
    
    INPUT:
    
    - ``L`` -- a conditional oriented matroid
    
    OUTPUT: the Varchenko determinant of the simplification of L
    
    .. MATH::
    
        `\det \mathrm{V}(\mathcal{L}) = \prod_{X \in \mathcal{L} \setminus \mathcal{T}} \Big(1 - \prod_{e \in X^0} q_{e,-1} q_{e,1}\Big)^{\theta(X)}`
        where `\theta(X) = \frac{\#\big\{T \in \mathcal{T}\ \big|\ \mathrm{Max}\,\{Y \in \mathcal{L}\ |\ Y \prec T,\, Y_f = 0\} = \{X\}\big\}}{2}` and `f \in X^0`.
        
    EXAMPLES::
    
        sage: L2 = Set([(1, 1, 1, 1, 0), (1, 0, 1, 1, 0), (1, -1, 1, 1, 0), (1, -1,
        ....: 0, 1, 0), (1, -1, 0, 0, 0), (1, -1, 1, 0, 0), (1, -1, -1, 1, 0), (1,
        ....: -1, -1, -1, 0), (1, -1, 1, -1, 0), (1, -1, -1, 0, 0), (1, -1, 0, -1,
        ....: 0)])
        sage: varchenko(L2)
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
        for i in zero(X):
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
        i = zero(X).random_element()
        for Y in tope(L):
            if is_max(iBoundary(L, i, Y), X):
                M.append(Y)
        return len(M)/2
    M = simplification(L)
    return prod([(1-Weight(X))^(Theta(M, X)) for X in M.difference(tope(M))])


# 3.2 The Aguiar-Mahajan Equation System


def aguiar_mahajan(L, o):
    r"""
    Return the solution of an Aguiar-Mahajan linear system.
    Let `(E,\, \mathcal{L})` be an oriented matroid. Assign a variable `x_X` to each element `X \in \mathcal{L}`.
    The Aguiar-Mahajan system for `(E,\,\mathcal{L})` is the linear equation system `\sum_{\substack{X \in \mathcal{L} \\ Y \circ X = Y}} x_X\,\mathrm{v}(X,\,Y) = 0` indexed by `Y \in \mathcal{L} \setminus (0,\dots,0)`.
    It has a unique solution which can be computed recursively with the formula `x_Y = \frac{-1}{1 - \mathrm{v}(Y,\,-Y) \, \mathrm{v}(-Y,\,Y)} \sum_{\substack{X \in \mathcal{L} \\ X \prec Y}} \big(x_X + (-1)^{\mathrm{drk}\,Y}x_{-X} \, \mathrm{v}(-Y,\,Y)\big` with `\mathrm{drk}\,Y := \mathrm{rank}\,\mathcal{L} - \mathrm{corank}\,Y`.
    
    INPUT:
    
    - ``L`` -- an oriented matroid
    - ``o`` -- an initial value associated to the zero vector
    
    OUTPUT: it gives X --> solution corresponding to X, for each element X of L
    
    EXAMPLES::
    
        sage: L3 = Set([(0, 0, 0, 0, 0), (0, 0, -1, -1, 0), (0, 1, -1, -1, 0), (0,
        ....: 1, 0, -1, 0), (0, 1, 1, -1, 0), (0, 1, 1, 0, 0), (0, 1, 1, 1, 0), (0,
        ....: 0, 1, 1, 0), (0, -1, 1, 1, 0), (0, -1, 0, 1, 0), (0, -1, -1, 1, 0),
        ....: (0, -1, -1, 0, 0), (0, -1, -1, -1, 0)])
        sage: aguiar_mahajan(L3, 1)
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
    def min_system(L):
        X = L.random_element()
        for Y in L:
            if prec(Y, X):
                X=Y
        return X
    def crk(L, X):
        k=0
        F = face(L, X).difference(Set([X]))
        while not (F == Set()):
            k = k+1
            Y = min_system(F)
            F = face(L, Y).difference(Set([Y]))
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
    def am(L, X, o):
        S = Inf(L, X)
        if X == min_system(L):
            return o
        else:
            return (-1/(1-v(X, neg(X))*v(neg(X), X))) * sum(am(L, Y, o)+
            (-1)^(drk(L, X))*v(neg(X), X)*am(L, neg(Y), o) for Y in S)
    M = simplification(L)
    for X in M:
        print (X, " --> ", am(M, X, o))
