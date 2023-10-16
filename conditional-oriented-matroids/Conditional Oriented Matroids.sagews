#Conditional Oriented Matroids in SageMath
#Hery Randriamaro
#Institut für Mathematik, Universität Kassel
#hery.randriamaro@mathematik.uni-kassel.de


# 1 Conditional Oriented Matroids


# 1.1 Sign Systems


#Zero Set#
def Zero(X):
    z = Set()
    for i in range(len(X)):
        if X[i] == 0:
            z = z.union(Set([i]))
    return z


#Element Support#
def Support(X):
    s = Set()
    for i in range(len(X)):
        if X[i] in Set([-1, 1]):
            s = s.union(Set([i]))
    return s


#Separation Set#
def Separation(X, Y):
    s = Set()
    for i in range(len(X)):
        if X[i] in Set([-1, 1]):
            if X[i] == -Y[i]:
                s = s.union(Set([i]))
    return s


#Element Composition#
def sigma(a, b):
    if a == 0:
        return b
    else:
        return a

def Composition(X, Y):
    return tuple(sigma(X[i], Y[i]) for i in range(len(X)))


#Sign System Fiber#
def Restriction(X, A):
    return tuple(X[i] for i in Set(range(len(X))).difference(A))

def Fiber(L, X, A):
    F = Set()
    for Y in L:
        if Restriction(Y, A) == Restriction(X, A):
            F = F.union(Set([Y]))
    return F


# 1.2 Conditional Oriented Matroids


#Face Symmetry Condition#
def neg(X):
    return tuple(-X[i] for i in range(len(X)))

def FS(L):
    a = True
    for X in L:
        for Y in L:
            a = a & (Composition(X, neg(Y)) in L)
    return a


#Strong Elimination Condition#
def SE(L):
    a = True
    for X in L:
        for Y in L:
            S = Set()
            for Z in Fiber(L, Composition(X, Y), Separation(X, Y)):
                S = S.union(Zero(Z))
            a = a & (S.intersection(Separation(X, Y)) == Separation(X, Y))
    return a


#Conditional Oriented Matroid#
def COM(L):
    return FS(L) & SE(L)


# 1.3 Oriented Matroids


#Zero Vector Condition#
def Z(L):
    return (tuple(0 for i in range(len(L.random_element()))) in L)


#Oriented Matroid#
def OM(L):
    return COM(L) & Z(L)


# 1.4 Lopsided Systems


#Tope Composition Condition#
def TC(L):
    a = True
    for X in L:
        for Y in Tuples([-1, 1], len(L.random_element())):
            a = a & (Composition(X, Y) in L)
    return a


#Lopsided System#
def LS(L):
    return COM(L) & TC(L)


# 2 Construction of Conditional Oriented Matroids


# 2.1 Deletion and Contraction


#Conditional Oriented Matroid Deletion#
def Deletion(L, A):
    return Set([Restriction(X, A) for X in L])


#Conditional Oriented Matroid Contraction#
def Contraction(L, A):
    C = Set()
    for X in L:
        if Support(X).intersection(A) == Set():
            C = C.union(Set([Restriction(X, A)]))
    return C


# 2.2 Simplification


#Coloop Set#
def Coloop(L):
    C = Set()
    for e in range(len(L.random_element())):
        Xe = Set([X[e] for X in L])
        if len(Xe) == 1:
            C = C.union(Set([e]))
    return C


#Parallel Element Classes#
def parallel(L, e, f):
    Le = tuple(X[e] for X in L)
    Lf = tuple(X[f] for X in L)
    return (Le == Lf) | (Le == neg(Lf))

def Parallel_Element(L):
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


#Simplification#
def Simplification(L):
    P = Parallel_Element(L)
    F = Set()
    for Q in P:
        R = Q.difference(Set([Q.random_element()]))
        F = F.union(R)
    E = Set(range(len(L.random_element())))
    A = F.union(Coloop(L))
    return Deletion(L, A)


# 2.3 Generating from Topes


#Generalized Theorem of Mandel#
def Generalized_Mandel(T):
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


#Variables#
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


#Covector Distance Function#
def v(X, Y):
    x=1
    for i in Separation(X, Y):
        x = x*q(i, X[i])
    return x


#Varchenko Matrix#
def prec(X, Y):
    a = True
    for i in range(len(X)):
        a = a & (X[i] in Set([0, Y[i]]))
    return a

def IsMax(L, X):
    a = True
    if X not in L:
        a = False
    for Y in L.difference(Set([X])):
        a = a & (not prec(X, Y))
    return a

def Tope(L):
    S = Set()
    for X in L:
        if IsMax(L, X):
            S = S.union(Set([X]))
    return S

def V(L):
    return matrix([[v(X, Y) for Y in Tope(L)] for X in Tope(L)])


#Varchenko Determinant#
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

def Varchenko(L):
    return prod([(1-Weight(X))^(Theta(L, X)) for X in L.difference(Tope(L))])


# 3.2 The Aguiar-Mahajan Equation System


#Covector Face#
def Face(L, X):
    S = Set()
    for Y in L:
        if prec(X, Y):
            S = S.union(Set([Y]))
    return S


#Element Corank#
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


#Element Dimensional Rank#
def rk(L):
    k=0
    for X in L:
        k = max(k, crk(L, X))
    return k

def drk(L, X):
    return rk(L) - crk(L, X)


#Aguiar-Mahajan Equation System Solution#
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
        return (-1/(1 - v(X, neg(X))*v(neg(X), X)))*sum(AM(L, Y, o) +
         (-1)^(drk(L, X))*v(neg(X), X)*AM(L, neg(Y), o) for Y in S)

def AguiarMahajan(L, o):
    for X in L:
        print (X, " --> ", AM(L, X, o))