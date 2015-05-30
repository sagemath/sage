## -*- encoding: utf-8 -*-
"""
Doctests from French Sage book
Test file for chapter "Programmation" ("Programming")

Tests extracted from ./programmation.tex.

Sage example in ./programmation.tex, line 67::

  sage: 2*3; 3*4; 4*5          # un commentaire, 3 résultats
  6
  12
  20

Sage example in ./programmation.tex, line 78::

  sage: 123 + \
  ....: 345
  468

Sage example in ./programmation.tex, line 113::

  sage: import keyword; keyword.kwlist
  ['and', 'as', 'assert', 'break', 'class', 'continue', 'def', 'del',
  'elif', 'else', 'except', 'exec', 'finally', 'for', 'from',
  'global', 'if', 'import', 'in', 'is', 'lambda', 'not', 'or', 'pass',
  'print', 'raise', 'return', 'try', 'while', 'with', 'yield']

Sage example in ./programmation.tex, line 189::

  sage: y = 3; y = 3 * y + 1; y = 3 * y + 1; y
  31

Sage example in ./programmation.tex, line 205::

  sage: a, b = 10, 20 # (a, b) = (10, 20) et [10, 20] possibles
  sage: a, b = b, a
  sage: a, b
  (20, 10)

Sage example in ./programmation.tex, line 214::

  sage: temp = a; a = b; b = temp # est équivalent à : a, b = b, a

Sage example in ./programmation.tex, line 221::

  sage: x, y = var('x, y'); a = x ; b = y
  sage: a, b
  (x, y)
  sage: a = a + b ; b = a - b ; a = a - b
  sage: a, b
  (y, x)

Sage example in ./programmation.tex, line 237::

  sage: 2 + 2 == 2^2, 3 * 3 == 3^3
  (True, False)

Sage example in ./programmation.tex, line 276::

  sage: for k in [1..5]:
  ....:    print 7*k  # bloc qui contient une seule instruction
  7
  14
  21
  28
  35

Sage example in ./programmation.tex, line 372::

  sage: S = 0 ; k = 0         #        La somme S commence à 0
  sage: while e^k <= 10^6:    #            e^13 <= 10^6 < e^14
  ....:     S = S + k^2       #           ajout des carrés k^2
  ....:     k = k + 1

Sage example in ./programmation.tex, line 378::

  sage: S
  819

Sage example in ./programmation.tex, line 403::

  sage: x = 10^4; u = 1; n = 0          # invariant : u = 2^n
  sage: while u <= x: n = n+1; u = 2*u  # ou n += 1; u *= 2
  sage: n
  14

Sage example in ./programmation.tex, line 449::

  sage: U = 1.0              # ou U = 1. ou U = 1.000
  sage: for n in [1..20]:
  ....:   U = 1 / (1 + U^2)
  sage: U
  0.682360434761105

Sage example in ./programmation.tex, line 492::

  sage: S = 0 ; n = 10
  sage: for k in [1..n]:
  ....:     S = S + (2*k) * (2*k+1)
  sage: S
  1650

Sage example in ./programmation.tex, line 502::

  sage: n, k = var('n, k') ; res = sum(2*k*(2*k+1), k, 1, n)
  sage: res, factor(res)    # résultat développé puis factorisé
  (4/3*n^3 + 3*n^2 + 5/3*n, 1/3*(4*n + 5)*(n + 1)*n)

Sage example in ./programmation.tex, line 574::

  sage: U = 2.0; V = 50.0;
  sage: while V-U >= 1.0e-6:      # 1.0e-6 signifie 1.0*10^-6
  ....:   temp = U
  ....:   U = 2 * U * V / (U + V)
  ....:   V = (temp + V) / 2
  sage: U, V
  (9.99999999989..., 10.0000000001...)

Sage example in ./programmation.tex, line 635::

  sage: U = 0.0  # la somme S0 est vide, de valeur nulle
  sage: V = -1.0 # S1 = -1/1^3
  sage: n = 0      # U et V contiennent S(2n) et S(2n+1)
  sage: while U-V >= 1.0e-6:
  ....:   n = n+1             # n += 1 est équivalent
  ....:   U = V + 1/(2*n)^3   # passage de S(2n-1) à S(2n)
  ....:   V = U - 1/(2*n+1)^3 # passage de S(2n) à S(2n+1)
  sage: V, U
  (-0.901543155458595, -0.901542184868447)

Sage example in ./programmation.tex, line 807::

  sage: u = 6 ; n = 0
  sage: while u != 1:    # test "différent de"
  ....:   if u % 2 == 0: # l'opérateur % donne le reste euclidien
  ....:     u = u//2     # // : quotient de la division euclidienne
  ....:   else:
  ....:     u = 3*u+1
  ....:   n = n+1
  sage: n
  8

Sage example in ./programmation.tex, line 883::

  sage: def fct2 (x, y):
  ....:   return x^2 + y^2
  sage: a = var('a')
  sage: fct2 (a, 2*a)
  5*a^2

Sage example in ./programmation.tex, line 908::

  sage: def essai (u):
  ....:   t = u^2
  ....:   return t*(t+1)
  sage: t = 1 ; u = 2
  sage: essai(3), t, u
  (90, 1, 2)

Sage example in ./programmation.tex, line 918::

  sage: a = b = 1
  sage: def f(): global a; a = b = 2
  sage: f(); a, b
  (2, 1)

Sage example in ./programmation.tex, line 928::

  sage: def MoyAH (u, v):
  ....:    u, v = min(u, v), max(u, v)
  ....:    while v-u > 2.0e-8:
  ....:       u, v = 2*u*v/(u+v), (u+v)/2
  ....:    return (u+v) / 2

Sage example in ./programmation.tex, line 935::

  sage: MoyAH (1., 2.)
  1.41421...
  sage: MoyAH                       # correspond à une fonction
    <function MoyAH at ...>

Sage example in ./programmation.tex, line 990::

  sage: def fact1 (n):
  ....:    res = 1
  ....:    for k in [1..n]: res = res*k
  ....:    return res

Sage example in ./programmation.tex, line 996::

  sage: def fact2 (n):
  ....:    if n == 0: return 1
  ....:    else: return n*fact2(n-1)

Sage example in ./programmation.tex, line 1013::

  sage: def fib1 (n):
  ....:   if n == 0 or n == 1: return n
  ....:   else:
  ....:     U = 0 ; V = 1 # les termes initiaux u0 et u1
  ....:     for k in [2..n]: W = U+V ; U = V ; V = W
  ....:     return V
  sage: fib1(8)
  21

Sage example in ./programmation.tex, line 1037::

  sage: def fib2 (n):
  ....:   if 0 <= n <= 1: return n      # pour n = 0 ou n = 1
  ....:   else: return fib2(n-1) + fib2(n-2)

Sage example in ./programmation.tex, line 1079::

  sage: a = 2; n = 6; res = 1      # 1 est le neutre du produit
  sage: for k in [1..n]: res = res*a
  sage: res                  # La valeur de res est 2^6
  64

Sage example in ./programmation.tex, line 1150::

  sage: def puiss1 (a, n):
  ....:   if n == 0: return 1
  ....:   elif n % 2 == 0: b = puiss1 (a, n//2); return b*b
  ....:   else: return a * puiss1(a, n-1)

Sage example in ./programmation.tex, line 1156::

  sage: puiss1 (2, 11)               # a pour résultat 2^11
  2048

Sage example in ./programmation.tex, line 1177::

  sage: def puiss2 (u, k):
  ....:   v = 1
  ....:   while k != 0:
  ....:     if k % 2 == 0: u = u*u ; k = k//2
  ....:     else: v = v*u ; k = k-1
  ....:   return v

Sage example in ./programmation.tex, line 1185::

  sage: puiss2 (2, 10)               # a pour résultat 2^10
  1024

Sage example in ./programmation.tex, line 1238::

  sage: def fib3 (n):
  ....:   A = matrix ([[0, 1], [1, 1]]) ; X0 = vector ([0, 1])
  ....:   return (A^n*X0)[0]

Sage example in ./programmation.tex, line 1243::

  sage: def fib4 (n):
  ....:   return (matrix([[0,1], [1,1]])^n * vector([0,1]))[0]

Sage example in ./programmation.tex, line 1257::

  sage: print 2^2, 3^3, 4^4 ; print 5^5, 6^6
  4 27 256
  3125 46656

Sage example in ./programmation.tex, line 1265::

  sage: for k in [1..10]: print '+', k,
  + 1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9 + 10

Sage example in ./programmation.tex, line 1273::

  sage: print 10, 0.5 ; print 10+0.5 ; print 10.0, 5
  10 0.500000000000000
  10.5000000000000
  10.0000000000000 5
  sage: print 10+0, 5 ; print str(10)+str(0.5)
  10 5
  100.500000000000000

Sage example in ./programmation.tex, line 1294::

  sage: for k in [1..6]: print '%2d^4 = %4d' % (k, k^4)
   1^4 =    1
   2^4 =   16
   3^4 =   81
   4^4 =  256
   5^4 =  625
   6^4 = 1296

Sage example in ./programmation.tex, line 1344::

  sage: L = [10, 20, 30]
  sage: L
  [10, 20, 30]
  sage: []                    # [] est la liste vide
  []

Sage example in ./programmation.tex, line 1361::

  sage: L[1], len(L), len([])
  (20, 3, 0)

Sage example in ./programmation.tex, line 1368::

  sage: L[2] = 33
  sage: L
  [10, 20, 33]

Sage example in ./programmation.tex, line 1375::

  sage: L = [11, 22, 33]
  sage: L[-1], L[-2], L[-3]
  (33, 22, 11)

Sage example in ./programmation.tex, line 1388::

  sage: L = [0, 11, 22, 33, 44, 55]
  sage: L[2:4]
  [22, 33]
  sage: L[-4:4]
  [22, 33]
  sage: L[2:-2]
  [22, 33]
  sage: L[:4]
  [0, 11, 22, 33]
  sage: L[4:]
  [44, 55]

Sage example in ./programmation.tex, line 1404::

  sage: L = [0, 11, 22, 33, 44, 55, 66, 77]
  sage: L[2:6] = [12, 13, 14]        # remplace [22, 33, 44, 55]

Sage example in ./programmation.tex, line 1429::

  sage: L = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19]
  sage: L[3:len(L)-5] == L[3-len(L):-5]
  True
  sage: [5 in L, 6 in L]
  [True, False]

Sage example in ./programmation.tex, line 1448::

  sage: L = [1, 2, 3] ; L + [10, 20, 30]
  [1, 2, 3, 10, 20, 30]
  sage: 4 * [1, 2, 3]
  [1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3]

Sage example in ./programmation.tex, line 1467::

  sage: L = 5*[10, 20, 30] ; L[:3]+L[3:] == L
  True

Sage example in ./programmation.tex, line 1476::

  sage: [1..3, 7, 10..13]
  [1, 2, 3, 7, 10, 11, 12, 13]

Sage example in ./programmation.tex, line 1493::

  sage: map (cos, [0, pi/6, pi/4, pi/3, pi/2])
  [1, 1/2*sqrt(3), 1/2*sqrt(2), 1/2, 0]

Sage example in ./programmation.tex, line 1501::

  sage: map (lambda t: cos(t), [0, pi/6, pi/4, pi/3, pi/2])
  [1, 1/2*sqrt(3), 1/2*sqrt(2), 1/2, 0]

Sage example in ./programmation.tex, line 1526::

  sage: map (lambda t: N(cos(t)), [0, pi/6, pi/4, pi/3, pi/2])
  [1.00000000000000, 0.866025403784439, 0.707106781186548,
  0.500000000000000, 0.000000000000000]

Sage example in ./programmation.tex, line 1538::

  sage: map (N, map (cos, [0, pi/6, pi/4, pi/3, pi/2]))
  [1.00000000000000, 0.866025403784439, 0.707106781186548,
  0.500000000000000, 0.000000000000000]

Sage example in ./programmation.tex, line 1543::

  sage: map (compose(N, cos), [0, pi/6, pi/4, pi/3, pi/2])
  [1.00000000000000, 0.866025403784439, 0.707106781186548,
  0.500000000000000, 0.000000000000000]

Sage example in ./programmation.tex, line 1552::

  sage: filter (is_prime, [1..55])
  [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]

Sage example in ./programmation.tex, line 1562::

  sage: p = 37 ; filter (lambda n: n^4 % p == 7, [0..p-1])
  [3, 18, 19, 34]

Sage example in ./programmation.tex, line 1571::

  sage: map(lambda n:2*n+1, [0..15])
  [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31]
  sage: [2*n+1 for n in [0..15]]
  [1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31]

Sage example in ./programmation.tex, line 1580::

  sage: filter (is_prime, [1..55])
  [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]
  sage: [p for p in [1..55] if is_prime(p)]
  [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]

Sage example in ./programmation.tex, line 1590::

  sage: filter (is_prime, [4*n+1 for n in [0..20]])
  [5, 13, 17, 29, 37, 41, 53, 61, 73]
  sage: [n^2 for n in [1..20] if is_prime(n)]
  [4, 9, 25, 49, 121, 169, 289, 361]

Sage example in ./programmation.tex, line 1609::

  sage: reduce (lambda x, y: 10*x+y, [1, 2, 3, 4])
  1234

Sage example in ./programmation.tex, line 1614::

  sage: reduce (lambda x, y: 10*x+y, [9, 8, 7, 6], 1)
  19876

Sage example in ./programmation.tex, line 1621::

  sage: L = [2*n+1 for n in [0..9]]
  sage: reduce (lambda x, y: x*y, L, 1)
  654729075

Sage example in ./programmation.tex, line 1672::

  sage: prod ([2*n+1 for n in [0..9]], 1) # une liste avec for
  654729075
  sage: prod ( 2*n+1 for n in [0..9])     # sans liste
  654729075
  sage: prod (n for n in [0..19] if n%2 == 1)
  654729075

Sage example in ./programmation.tex, line 1686::

  sage: def fct (x): return 4/x == 2
  sage: all (fct(x) for x in [2, 1, 0])
  False
  sage: any (fct(x) for x in [2, 1, 0])
  True

Sage example in ./programmation.tex, line 1704::

  sage: [[x, y] for x in [1..2] for y in [6..8]]
  [[1, 6], [1, 7], [1, 8], [2, 6], [2, 7], [2, 8]]

Sage example in ./programmation.tex, line 1709::

  sage: [[[x, y] for x in [1..2]] for y in [6..8]]
  [[[1, 6], [2, 6]], [[1, 7], [2, 7]], [[1, 8], [2, 8]]]

Sage example in ./programmation.tex, line 1716::

  sage: map (lambda x, y: [x, y], [1..3], [6..8])
  [[1, 6], [2, 7], [3, 8]]

Sage example in ./programmation.tex, line 1723::

  sage: L = [[1, 2, [3]], [4, [5, 6]], [7, [8, [9]]]]
  sage: flatten (L, max_level = 1)
  [1, 2, [3], 4, [5, 6], 7, [8, [9]]]
  sage: flatten (L, max_level = 2)
  [1, 2, 3, 4, 5, 6, 7, 8, [9]]
  sage: flatten (L)             # équivaut à flatten (L, max_level = 3)
  [1, 2, 3, 4, 5, 6, 7, 8, 9]

Sage example in ./programmation.tex, line 1740::

  sage: x = var('x')
  sage: factor(diff(x*exp(x), [x, x]))
  (x + 2)*e^x
  sage: map(lambda n: factor(diff(x*exp(x), n*[x])), [0..6])
  [x*e^x, (x + 1)*e^x, (x + 2)*e^x, (x + 3)*e^x, (x + 4)*e^x,
  (x + 5)*e^x, (x + 6)*e^x]
  sage: [factor (diff (x*exp(x), n*[x])) for n in [0..6]]
  [x*e^x, (x + 1)*e^x, (x + 2)*e^x, (x + 3)*e^x, (x + 4)*e^x,
  (x + 5)*e^x, (x + 6)*e^x]

Sage example in ./programmation.tex, line 1774::

  sage: L = [1, 8, 5, 2, 9] ; L.reverse() ; L
  [9, 2, 5, 8, 1]
  sage: L.sort() ; L
  [1, 2, 5, 8, 9]
  sage: L.sort(reverse = True) ; L
  [9, 8, 5, 2, 1]

Sage example in ./programmation.tex, line 1818::

  sage: def alpha (P, Q):     # len(P) = len(Q) par hypothèse
  ....:   i = 0
  ....:   while True:
  ....:     if i == len(P): return int(0)
  ....:     elif P[i] < Q[i]: return int(-1)
  ....:     elif P[i] > Q[i]: return int(1)
  ....:     else: i = i+1
  sage: alpha ([2, 3, 4, 6, 5], [2, 3, 4, 5, 6])
  1

Sage example in ./programmation.tex, line 1835::

  sage: L = [[2, 2, 5], [2, 3, 4], [3, 2, 4], [3, 3, 3],\
  ....: [1, 1, 2], [1, 2, 7]]
  sage: L.sort (cmp = alpha) ; L
  [[1, 1, 2], [1, 2, 7], [2, 2, 5], [2, 3, 4], [3, 2, 4], [3, 3, 3]]

Sage example in ./programmation.tex, line 1856::

  sage: def homogLex (P, Q):
  ....:  sp = sum (P) ; sq = sum (Q)
  ....:  if sp < sq: return int(-1)
  ....:  elif sp > sq: return int(1)
  ....:  else: return alpha (P, Q)

Sage example in ./programmation.tex, line 1863::

  sage: homogLex ([2, 3, 4, 6, 4], [2, 3, 4, 5, 6])
  -1

Sage example in ./programmation.tex, line 1914::

  sage: def fct1(L):
  ....:   return [filter (lambda n: n % 2 == 0, L),
  ....:           filter (lambda n: n % 2 == 1, L)]

Sage example in ./programmation.tex, line 1919::

  sage: fct1([1..10])
  [[2, 4, 6, 8, 10], [1, 3, 5, 7, 9]]

Sage example in ./programmation.tex, line 1926::

  sage: def fct2 (L):
  ....:   res0 = [] ; res1 = []
  ....:   for k in L:
  ....:     if k%2 == 0: res0.append(k) # ou res0[len(res0):] = [k]
  ....:     else: res1.append(k)        # ou res1[len(res1):] = [k]
  ....:   return [res0, res1]

Sage example in ./programmation.tex, line 1936::

  sage: def fct3a (L, res0, res1):
  ....:   if L == []: return [res0, res1]
  ....:   elif L[0]%2 == 0: return fct3a(L[1:], res0+[L[0]], res1)
  ....:   else: return fct3a (L[1:], res0, res1+[L[0]])

Sage example in ./programmation.tex, line 1942::

  sage: def fct3 (L): return fct3a (L, [], [])

Sage example in ./programmation.tex, line 1955::

  sage: def sousSuites (L):
  ....:   if L == []: return []
  ....:   res = [] ; debut = 0 ; k = 1
  ....:   while k < len(L):    # 2 termes consécutifs sont définis
  ....:     if L[k-1] > L[k]:
  ....:       res.append (L[debut:k]) ; debut = k
  ....:     k = k+1
  ....:   res.append (L[debut:k])
  ....:   return res

Sage example in ./programmation.tex, line 1966::

  sage: sousSuites([1, 4, 1, 5])
  [[1, 4], [1, 5]]
  sage: sousSuites([4, 1, 5, 1])
  [[4], [1, 5], [1]]

Sage example in ./programmation.tex, line 1991::

  sage: S = 'Ceci est une chaîne de caractères.'

Sage example in ./programmation.tex, line 2001::

  sage: S = 'Ceci est une chaîne de caractères.'; S
  'Ceci est une cha\xc3\xaene de caract\xc3\xa8res.'
  sage: print S
  Ceci est une chaîne de caractères.

Sage example in ./programmation.tex, line 2026::

  sage: S='un deux trois quatre cinq six sept'; L=S.split(); L
  ['un', 'deux', 'trois', 'quatre', 'cinq', 'six', 'sept']

Sage example in ./programmation.tex, line 2052::

  sage: L1 = [11, 22, 33] ; L2 = L1
  sage: L1[1] = 222 ; L2.sort() ; L1, L2
  ([11, 33, 222], [11, 33, 222])
  sage: L1[2:3] = []; L2[0:0] = [6, 7, 8]
  sage: L1, L2
  ([6, 7, 8, 11, 33], [6, 7, 8, 11, 33])

Sage example in ./programmation.tex, line 2088::

  sage: L1 = [11, 22, 33] ; L2 = L1 ; L3 = L1[:]
  sage: [L1 is L2, L2 is L1, L1 is L3, L1 == L3]
  [True, True, False, True]

Sage example in ./programmation.tex, line 2096::

  sage: La = [1, 2, 3] ; L1 = [1, La] ; L2 = copy(L1)
  sage: L1[1][0] = 5         # [1, [5, 2, 3]] pour L1 et L2
  sage: [L1 == L2, L1 is L2, L1[1] is L2[1]]
  [True, False, True]

Sage example in ./programmation.tex, line 2141::

  sage: def lexInverse (P, Q):
  ....:   P1 = copy(P) ; P1.reverse()
  ....:   Q1 = copy(Q) ; Q1.reverse()
  ....:   return - alpha (P1, Q1)

Sage example in ./programmation.tex, line 2202::

  sage: S0 = (); S1 = (1, ); S2 = (1, 2)
  sage: [1 in S1, 1 == (1)]
  [True, True]

Sage example in ./programmation.tex, line 2214::

  sage: S1 = (1, 4, 9, 16, 25); [k for k in S1]
  [1, 4, 9, 16, 25]

Sage example in ./programmation.tex, line 2221::

  sage: L1 = [0..4]; L2 = [5..9]
  sage: zip(L1, L2)
  [(0, 5), (1, 6), (2, 7), (3, 8), (4, 9)]
  sage: map(lambda x, y:(x, y), L1, L2)
  [(0, 5), (1, 6), (2, 7), (3, 8), (4, 9)]

Sage example in ./programmation.tex, line 2238::

  sage: E=Set([1, 2, 4, 8, 2, 2, 2]); F=Set([7, 5, 3, 1]); E, F
  ({8, 1, 2, 4}, {1, 3, 5, 7})

Sage example in ./programmation.tex, line 2252::

  sage: E = Set([1, 2, 4, 8, 2, 2, 2]); F = Set([7, 5, 3, 1])
  sage: 5 in E, 5 in F, E + F == F | E
  (False, True, True)
  sage: E & F, E - F, E ^^ F
  ({1}, {8, 2, 4}, {2, 3, 4, 5, 7, 8})

Sage example in ./programmation.tex, line 2267::

  sage: E = Set([1, 2, 4, 8, 2, 2, 2])
  sage: [E[k] for k in [0..len(E)-1]], [t for t in E]
  ([8, 1, 2, 4], [8, 1, 2, 4])

Sage example in ./programmation.tex, line 2275::

  sage: def inclus (E, F): return E+F == F

Sage example in ./programmation.tex, line 2283::

  sage: Set([Set([]), Set([1]), Set([2]), Set([1, 2])])
  {{1, 2}, {}, {2}, {1}}
  sage: Set([    (),     (1, ),    (2, ),    (1, 2) ])
  {(1, 2), (2,), (), (1,)}

Sage example in ./programmation.tex, line 2292::

  sage: def Parties (EE):
  ....:   if EE == Set([]): return Set([EE])
  ....:   else:
  ....:     return avecOuSansElt (EE[0], Parties(Set(EE[1:])))

Sage example in ./programmation.tex, line 2298::

  sage: def avecOuSansElt (a, E):
  ....:   return Set (map (lambda F: Set([a])+F, E)) + E

Sage example in ./programmation.tex, line 2302::

  sage: Parties(Set([1, 2, 3]))
  {{3}, {1, 2}, {}, {2, 3}, {1}, {1, 3}, {1, 2, 3}, {2}}

Sage example in ./programmation.tex, line 2324::

  sage: D={}; D['un']=1; D['deux']=2; D['trois']=3; D['dix']=10
  sage: D['deux'] + D['trois']
  5

Sage example in ./programmation.tex, line 2354::

  sage: D = {'a0':'b0', 'a1':'b1', 'a2':'b2', 'a3':'b0',\
  ....: 'a4':'b3', 'a5':'b3'}
  sage: E  = Set(D.keys()) ; Imf = Set(D.values())
  sage: Imf == Set(map (lambda t:D[t], E))     # est équivalent
  True

Sage example in ./programmation.tex, line 2381::

  sage: def injective(D):
  ....:   return len(D) == len (Set(D.values()))

"""
