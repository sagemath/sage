# sage.doctest: needs sage.combinat sage.graphs sage.groups
r"""
This file contains doctests for the Chapter "k-Schur function primer"
for the book "k-Schur functions and affine Schubert calculus"
by Thomas Lam, Luc Lapointe, Jennifer Morse, Anne Schilling, Mark Shimozono,
and Mike Zabrocki, :arxiv:`1301.3569`.
The code was written by Anne Schilling and Mike Zabrocki, 2012 and 2013.

IF IT BECOMES NECESSARY TO CHANGE ANY TESTS IN THIS FILE, THERE
NEEDS TO BE A ONE-YEAR DEPRECATION PERIOD. ALSO, PLEASE IN THIS CASE
CONTACT Anne Schilling (anne@math.ucdavis.edu) AND Mike Zabrocki
(zabrocki@mathstat.yorku.ca) REGARDING THE CHANGES!
"""
"""
Sage example in ./kschurnotes/notes-mike-anne.tex, line 198::

    sage: P = Partitions(4); P
    Partitions of the integer 4
    sage: P.list()
    [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 205::

    sage: [p for p in P]
    [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 210::

    sage: la=Partition([2,2]); mu=Partition([3,1])
    sage: mu.dominates(la)
    True

Sage example in ./kschurnotes/notes-mike-anne.tex, line 216::

    sage: ord = lambda x,y: y.dominates(x)
    sage: P = Poset([Partitions(6), ord], facade=True)
    sage: H = P.hasse_diagram()

Sage example in ./kschurnotes/notes-mike-anne.tex, line 228::

    sage: la=Partition([4,3,3,3,2,2,1])
    sage: la.conjugate()
    [7, 6, 4, 1]
    sage: la.k_split(4)
    [[4], [3, 3], [3, 2], [2, 1]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 236::

    sage: p = SkewPartition([[2,1],[1]])
    sage: p.is_connected()
    False

Sage example in ./kschurnotes/notes-mike-anne.tex, line 334::

    sage: la = Partition([4,3,3,3,2,2,1])
    sage: kappa = la.k_skew(4); kappa
    [12, 8, 5, 5, 2, 2, 1] / [8, 5, 2, 2]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 340::

    sage: kappa.row_lengths()
    [4, 3, 3, 3, 2, 2, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 345::

    sage: tau = Core([12,8,5,5,2,2,1],5)
    sage: mu = tau.to_bounded_partition(); mu
    [4, 3, 3, 3, 2, 2, 1]
    sage: mu.to_core(4)
    [12, 8, 5, 5, 2, 2, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 353::

    sage: Cores(3,6).list()
    [[6, 4, 2], [5, 3, 1, 1], [4, 2, 2, 1, 1], [3, 3, 2, 2, 1, 1]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 398::

    sage: W = WeylGroup(['A',4,1])      # long time (5.47 s, 2013)
    sage: S = W.simple_reflections()    # long time
    sage: [s.reduced_word() for s in S] # long time
    [[0], [1], [2], [3], [4]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 406::

    sage: # long time
    sage: w = W.an_element(); w
    [ 2  0  0  1 -2]
    [ 2  0  0  0 -1]
    [ 1  1  0  0 -1]
    [ 1  0  1  0 -1]
    [ 1  0  0  1 -1]
    sage: w.reduced_word()
    [0, 1, 2, 3, 4]
    sage: w = W.from_reduced_word([2,1,0])
    sage: w.is_affine_grassmannian()
    True

Sage example in ./kschurnotes/notes-mike-anne.tex, line 464::

    sage: c = Core([7,3,1],5)
    sage: c.affine_symmetric_group_simple_action(2)
    [8, 4, 1, 1]
    sage: c.affine_symmetric_group_simple_action(0)
    [7, 3, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 474::

    sage: k=4; length=3
    sage: W = WeylGroup(['A',k,1])
    sage: G = W.affine_grassmannian_elements_of_given_length(length)
    sage: [w.reduced_word() for w in G]
    [[2, 1, 0], [4, 1, 0], [3, 4, 0]]

    sage: C = Cores(k+1,length)
    sage: [c.to_grassmannian().reduced_word() for c in C]
    [[2, 1, 0], [4, 1, 0], [3, 4, 0]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 543::

    sage: la = Partition([4,3,3,3,2,2,1])
    sage: c = la.to_core(4); c
    [12, 8, 5, 5, 2, 2, 1]
    sage: W = WeylGroup(['A',4,1])
    sage: w = W.from_reduced_word([4,1,0,2,1,4,3,2,0,4,3,1,0,4,3,2,1,0])
    sage: c.to_grassmannian() == w
    True

Sage example in ./kschurnotes/notes-mike-anne.tex, line 643::

sage: A = AffinePermutationGroup(['A',2,1])
sage: w = A([-2,0,8])
sage: w.reduced_word()
[1, 0, 2, 1, 0]
sage: w.to_core()
[5, 3, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 767::

    sage: la = Partition([4,3,3,3,2,2,1])
    sage: la.k_conjugate(4)
    [3, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1065::

    sage: c = Core([3,1,1],3)
    sage: c.weak_covers()
    [[4, 2, 1, 1]]
    sage: c.strong_covers()
    [[5, 3, 1], [4, 2, 1, 1], [3, 2, 2, 1, 1]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1073::

    sage: kappa = Core([4,1],4)
    sage: tau = Core([2,1],4)
    sage: tau.weak_le(kappa)
    False
    sage: tau.strong_le(kappa)
    True

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1082::

    sage: C = sum(([c for c in Cores(4,m)] for m in range(7)),[])
    sage: ord = lambda x,y: x.weak_le(y)
    sage: P = Poset([C, ord], cover_relations = False)  # long time (3.99 s, 2013)
    sage: H = P.hasse_diagram()                         # long time
    sage: view(H)  # not tested

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1253::

    sage: Sym = SymmetricFunctions(QQ)
    sage: h = Sym.homogeneous()
    sage: m = Sym.monomial()

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1259::

    sage: f = h[3,1]+h[2,2]
    sage: m(f)
    10*m[1, 1, 1, 1] + 7*m[2, 1, 1] + 5*m[2, 2] + 4*m[3, 1] + 2*m[4]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1266::

    sage: f.scalar(h[2,1,1])
    7
    sage: m(f).coefficient([2,1,1])
    7

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1274::

    sage: p = Sym.power()
    sage: e = Sym.elementary()
    sage: sum( (-1)**(i-1)*e[4-i]*p[i] for i in range(1,4) ) - p[4]
    4*e[4]
    sage: sum( (-1)**(i-1)*p[i]*e[4-i] for i in range(1,4) ) - p[4]
    1/6*p[1, 1, 1, 1] - p[2, 1, 1] + 1/2*p[2, 2] + 4/3*p[3, 1] - p[4]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1327::

    sage: Sym = SymmetricFunctions(QQ)
    sage: s = Sym.schur()
    sage: m = Sym.monomial()
    sage: h = Sym.homogeneous()
    sage: m(s[1,1,1])
    m[1, 1, 1]
    sage: h(s[1,1,1])
    h[1, 1, 1] - 2*h[2, 1] + h[3]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1338::

    sage: p = Sym.power()
    sage: s = Sym.schur()
    sage: p(s[1,1,1])
    1/6*p[1, 1, 1] - 1/2*p[2, 1] + 1/3*p[3]
    sage: p(s[2,1])
    1/3*p[1, 1, 1] - 1/3*p[3]
    sage: p(s[3])
    1/6*p[1, 1, 1] + 1/2*p[2, 1] + 1/3*p[3]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1349::

    sage: s[2,1].scalar(s[1,1,1])
    0
    sage: s[2,1].scalar(s[2,1])
    1

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1517::

  sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
  sage: Qp = Sym.hall_littlewood().Qp()
  sage: Qp.base_ring()
  Fraction Field of Univariate Polynomial Ring in t over Rational Field
  sage: s = Sym.schur()
  sage: s(Qp[1,1,1])
  s[1, 1, 1] + (t^2+t)*s[2, 1] + t^3*s[3]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1530::

  sage: t = Qp.t
  sage: s[2,1].scalar(s[3].theta_qt(t,0))
  t^2 - t

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1536::

  sage: s(Qp([1,1])).hl_creation_operator([3])
  s[3, 1, 1] + t*s[3, 2] + (t^2+t)*s[4, 1] + t^3*s[5]
  sage: s(Qp([3,1,1]))
  s[3, 1, 1] + t*s[3, 2] + (t^2+t)*s[4, 1] + t^3*s[5]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1568::

  sage: Sym = SymmetricFunctions(FractionField(QQ['q,t']))
  sage: Mac = Sym.macdonald()
  sage: H = Mac.H()
  sage: s = Sym.schur()
  sage: for la in Partitions(3):
  ....:     print("H {} = {}".format(la, s(H(la))))
  H [3] = q^3*s[1, 1, 1] + (q^2+q)*s[2, 1] + s[3]
  H [2, 1] = q*s[1, 1, 1] + (q*t+1)*s[2, 1] + t*s[3]
  H [1, 1, 1] = s[1, 1, 1] + (t^2+t)*s[2, 1] + t^3*s[3]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1581::

  sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
  sage: Mac = Sym.macdonald(q=0)
  sage: H = Mac.H()
  sage: s = Sym.schur()
  sage: for la in Partitions(3):
  ....:     print("H {} = {}".format(la, s(H(la))))
  H [3] = s[3]
  H [2, 1] = s[2, 1] + t*s[3]
  H [1, 1, 1] = s[1, 1, 1] + (t^2+t)*s[2, 1] + t^3*s[3]
  sage: Qp = Sym.hall_littlewood().Qp()
  sage: s(Qp[1, 1, 1])
  s[1, 1, 1] + (t^2+t)*s[2, 1] + t^3*s[3]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1596::

  sage: Sym = SymmetricFunctions(FractionField(QQ['q']))
  sage: Mac = Sym.macdonald(t=0)
  sage: H = Mac.H()
  sage: s = Sym.schur()
  sage: for la in Partitions(3):
  ....:     print("H {} = {}".format(la, s(H(la))))
  H [3] = q^3*s[1, 1, 1] + (q^2+q)*s[2, 1] + s[3]
  H [2, 1] = q*s[1, 1, 1] + s[2, 1]
  H [1, 1, 1] = s[1, 1, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1717::

    sage: t = Tableau([[1,1,1,2,3,7],[2,2,3,5],[3,4],[4,5],[6]])
    sage: t.charge()
    9

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1881::

    sage: la = Partition([2,2])
    sage: la.k_conjugate(2).conjugate()
    [4]
    sage: la = Partition([2,1,1])
    sage: la.k_conjugate(2).conjugate()
    [3, 1]
    sage: la = Partition([1,1,1,1])
    sage: la.k_conjugate(2).conjugate()
    [2, 2]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1893::

    sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
    sage: ks = Sym.kschur(2)
    sage: ks[2,2].omega_t_inverse()
    1/t^2*ks2[1, 1, 1, 1]
    sage: ks[2,1,1].omega_t_inverse()
    1/t*ks2[2, 1, 1]
    sage: ks[1,1,1,1].omega_t_inverse()
    1/t^2*ks2[2, 2]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 1904::

    sage: Sym = SymmetricFunctions(FractionField(QQ['q,t']))
    sage: H = Sym.macdonald().H()
    sage: ks = Sym.kschur(2)
    sage: ks(H[2,2])
    q^2*ks2[1, 1, 1, 1] + (q*t+q)*ks2[2, 1, 1] + ks2[2, 2]
    sage: ks(H[2,1,1])
    q*ks2[1, 1, 1, 1] + (q*t^2+1)*ks2[2, 1, 1] + t*ks2[2, 2]
    sage: ks(H[1,1,1,1])
    ks2[1, 1, 1, 1] + (t^3+t^2)*ks2[2, 1, 1] + t^4*ks2[2, 2]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 2174::

    sage: SemistandardTableaux([5,2],[4,2,1]).list()
    [[[1, 1, 1, 1, 2], [2, 3]], [[1, 1, 1, 1, 3], [2, 2]]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 2179::

    sage: P = Partitions(4)
    sage: P.list()
    [[4], [3, 1], [2, 2], [2, 1, 1], [1, 1, 1, 1]]
    sage: n = P.cardinality(); n
    5
    sage: K = matrix(QQ,n,n,
    ....:           [[SemistandardTableaux(la,mu).cardinality()
    ....:            for mu in P] for la in P])
    sage: K
    [1 1 1 1 1]
    [0 1 1 2 3]
    [0 0 1 1 2]
    [0 0 0 1 3]
    [0 0 0 0 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 2376::

    sage: T = WeakTableaux(6, [5,3], [4,3,1])
    sage: T.list()
    [[[1, 1, 1, 1, 3], [2, 2, 2]], [[1, 1, 1, 1, 2], [2, 2, 3]]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 2382::

    sage: k = 3
    sage: c = Core([5,2,1], k+1)
    sage: la = c.to_bounded_partition(); la
    [3, 2, 1]
    sage: for mu in Partitions(la.size(), max_part = 3):
    ....:     T = WeakTableaux(k, c, mu)
    ....:     print("weight {}".format(mu))
    ....:     print(T.list())
    weight [3, 3]
    []
    weight [3, 2, 1]
    [[[1, 1, 1, 2, 2], [2, 2], [3]]]
    weight [3, 1, 1, 1]
    [[[1, 1, 1, 2, 4], [2, 4], [3]], [[1, 1, 1, 2, 3], [2, 3], [4]]]
    weight [2, 2, 2]
    [[[1, 1, 2, 2, 3], [2, 3], [3]]]
    weight [2, 2, 1, 1]
    [[[1, 1, 2, 2, 4], [2, 4], [3]], [[1, 1, 2, 2, 3], [2, 3], [4]]]
    weight [2, 1, 1, 1, 1]
    [[[1, 1, 3, 4, 5], [2, 5], [3]], [[1, 1, 2, 3, 5], [3, 5], [4]],
     [[1, 1, 2, 3, 4], [3, 4], [5]]]
    weight [1, 1, 1, 1, 1, 1]
    [[[1, 3, 4, 5, 6], [2, 6], [4]], [[1, 2, 4, 5, 6], [3, 6], [4]],
     [[1, 2, 3, 4, 6], [4, 6], [5]], [[1, 2, 3, 4, 5], [4, 5], [6]]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 2487::

    sage: Sym = SymmetricFunctions(QQ)
    sage: ks = Sym.kschur(3, t=1)
    sage: h = Sym.homogeneous()
    sage: for mu in Partitions(7, max_part=3):
    ....:     print(h(ks(mu)))
    h[3, 3, 1]
    h[3, 2, 2] - h[3, 3, 1]
    h[3, 2, 1, 1] - h[3, 2, 2]
    h[3, 1, 1, 1, 1] - 2*h[3, 2, 1, 1] + h[3, 3, 1]
    h[2, 2, 2, 1] - h[3, 2, 1, 1] - h[3, 2, 2] + h[3, 3, 1]
    h[2, 2, 1, 1, 1] - 2*h[2, 2, 2, 1] - h[3, 1, 1, 1, 1]
        + 2*h[3, 2, 1, 1] + h[3, 2, 2] - h[3, 3, 1]
    h[2, 1, 1, 1, 1, 1] - 3*h[2, 2, 1, 1, 1] + 2*h[2, 2, 2, 1]
        + h[3, 2, 1, 1] - h[3, 2, 2]
    h[1, 1, 1, 1, 1, 1, 1] - 4*h[2, 1, 1, 1, 1, 1] + 4*h[2, 2, 1, 1, 1]
        + 2*h[3, 1, 1, 1, 1] - 4*h[3, 2, 1, 1] + h[3, 3, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 2608::

    sage: ks6 = Sym.kschur(6,t=1)
    sage: ks6(h[4,3,1])
    ks6[4, 3, 1] + ks6[4, 4] + ks6[5, 2, 1] + 2*ks6[5, 3]
      + ks6[6, 1, 1] + ks6[6, 2]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 2617::

    sage: Sym = SymmetricFunctions(QQ)
    sage: ks = Sym.kschur(3,t=1)
    sage: ks.realization_of()
    3-bounded Symmetric Functions over Rational Field with t=1
    sage: s = Sym.schur()
    sage: s.realization_of()
    Symmetric Functions over Rational Field

Sage example in ./kschurnotes/notes-mike-anne.tex, line 2658::

    sage: k = 6
    sage: weight = Partition([4,3,1])
    sage: for la in Partitions(weight.size(), max_part = k):
    ....:     if la.dominates(weight):
    ....:         print(la)
    ....:         T = WeakTableaux(k, la, weight, representation = 'bounded')
    ....:         print(T.list())
    [6, 2]
    [[[1, 1, 1, 1, 2, 2], [2, 3]]]
    [6, 1, 1]
    [[[1, 1, 1, 1, 2, 2], [2], [3]]]
    [5, 3]
    [[[1, 1, 1, 1, 3], [2, 2, 2]], [[1, 1, 1, 1, 2], [2, 2, 3]]]
    [5, 2, 1]
    [[[1, 1, 1, 1, 2], [2, 2], [3]]]
    [4, 4]
    [[[1, 1, 1, 1], [2, 2, 2, 3]]]
    [4, 3, 1]
    [[[1, 1, 1, 1], [2, 2, 2], [3]]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 2740::

    sage: mu = Partition([3,2,1])
    sage: c = mu.to_core(3)
    sage: w = c.to_grassmannian()
    sage: w.stanley_symmetric_function()
    4*m[1, 1, 1, 1, 1, 1] + 3*m[2, 1, 1, 1, 1] + 2*m[2, 2, 1, 1]
        + m[2, 2, 2] + 2*m[3, 1, 1, 1] + m[3, 2, 1]
    sage: sorted(w.reduced_words())
    [[0, 2, 3, 2, 1, 0],
     [0, 3, 2, 1, 3, 0],
     [0, 3, 2, 3, 1, 0],
     [2, 0, 3, 2, 1, 0]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 2752::

    sage: Sym = SymmetricFunctions(QQ)
    sage: Q3 = Sym.kBoundedQuotient(3,t=1)
    sage: F3 = Q3.affineSchur()
    sage: m = Q3.kmonomial()
    sage: m(F3([3,2,1]))
    4*m3[1, 1, 1, 1, 1, 1] + 3*m3[2, 1, 1, 1, 1] + 2*m3[2, 2, 1, 1]
        + m3[2, 2, 2] + 2*m3[3, 1, 1, 1] + m3[3, 2, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 2799::

    sage: Sym = SymmetricFunctions(QQ)
    sage: Q3 = Sym.kBoundedQuotient(3,t=1)
    sage: F3 = Q3.affineSchur()
    sage: h = Sym.homogeneous()
    sage: f = F3[3,2,1]*h[1]; f
    F3[3, 1, 1, 1, 1] + 3*F3[3, 2, 1, 1] + F3[3, 2, 2] + 2*F3[3, 3, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 2810::

    sage: c = Partition([3,2,1]).to_core(3)
    sage: for p in sorted(f.support()):   # Sorted for consistant doctest ordering
    ....:   print("{} {}".format(p, SkewPartition([p.to_core(3).to_partition(),c.to_partition()])))
    [3, 1, 1, 1, 1] [[5, 2, 1, 1, 1], [5, 2, 1]]
    [3, 2, 1, 1] [[6, 3, 1, 1], [5, 2, 1]]
    [3, 2, 2] [[5, 2, 2], [5, 2, 1]]
    [3, 3, 1] [[7, 4, 1], [5, 2, 1]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 2976::

sage: T = StrongTableau([[-1,-1,-2,-3],[-2,3,-3,4],[2,3],[-3,-4]], 3)
sage: T.to_transposition_sequence()
[[-2, -1], [3, 4], [0, 2], [-3, -2], [2, 3], [-1, 0], [1, 2], [0, 1]]
sage: T.intermediate_shapes()
[[], [2], [3, 1, 1], [4, 3, 2, 1], [4, 4, 2, 2]]
sage: [T.content_of_marked_head(v+1) for v in range(8)]
[0, 1, -1, 2, -3, 1, 3, -2]
sage: T.left_action([0,1])
[[-1, -1, -2, -3, 5], [-2, 3, -3, 4], [2, 3, -5], [-3, -4], [5]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 2999::

sage: ST = StrongTableaux(3, [6,3,1,1], [4,2,1]); ST
Set of strong 3-tableaux of shape [6, 3, 1, 1] and of weight (4, 2, 1)
sage: ST.list()
[[[-1, -1, -1, -1, 2, 2], [1, -2, -2], [-3], [3]],
 [[-1, -1, -1, -1, 2, -2], [1, -2, 2], [-3], [3]],
 [[-1, -1, -1, -1, -2, -2], [1, 2, 2], [-3], [3]],
 [[-1, -1, -1, -1, 2, 3], [1, -2, 3], [-2], [-3]],
 [[-1, -1, -1, -1, 2, 3], [1, -2, -3], [-2], [3]],
 [[-1, -1, -1, -1, 2, -3], [1, -2, 3], [-2], [3]],
 [[-1, -1, -1, -1, -2, 3], [1, 2, 3], [-2], [-3]],
 [[-1, -1, -1, -1, -2, 3], [1, 2, -3], [-2], [3]],
 [[-1, -1, -1, -1, -2, -3], [1, 2, 3], [-2], [3]]]
sage: ks = SymmetricFunctions(QQ).kschur(3,1)
sage: m = SymmetricFunctions(QQ).m()
sage: m(ks[3,2,1,1]).coefficient([4,2,1])
9

Sage example in ./kschurnotes/notes-mike-anne.tex, line 3243::

  sage: W = WeylGroup(['A',3,1])
  sage: [w.reduced_word() for w in W.pieri_factors()]
  [[], [0], [1], [2], [3], [1, 0], [2, 0], [0, 3], [2, 1], [3, 1], [3, 2],
   [2, 1, 0], [1, 0, 3], [0, 3, 2], [3, 2, 1]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 3251::

  sage: A = NilCoxeterAlgebra(WeylGroup(['A',3,1]), prefix = 'A')
  sage: A.homogeneous_noncommutative_variables([2])
  A[1,0] + A[2,0] + A[0,3] + A[3,2] + A[3,1] + A[2,1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 3258::

  sage: A.k_schur_noncommutative_variables([2,2])
  A[0,3,1,0] + A[3,1,2,0] + A[1,2,0,1] + A[3,2,0,3] + A[2,0,3,1] + A[2,3,1,2]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 3265::

  sage: Sym = SymmetricFunctions(ZZ)
  sage: ks = Sym.kschur(5,t=1)
  sage: ks[2,1]*ks[2,1]
  ks5[2, 2, 1, 1] + ks5[2, 2, 2] + ks5[3, 1, 1, 1] + 2*ks5[3, 2, 1]
  + ks5[3, 3] + ks5[4, 2]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 3536::

    sage: la = Partition([3,2,1,1])
    sage: la.k_atom(4)
    [[[1, 1, 1, 4], [2, 2], [3]], [[1, 1, 1], [2, 2], [3], [4]]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 3639::

    sage: s = SymmetricFunctions(QQ['t']).schur()
    sage: G1 = s[1]
    sage: G211 = G1.hl_creation_operator([2,1]); G211
    s[2, 1, 1] + t*s[2, 2] + t*s[3, 1]
    sage: G3211 = G211.hl_creation_operator([3]); G3211
    s[3, 2, 1, 1] + t*s[3, 2, 2] + t*s[3, 3, 1] + t*s[4, 1, 1, 1]
     + (2*t^2+t)*s[4, 2, 1] + t^2*s[4, 3] + (t^3+t^2)*s[5, 1, 1]
     + 2*t^3*s[5, 2] + t^4*s[6, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 3954::

  sage: T = WeakTableau([[1,1,2,3,4,4,5,5,6],[2,3,5,5,6],[3,4,7],
  ....:  [5,6],[6],[7]],4)
  sage: T.k_charge()
  12

Sage example in ./kschurnotes/notes-mike-anne.tex, line 3962::

  sage: Sym = SymmetricFunctions(QQ['t'].fraction_field())
  sage: Qp = Sym.hall_littlewood().Qp()
  sage: ks = Sym.kBoundedSubspace(3).kschur()
  sage: t = ks.base_ring().gen()
  sage: ks(Qp[3,2,2,1])
  ks3[3, 2, 2, 1] + t*ks3[3, 3, 1, 1] + t^2*ks3[3, 3, 2]
  sage: sum(t^T.k_charge()*ks(la) for la in Partitions(8, max_part=3)
  ....:  for T in WeakTableaux(3,la,[3,2,2,1],representation = 'bounded'))
  ks3[3, 2, 2, 1] + t*ks3[3, 3, 1, 1] + t^2*ks3[3, 3, 2]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 4055::

    sage: t = var('t')                                                                  # needs sage.symbolic
    sage: for mu in Partitions(5):                                                      # needs sage.symbolic
    ....:     print("{} {}".format(mu, sum(t^T.spin() for T in StrongTableaux(3,[4,1,1],mu))))
    [5] 0
    [4, 1] t
    [3, 2] t
    [3, 1, 1] 2*t + 1
    [2, 2, 1] 2*t + 1
    [2, 1, 1, 1] 3*t + 3
    [1, 1, 1, 1, 1] 4*t + 6
    sage: StrongTableaux( 3, [4,1,1], (1,)*5 ).cardinality()
    10
    sage: StrongTableaux( 3, [4,1,1], (1,)*5 ).list()
    [[[-1, -2, -3, 4], [-4], [-5]],
     [[-1, -2, -3, -4], [4], [-5]],
     [[-1, -2, -3, -5], [-4], [4]],
     [[-1, -2, 4, -4], [-3], [-5]],
     [[-1, -2, 4, -5], [-3], [-4]],
     [[-1, -2, -4, -5], [-3], [4]],
     [[-1, -3, 4, -4], [-2], [-5]],
     [[-1, -3, 4, -5], [-2], [-4]],
     [[-1, -3, -4, -5], [-2], [4]],
     [[-1, 4, -4, -5], [-2], [-3]]]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 4385::

    sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
    sage: ks4 = Sym.kschur(4)
    sage: ks4([3, 1, 1]).hl_creation_operator([1])
    (t-1)*ks4[2, 2, 1, 1] + t^2*ks4[3, 1, 1, 1] + t^3*ks4[3, 2, 1]
     + (t^3-t^2)*ks4[3, 3] + t^4*ks4[4, 1, 1]
    sage: ks4([3, 1, 1]).hl_creation_operator([2])
    t*ks4[3, 2, 1, 1] + t^2*ks4[3, 3, 1] + t^2*ks4[4, 1, 1, 1]
     + t^3*ks4[4, 2, 1]
    sage: ks4([3, 1, 1]).hl_creation_operator([3])
    ks4[3, 3, 1, 1] + t*ks4[4, 2, 1, 1] + t^2*ks4[4, 3, 1]
    sage: ks4([3, 1, 1]).hl_creation_operator([4])
    ks4[4, 3, 1, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 4456::

    sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
    sage: ks3 = Sym.kschur(3)
    sage: ks3([3,2]).omega()
    Traceback (most recent call last):
    ...
    ValueError: t^2*s[1, 1, 1, 1, 1] + t*s[2, 1, 1, 1] + s[2, 2, 1] is not
    in the image

    sage: s = Sym.schur()
    sage: s(ks3[3,2])
    s[3, 2] + t*s[4, 1] + t^2*s[5]
    sage: t = s.base_ring().gen()
    sage: invert = lambda x: s.base_ring()(x.subs(t=1/t))
    sage: ks3(s(ks3([3,2])).omega().map_coefficients(invert))
    1/t^2*ks3[1, 1, 1, 1, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 4478::

    sage: ks3[3,2].omega_t_inverse()
    1/t^2*ks3[1, 1, 1, 1, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 4686::

    sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
    sage: ks3 = Sym.kschur(3)
    sage: ks3[3,1].coproduct()
    ks3[] # ks3[3, 1] + ks3[1] # ks3[2, 1] + (t+1)*ks3[1] # ks3[3]
    + ks3[1, 1] # ks3[2] + ks3[2] # ks3[1, 1] + (t+1)*ks3[2] # ks3[2]
    + ks3[2, 1] # ks3[1] + (t+1)*ks3[3] # ks3[1]  + ks3[3, 1] # ks3[]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 4720::

    sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
    sage: ks2 = Sym.kschur(2)
    sage: ks3 = Sym.kschur(3)
    sage: ks5 = Sym.kschur(5)
    sage: ks5(ks3[2])*ks5(ks2[1])
    ks5[2, 1] + ks5[3]
    sage: ks5(ks3[2])*ks5(ks2[2,1])
    ks5[2, 2, 1] + ks5[3, 1, 1] + (t+1)*ks5[3, 2] + (t+1)*ks5[4, 1]
      + t*ks5[5]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 4779::

    sage: Sym = SymmetricFunctions(FractionField(QQ['t']))
    sage: ks3 = Sym.kschur(3)
    sage: ks4 = Sym.kschur(4)
    sage: ks5 = Sym.kschur(5)
    sage: ks4(ks3[3,2,1,1])
    ks4[3, 2, 1, 1] + t*ks4[3, 3, 1] + t*ks4[4, 1, 1, 1] + t^2*ks4[4, 2, 1]
    sage: ks5(ks3[3,2,1,1])
    ks5[3, 2, 1, 1] + t*ks5[3, 3, 1] + t*ks5[4, 1, 1, 1] + t^2*ks5[4, 2, 1]
     + t^2*ks5[4, 3] + t^3*ks5[5, 1, 1]

    sage: ks5(ks4[3,2,1,1])
    ks5[3, 2, 1, 1]
    sage: ks5(ks4[4,3,3,2,1,1])
    ks5[4, 3, 3, 2, 1, 1] + t*ks5[4, 4, 3, 1, 1, 1]
     + t^2*ks5[5, 3, 3, 1, 1, 1]
    sage: ks5(ks4[4,3,3,2,1,1,1])
    ks5[4, 3, 3, 2, 1, 1, 1] + t*ks5[4, 3, 3, 3, 1, 1]
     + t*ks5[4, 4, 3, 1, 1, 1, 1] + t^2*ks5[4, 4, 3, 2, 1, 1]
     + t^2*ks5[5, 3, 3, 1, 1, 1, 1] + t^3*ks5[5, 3, 3, 2, 1, 1]
     + t^4*ks5[5, 4, 3, 1, 1, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 4858::

    sage: Sym = SymmetricFunctions(FractionField(QQ['q,t']))
    sage: H = Sym.macdonald().H()
    sage: ks = Sym.kschur(3)
    sage: ks(H[3])
    q^3*ks3[1, 1, 1] + (q^2+q)*ks3[2, 1] + ks3[3]
    sage: ks(H[3,2])                          # long time (2.11 s, 2013)
    q^4*ks3[1, 1, 1, 1, 1] + (q^3*t+q^3+q^2)*ks3[2, 1, 1, 1]
     + (q^3*t+q^2*t+q^2+q)*ks3[2, 2, 1]
     + (q^2*t+q*t+q)*ks3[3, 1, 1] + ks3[3, 2]
    sage: ks(H[3,1,1])
    q^3*ks3[1, 1, 1, 1, 1] + (q^3*t^2+q^2+q)*ks3[2, 1, 1, 1]
     + (q^2*t^2+q^2*t+q*t+q)*ks3[2, 2, 1]
     + (q^2*t^2+q*t^2+1)*ks3[3, 1, 1] + t*ks3[3, 2]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 4965::

    sage: Sym = SymmetricFunctions(QQ)
    sage: Q3 = Sym.kBoundedQuotient(3,t=1)
    sage: F = Q3.affineSchur()
    sage: p = Sym.power()
    sage: F[2,1]*p[2]
    -F3[1, 1, 1, 1, 1] - F3[2, 1, 1, 1] + F3[3, 1, 1] + F3[3, 2]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 5028::

    sage: R = QQ[I]; z4 = R.zeta(4)
    sage: Sym = SymmetricFunctions(R)
    sage: ks3z = Sym.kschur(3,t=z4)
    sage: ks3 = Sym.kschur(3,t=1)
    sage: p = Sym.p()
    sage: p(ks3z[2, 2, 2, 2, 2, 2, 2, 2])     # long time (17s on sage.math, 2013)
    1/12*p[4, 4, 4, 4] + 1/4*p[8, 8] - 1/3*p[12, 4]
    sage: p(ks3[2,2])
    1/12*p[1, 1, 1, 1] + 1/4*p[2, 2] - 1/3*p[3, 1]
    sage: p(ks3[2,2]).plethysm(p[4])
    1/12*p[4, 4, 4, 4] + 1/4*p[8, 8] - 1/3*p[12, 4]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 5042::

    sage: ks3z[3, 3, 3, 3]*ks3z[2, 1]   # long time (10s on sage.math, 2013)
    ks3[3, 3, 3, 3, 2, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 5197::

    sage: Sym = SymmetricFunctions(QQ['t'].fraction_field())
    sage: Q3 = Sym.kBoundedQuotient(3)
    sage: dks = Q3.dual_k_Schur()
    sage: dks[2, 1, 1]*dks[3, 2, 1]          # long time (25.7s, 2013)
    (t^7+t^6)*dks3[2, 1, 1, 1, 1, 1, 1, 1, 1]
    + (t^4+t^3+t^2)*dks3[2, 2, 2, 1, 1, 1, 1]
    + (t^3+t^2)*dks3[2, 2, 2, 2, 1, 1]
    + (t^5+2*t^4+2*t^3+t^2)*dks3[2, 2, 2, 2, 2]
    + (t^5+2*t^4+t^3)*dks3[3, 1, 1, 1, 1, 1, 1, 1]
    + (2*t^5+3*t^4+4*t^3+3*t^2+t)*dks3[3, 2, 1, 1, 1, 1, 1]
    + (2*t^2+t+1)*dks3[3, 2, 2, 1, 1, 1]
    + (t^4+3*t^3+4*t^2+3*t+1)*dks3[3, 2, 2, 2, 1]
    + (t^5+t^4+4*t^3+4*t^2+3*t+1)*dks3[3, 3, 1, 1, 1, 1]
    + (2*t^5+3*t^4+5*t^3+6*t^2+4*t+2)*dks3[3, 3, 2, 1, 1]
    + (t^4+t^3+3*t^2+2*t+1)*dks3[3, 3, 2, 2]
    + (t^5+3*t^4+3*t^3+4*t^2+2*t+1)*dks3[3, 3, 3, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 5294::

    sage: ks2 = SymmetricFunctions(QQ['t']).kschur(2)
    sage: HLQp = SymmetricFunctions(QQ['t']).hall_littlewood().Qp()
    sage: ks2( (HLQp(ks2[1,1])*HLQp(ks2[1])).restrict_parts(2) )
    ks2[1, 1, 1] + (-t+1)*ks2[2, 1]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 5303::

    sage: dks = SymmetricFunctions(QQ['t']).kBoundedQuotient(2).dks()
    sage: dks[2,1].coproduct()
    dks2[] # dks2[2, 1] + (-t+1)*dks2[1] # dks2[1, 1] +
     dks2[1] # dks2[2] + (-t+1)*dks2[1, 1] # dks2[1] +
     dks2[2] # dks2[1] + dks2[2, 1] # dks2[]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 5560::

sage: Sym = SymmetricFunctions(QQ)
sage: Sym3 = Sym.kBoundedSubspace(3,t=1)
sage: Kks3 = Sym3.K_kschur()
sage: s = Sym.s()
sage: m = Sym.m()
sage: s(Kks3[3,1])
s[3] + s[3, 1] + s[4]
sage: m(Kks3[3,1])
m[1, 1, 1] + 4*m[1, 1, 1, 1] + m[2, 1] + 3*m[2, 1, 1] + 2*m[2, 2]
+ m[3] + 2*m[3, 1] + m[4]
sage: ks3 = Sym3.kschur()
sage: ks3(Kks3[3,1])
ks3[3] + ks3[3, 1]
sage: Kks3[3,1]*Kks3[2]            # long time (11.85 s, 2013)
-Kks3[3, 1, 1] - Kks3[3, 2] + Kks3[3, 2, 1] + Kks3[3, 3]
sage: Kks3[3,1].coproduct()
Kks3[] # Kks3[3, 1] - Kks3[1] # Kks3[2] + Kks3[1] # Kks3[2, 1]
+ 2*Kks3[1] # Kks3[3] + Kks3[1, 1] # Kks3[2] - Kks3[2] # Kks3[1]
+ Kks3[2] # Kks3[1, 1] + 2*Kks3[2] # Kks3[2] + Kks3[2, 1] # Kks3[1]
+ 2*Kks3[3] # Kks3[1] + Kks3[3, 1] # Kks3[]

Sage example in ./kschurnotes/notes-mike-anne.tex, line 5588::

sage: SymQ3 = Sym.kBoundedQuotient(3,t=1)
sage: G1 = SymQ3.AffineGrothendieckPolynomial([1],6)
sage: G2 = SymQ3.AffineGrothendieckPolynomial([2],6)
sage: (G1*G2).lift().scalar(Kks3[3,1])
-1
"""
