# sage_setup: distribution = sagemath-repl
# sage.doctest: needs sage.libs.flint sage.modular
"""
This file contains a bunch of tests extracted from the published book
'Modular Forms: a Computational Approach' by William Stein, AMS 2007.

TESTS::

    sage: # needs sage.libs.gap
    sage: G = SL(2,ZZ); G
    Special Linear Group of degree 2 over Integer Ring
    sage: S, T = G.gens()
    sage: S
    [ 0  1]
    [-1  0]
    sage: T
    [1 1]
    [0 1]

    sage: delta_qexp(6)
    q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 + O(q^6)
    sage: bernoulli(12)
    -691/2730
    sage: bernoulli(50)
    495057205241079648212477525/66
    sage: len(str(bernoulli(10000)))
    27706
    sage: E4 = eisenstein_series_qexp(4, 3)
    sage: E6 = eisenstein_series_qexp(6, 3)
    sage: E4^6
    1/191102976000000 + 1/132710400000*q + 203/44236800000*q^2 + O(q^3)
    sage: E4^3*E6^2
    1/3511517184000 - 1/12192768000*q - 377/4064256000*q^2 + O(q^3)
    sage: E6^4
    1/64524128256 - 1/32006016*q + 241/10668672*q^2 + O(q^3)
    sage: victor_miller_basis(28,5)
    [
    1 + 15590400*q^3 + 36957286800*q^4 + O(q^5),
    q + 151740*q^3 + 61032448*q^4 + O(q^5),
    q^2 + 192*q^3 - 8280*q^4 + O(q^5)
    ]
    sage: R.<q> = QQ[['q']]
    sage: F4 =  240 * eisenstein_series_qexp(4,3)
    sage: F6 = -504 * eisenstein_series_qexp(6,3)
    sage: F4^3
    1 + 720*q + 179280*q^2 + O(q^3)
    sage: Delta = (F4^3 - F6^2)/1728; Delta
    q - 24*q^2 + O(q^3)
    sage: F4^3 - 720*Delta
    1 + 196560*q^2 + O(q^3)
    sage: M = ModularForms(1,36, prec=6).echelon_form()
    sage: M.basis()
    [
    1 + 6218175600*q^4 + 15281788354560*q^5 + O(q^6),
    q + 57093088*q^4 + 37927345230*q^5 + O(q^6),
    q^2 + 194184*q^4 + 7442432*q^5 + O(q^6),
    q^3 - 72*q^4 + 2484*q^5 + O(q^6)
    ]
    sage: T2 = M.hecke_matrix(2); T2
    [     34359738369                0       6218175600 9026867482214400]
    [               0                0      34416831456    5681332472832]
    [               0                1           194184       -197264484]
    [               0                0              -72           -54528]
    sage: T2.charpoly().factor()
    (x - 34359738369) * (x^3 - 139656*x^2 - 59208339456*x - 1467625047588864)
    sage: bernoulli_mod_p(23)
    [1, 4, 13, 17, 13, 6, 10, 5, 10, 9, 15]
    sage: set_modsym_print_mode ('modular')
    sage: M = ModularSymbols(11, 2)
    sage: M.basis()
    ({Infinity, 0}, {-1/8, 0}, {-1/9, 0})
    sage: S = M.cuspidal_submodule()
    sage: S.integral_basis()     # basis over ZZ.
    ({-1/8, 0}, {-1/9, 0})
    sage: set_modsym_print_mode ('manin')    # set it back
    sage: continued_fraction(4/7).convergents()
    [0, 1, 1/2, 4/7]
    sage: M = ModularSymbols(2,2)
    sage: M
    Modular Symbols space of dimension 1 for Gamma_0(2)
    of weight 2 with sign 0 over Rational Field
    sage: M.manin_generators()
    [(0,1), (1,0), (1,1)]

    sage: M = ModularSymbols(3,2)
    sage: M.manin_generators()
    [(0,1), (1,0), (1,1), (1,2)]

    sage: M = ModularSymbols(6,2)
    sage: M.manin_generators()
    [(0,1), (1,0), (1,1), (1,2), (1,3), (1,4), (1,5), (2,1),
     (2,3), (2,5), (3,1), (3,2)]
    sage: M = ModularSymbols(2,2)
    sage: [x.lift_to_sl2z(2) for x in M.manin_generators()]
    [[1, 0, 0, 1], [0, -1, 1, 0], [1, 0, 1, 1]]
    sage: M = ModularSymbols(6,2)
    sage: x = M.manin_generators()[9]
    sage: x
    (2,5)
    sage: x.lift_to_sl2z(6)
    [1, 2, 2, 5]
    sage: M = ModularSymbols(2,2)
    sage: M.manin_basis()
    [1]
    sage: [M.manin_generators()[i] for i in M.manin_basis()]
    [(1,0)]
    sage: M = ModularSymbols(6,2)
    sage: M.manin_basis()
    [1, 10, 11]
    sage: [M.manin_generators()[i] for i in M.manin_basis()]
    [(1,0), (3,1), (3,2)]
    sage: M.basis()
    ((1,0), (3,1), (3,2))
    sage: [x.modular_symbol_rep() for x in M.basis()]
    [{Infinity, 0}, {0, 1/3}, {-1/2, -1/3}]
    sage: M = ModularSymbols(2,2)
    sage: M.manin_gens_to_basis()
    [-1]
    [ 1]
    [ 0]
    sage: M = ModularSymbols(2,2)
    sage: x = (1,0); M(x)
    (1,0)
    sage: M( (3,1) )    # entries are reduced modulo 2 first
    0
    sage: M( (10,19) )
    -(1,0)
    sage: M = ModularSymbols(6,2)
    sage: M.manin_gens_to_basis()
    [-1  0  0]
    [ 1  0  0]
    [ 0  0  0]
    [ 0 -1  1]
    [ 0 -1  0]
    [ 0 -1  1]
    [ 0  0  0]
    [ 0  1 -1]
    [ 0  0 -1]
    [ 0  1 -1]
    [ 0  1  0]
    [ 0  0  1]
    sage: M = ModularSymbols(6,2)
    sage: M((0,1))
    -(1,0)
    sage: M((1,2))
    -(3,1) + (3,2)
    sage: HeilbronnCremona(2).to_list()
    [[1, 0, 0, 2], [2, 0, 0, 1], [2, 1, 0, 1], [1, 0, 1, 2]]
    sage: HeilbronnCremona(3).to_list()
    [[1, 0, 0, 3], [3, 1, 0, 1], [1, 0, 1, 3], [3, 0, 0, 1],
     [3, -1, 0, 1], [-1, 0, 1, -3]]
    sage: HeilbronnCremona(5).to_list()
    [[1, 0, 0, 5], [5, 2, 0, 1], [2, 1, 1, 3], [1, 0, 3, 5],
     [5, 1, 0, 1], [1, 0, 1, 5], [5, 0, 0, 1], [5, -1, 0, 1],
     [-1, 0, 1, -5], [5, -2, 0, 1], [-2, 1, 1, -3],
     [1, 0, -3, 5]]
    sage: len(HeilbronnCremona(37))
    128
    sage: len(HeilbronnCremona(389))
    1892
    sage: len(HeilbronnCremona(2003))
    11662
    sage: M = ModularSymbols(2,2)
    sage: M.T(2).matrix()
    [1]
    sage: M = ModularSymbols(6, 2)
    sage: M.T(2).matrix()
    [ 2  1 -1]
    [-1  0  1]
    [-1 -1  2]
    sage: M.T(3).matrix()
    [3 2 0]
    [0 1 0]
    [2 2 1]
    sage: M.T(3).fcp()  # factored characteristic polynomial
    (x - 3) * (x - 1)^2
    sage: M = ModularSymbols(39, 2)
    sage: T2 = M.T(2)
    sage: T2.matrix()
    [ 3  0 -1  0  0  1  1 -1  0]
    [ 0  0  2  0 -1  1  0  1 -1]
    [ 0  1  0 -1  1  1  0  1 -1]
    [ 0  0  1  0  0  1  0  1 -1]
    [ 0 -1  2  0  0  1  0  1 -1]
    [ 0  0  1  1  0  1  1 -1  0]
    [ 0  0  0 -1  0  1  1  2  0]
    [ 0  0  0  1  0  0  2  0  1]
    [ 0  0 -1  0  0  0  1  0  2]
    sage: T2.fcp()     # factored characteristic polynomial
    (x - 1)^2 * (x - 3)^3 * (x^2 + 2*x - 1)^2
    sage: T2 = M.T(2).matrix()
    sage: T5 = M.T(5).matrix()
    sage: T2*T5 - T5*T2 == 0
    True
    sage: T5.charpoly().factor()
    (x - 2)^2 * (x - 6)^3 * (x^2 - 8)^2
    sage: M = ModularSymbols(39, 2)
    sage: M.T(2).decomposition()
    [
    Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 9 for Gamma_0(39) of weight 2 with sign 0 over Rational Field,
    Modular Symbols subspace of dimension 3 of Modular Symbols space of dimension 9 for Gamma_0(39) of weight 2 with sign 0 over Rational Field,
    Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 9 for Gamma_0(39) of weight 2 with sign 0 over Rational Field
    ]
    sage: M = ModularSymbols(2, 2)
    sage: M.boundary_map()
    Hecke module morphism boundary map defined by the matrix
    [ 1 -1]
    Domain: Modular Symbols space of dimension 1 for
    Gamma_0(2) of weight ...
    Codomain: Space of Boundary Modular Symbols for
    Congruence Subgroup Gamma0(2) ...
    sage: M.cuspidal_submodule()
    Modular Symbols subspace of dimension 0 of Modular
    Symbols space of dimension 1 for Gamma_0(2) of weight
    2 with sign 0 over Rational Field
    sage: M = ModularSymbols(11, 2)
    sage: M.boundary_map().matrix()
    [ 1 -1]
    [ 0  0]
    [ 0  0]
    sage: M.cuspidal_submodule()
    Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
    sage: S = M.cuspidal_submodule(); S
    Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
    sage: S.basis()
    ((1,8), (1,9))
    sage: S.T(2).matrix()
    [-2  0]
    [ 0 -2]
    sage: S.T(3).matrix()
    [-1  0]
    [ 0 -1]
    sage: S.T(5).matrix()
    [1 0]
    [0 1]
    sage: E = EllipticCurve([0,-1,1,-10,-20])
    sage: 2 + 1 - E.Np(2)
    -2
    sage: 3 + 1 - E.Np(3)
    -1
    sage: 5 + 1 - E.Np(5)
    1
    sage: 7 + 1 - E.Np(7)
    -2
    sage: [S.T(p).matrix()[0,0] for p in [2,3,5,7]]
    [-2, -1, 1, -2]
    sage: M = ModularSymbols(11); M.basis()
    ((1,0), (1,8), (1,9))
    sage: S = M.cuspidal_submodule(); S
    Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
    sage: S.T(2).matrix()
    [-2  0]
    [ 0 -2]
    sage: S.T(3).matrix()
    [-1  0]
    [ 0 -1]
    sage: M = ModularSymbols(33)
    sage: S = M.cuspidal_submodule(); S
    Modular Symbols subspace of dimension 6 of Modular
    Symbols space of dimension 9 for Gamma_0(33) of weight
    2 with sign 0 over Rational Field
    sage: R.<q> = PowerSeriesRing(QQ)
    sage: v = [S.T(n).matrix()[0,0] for n in range(1,6)]
    sage: f00 = sum(v[n-1]*q^n for n in range(1,6)) + O(q^6)
    sage: f00
    q - q^2 - q^3 + q^4 + O(q^6)
    sage: v = [S.T(n).matrix()[0,1] for n in range(1,6)]
    sage: f01 = sum(v[n-1]*q^n for n in range(1,6)) + O(q^6)
    sage: f01
    -2*q^3 + O(q^6)
    sage: v = [S.T(n).matrix()[1,0] for n in range(1,6)]
    sage: f10 = sum(v[n-1]*q^n for n in range(1,6)) + O(q^6)
    sage: f10
    q^3 + O(q^6)
    sage: v = [S.T(n).matrix()[1,1] for n in range(1,6)]
    sage: f11 = sum(v[n-1]*q^n for n in range(1,6)) + O(q^6)
    sage: f11
    q - 2*q^2 + 2*q^4 + q^5 + O(q^6)
    sage: M = ModularSymbols(23)
    sage: S = M.cuspidal_submodule()
    sage: S
    Modular Symbols subspace of dimension 4 of Modular
    Symbols space of dimension 5 for Gamma_0(23) of weight
    2 with sign 0 over Rational Field
    sage: f = S.q_expansion_cuspforms(6)
    sage: f(0,0)
    q - 2/3*q^2 + 1/3*q^3 - 1/3*q^4 - 4/3*q^5 + O(q^6)
    sage: f(0,1)
    O(q^6)
    sage: f(1,0)
    -1/3*q^2 + 2/3*q^3 + 1/3*q^4 - 2/3*q^5 + O(q^6)
    sage: S.q_expansion_basis(6)
    [
    q - q^3 - q^4 + O(q^6),
    q^2 - 2*q^3 - q^4 + 2*q^5 + O(q^6)
    ]
    sage: R = Integers(49)
    sage: R
    Ring of integers modulo 49
    sage: R.unit_gens()
    (3,)
    sage: Integers(25).unit_gens()
    (2,)
    sage: Integers(100).unit_gens()
    (51, 77)
    sage: Integers(200).unit_gens()
    (151, 101, 177)
    sage: Integers(2005).unit_gens()
    (402, 1206)
    sage: Integers(200000000).unit_gens()
    (174218751, 51562501, 187109377)
    sage: list(DirichletGroup(5))
    [Dirichlet character modulo 5 of conductor 1 mapping 2 |--> 1,
    Dirichlet character modulo 5 of conductor 5 mapping 2 |--> zeta4,
    Dirichlet character modulo 5 of conductor 5 mapping 2 |--> -1,
    Dirichlet character modulo 5 of conductor 5 mapping 2 |--> -zeta4]
    sage: list(DirichletGroup(5, QQ))
    [Dirichlet character modulo 5 of conductor 1 mapping 2 |--> 1, Dirichlet character modulo 5 of conductor 5 mapping 2 |--> -1]
    sage: G = DirichletGroup(200)
    sage: G
    Group of Dirichlet characters modulo 200 with values in Cyclotomic Field of order 20 and degree 8
    sage: G.exponent()
    20
    sage: G.gens()
    (Dirichlet character modulo 200 of conductor 4 mapping 151 |--> -1, 101 |--> 1, 177 |--> 1,
    Dirichlet character modulo 200 of conductor 8 mapping 151 |--> 1, 101 |--> -1, 177 |--> 1,
    Dirichlet character modulo 200 of conductor 25 mapping 151 |--> 1, 101 |--> 1, 177 |--> zeta20)
    sage: K = G.base_ring()
    sage: zeta = K.0
    sage: eps = G([1,-1,zeta^5])
    sage: eps
    Dirichlet character modulo 200 of conductor 40 mapping 151 |--> 1, 101 |--> -1, 177 |--> zeta20^5
    sage: eps(3)
    zeta20^5
    sage: -zeta^15
    zeta20^5
    sage: kronecker(151,200)
    1
    sage: kronecker(101,200)
    -1
    sage: kronecker(177,200)
    1
    sage: G = DirichletGroup(20)
    sage: G.galois_orbits()
    [
    [Dirichlet character modulo 20 of conductor 20 mapping 11 |--> -1, 17 |--> -1],
    [Dirichlet character modulo 20 of conductor 20 mapping 11 |--> -1, 17 |--> -zeta4,
     Dirichlet character modulo 20 of conductor 20 mapping 11 |--> -1, 17 |--> zeta4],
    [Dirichlet character modulo 20 of conductor 4 mapping 11 |--> -1, 17 |--> 1],
    [Dirichlet character modulo 20 of conductor 5 mapping 11 |--> 1, 17 |--> -1],
    [Dirichlet character modulo 20 of conductor 5 mapping 11 |--> 1, 17 |--> -zeta4,
     Dirichlet character modulo 20 of conductor 5 mapping 11 |--> 1, 17 |--> zeta4],
    [Dirichlet character modulo 20 of conductor 1 mapping 11 |--> 1, 17 |--> 1]
    ]
    sage: G = DirichletGroup(11, QQ); G
    Group of Dirichlet characters modulo 11 with values in Rational Field
    sage: list(G)
    [Dirichlet character modulo 11 of conductor 1 mapping 2 |--> 1,
    Dirichlet character modulo 11 of conductor 11 mapping 2 |--> -1]
    sage: eps = G.0      # 0th generator for Dirichlet group
    sage: eps
    Dirichlet character modulo 11 of conductor 11 mapping 2 |--> -1
    sage: G.unit_gens()
    (2,)
    sage: eps(2)
    -1
    sage: eps(3)
    1
    sage: [eps(11*n) for n in range(10)]
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    sage: R = CyclotomicField(4)
    sage: CyclotomicField(4)
    Cyclotomic Field of order 4 and degree 2
    sage: G = DirichletGroup(15, R)
    sage: G
    Group of Dirichlet characters modulo 15 with values in Cyclotomic Field of order 4 and degree 2
    sage: list(G)
    [Dirichlet character modulo 15 of conductor 1 mapping 11 |--> 1, 7 |--> 1,
    Dirichlet character modulo 15 of conductor 3 mapping 11 |--> -1, 7 |--> 1,
    Dirichlet character modulo 15 of conductor 5 mapping 11 |--> 1, 7 |--> zeta4,
    Dirichlet character modulo 15 of conductor 15 mapping 11 |--> -1, 7 |--> zeta4,
    Dirichlet character modulo 15 of conductor 5 mapping 11 |--> 1, 7 |--> -1,
    Dirichlet character modulo 15 of conductor 15 mapping 11 |--> -1, 7 |--> -1,
    Dirichlet character modulo 15 of conductor 5 mapping 11 |--> 1, 7 |--> -zeta4,
    Dirichlet character modulo 15 of conductor 15 mapping 11 |--> -1, 7 |--> -zeta4]
    sage: e = G.1
    sage: e(4)
    -1
    sage: e(-1)
    -1
    sage: e(5)
    0
    sage: [e(n) for n in range(15)]
    [0, 1, zeta4, 0, -1, 0, 0, zeta4, -zeta4,
        0, 0, 1, 0, -zeta4, -1]
    sage: G = DirichletGroup(15, GF(5)); G
    Group of Dirichlet characters modulo 15 with values in Finite Field of size 5
    sage: list(G)
    [Dirichlet character modulo 15 of conductor 1 mapping 11 |--> 1, 7 |--> 1,
    Dirichlet character modulo 15 of conductor 3 mapping 11 |--> 4, 7 |--> 1,
    Dirichlet character modulo 15 of conductor 5 mapping 11 |--> 1, 7 |--> 2,
    Dirichlet character modulo 15 of conductor 15 mapping 11 |--> 4, 7 |--> 2,
    Dirichlet character modulo 15 of conductor 5 mapping 11 |--> 1, 7 |--> 4,
    Dirichlet character modulo 15 of conductor 15 mapping 11 |--> 4, 7 |--> 4,
    Dirichlet character modulo 15 of conductor 5 mapping 11 |--> 1, 7 |--> 3,
    Dirichlet character modulo 15 of conductor 15 mapping 11 |--> 4, 7 |--> 3]
    sage: e = G.1
    sage: e(-1)
    4
    sage: e(2)
    2
    sage: e(5)
    0
    sage: [e(n) for n in range(15)]
    [0, 1, 2, 0, 4, 0, 0, 2, 3, 0, 0, 1, 0, 3, 4]
    sage: G = DirichletGroup(5)
    sage: e = G.0
    sage: e(2)
    zeta4
    sage: e.bernoulli(1)
    -1/5*zeta4 - 3/5
    sage: e.bernoulli(9)
    -108846/5*zeta4 - 176868/5
    sage: E = EisensteinForms(Gamma1(13),2)
    sage: E.eisenstein_series()
    [
    1/2 + q + 3*q^2 + 4*q^3 + 7*q^4 + 6*q^5 + O(q^6),
    -7/13*zeta6 - 11/13 + q + (2*zeta6 + 1)*q^2 + (-3*zeta6 + 1)*q^3 + (6*zeta6 - 3)*q^4 - 4*q^5 + O(q^6),
    q + (zeta6 + 2)*q^2 + (-zeta6 + 3)*q^3 + (3*zeta6 + 3)*q^4 + 4*q^5 + O(q^6),
    -zeta6 + q + (2*zeta6 - 1)*q^2 + (3*zeta6 - 2)*q^3 + (-2*zeta6 - 1)*q^4 + 6*q^5 + O(q^6),
    q + (zeta6 + 1)*q^2 + (zeta6 + 2)*q^3 + (zeta6 + 2)*q^4 + 6*q^5 + O(q^6),
    -1 + q - q^2 + 4*q^3 + 3*q^4 - 4*q^5 + O(q^6),
    q + q^2 + 4*q^3 + 3*q^4 + 4*q^5 + O(q^6),
    zeta6 - 1 + q + (-2*zeta6 + 1)*q^2 + (-3*zeta6 + 1)*q^3 + (2*zeta6 - 3)*q^4 + 6*q^5 + O(q^6),
    q + (-zeta6 + 2)*q^2 + (-zeta6 + 3)*q^3 + (-zeta6 + 3)*q^4 + 6*q^5 + O(q^6),
    7/13*zeta6 - 18/13 + q + (-2*zeta6 + 3)*q^2 + (3*zeta6 - 2)*q^3 + (-6*zeta6 + 3)*q^4 - 4*q^5 + O(q^6),
    q + (-zeta6 + 3)*q^2 + (zeta6 + 2)*q^3 + (-3*zeta6 + 6)*q^4 + 4*q^5 + O(q^6)
    ]
    sage: e = E.eisenstein_series()
    sage: for e in E.eisenstein_series():
    ....:     print(e.parameters())
    (Dirichlet character modulo 13 of conductor 1 mapping 2 |--> 1, Dirichlet character modulo 13 of conductor 1 mapping 2 |--> 1, 13)
    (Dirichlet character modulo 13 of conductor 1 mapping 2 |--> 1, Dirichlet character modulo 13 of conductor 13 mapping 2 |--> zeta6, 1)
    (Dirichlet character modulo 13 of conductor 13 mapping 2 |--> zeta6, Dirichlet character modulo 13 of conductor 1 mapping 2 |--> 1, 1)
    (Dirichlet character modulo 13 of conductor 1 mapping 2 |--> 1, Dirichlet character modulo 13 of conductor 13 mapping 2 |--> zeta6 - 1, 1)
    (Dirichlet character modulo 13 of conductor 13 mapping 2 |--> zeta6 - 1, Dirichlet character modulo 13 of conductor 1 mapping 2 |--> 1, 1)
    (Dirichlet character modulo 13 of conductor 1 mapping 2 |--> 1, Dirichlet character modulo 13 of conductor 13 mapping 2 |--> -1, 1)
    (Dirichlet character modulo 13 of conductor 13 mapping 2 |--> -1, Dirichlet character modulo 13 of conductor 1 mapping 2 |--> 1, 1)
    (Dirichlet character modulo 13 of conductor 1 mapping 2 |--> 1, Dirichlet character modulo 13 of conductor 13 mapping 2 |--> -zeta6, 1)
    (Dirichlet character modulo 13 of conductor 13 mapping 2 |--> -zeta6, Dirichlet character modulo 13 of conductor 1 mapping 2 |--> 1, 1)
    (Dirichlet character modulo 13 of conductor 1 mapping 2 |--> 1, Dirichlet character modulo 13 of conductor 13 mapping 2 |--> -zeta6 + 1, 1)
    (Dirichlet character modulo 13 of conductor 13 mapping 2 |--> -zeta6 + 1, Dirichlet character modulo 13 of conductor 1 mapping 2 |--> 1, 1)
    sage: from sage.modular.dims import *
    sage: dimension_cusp_forms(Gamma0(2007),2)
    221
    sage: dimension_eis(Gamma0(2007),2)
    7
    sage: dimension_modular_forms(Gamma0(2007),2)
    228
    sage: dimension_new_cusp_forms(Gamma0(11),12)
    8
    sage: dimension_cusp_forms(Gamma0(11),12)
    10
    sage: dimension_new_cusp_forms(Gamma0(2007),12)
    1017
    sage: dimension_cusp_forms(Gamma0(2007),12)
    2460
    sage: dimension_cusp_forms(Gamma1(2007),2)
    147409
    sage: dimension_eis(Gamma1(2007),2)
    3551
    sage: dimension_modular_forms(Gamma1(2007),2)
    150960
    sage: G = DirichletGroup(2007)
    sage: e = prod(G.gens(), G(1))
    sage: dimension_cusp_forms(e,2)
    222
    sage: dimension_cusp_forms(e,3)
    0
    sage: dimension_cusp_forms(e,4)
    670
    sage: dimension_cusp_forms(e,24)
    5150
    sage: dimension_eis(e,2)
    4
    sage: dimension_eis(e,3)
    0
    sage: dimension_eis(e,4)
    4
    sage: dimension_eis(e,24)
    4
    sage: G = DirichletGroup(2007, QQ)
    sage: e = prod(G.gens(), G(1))
    sage: dimension_new_cusp_forms(e,2)
    76
    sage: p = 389
    sage: k = GF(p)
    sage: a = k(7/13); a
    210
    sage: a.rational_reconstruction()
    7/13
    sage: M = MatrixSpace(QQ,4,8)
    sage: A = M([[-9,6,7,3,1,0,0,0],[-10,3,8,2,0,1,0,0],[3,-6,2,8,0,0,1,0],[-8,-6,-8,6,0,0,0,1]])
    sage: A41 = MatrixSpace(GF(41),4,8)(A)
    sage: E41 = A41.echelon_form()
    sage: B = A.matrix_from_columns([0,1,2,4])
    sage: E = B^(-1)*A
    sage: B = A.matrix_from_columns([0,1,2,3])
    sage: M = ModularSymbols(Gamma0(6)); M
    Modular Symbols space of dimension 3 for Gamma_0(6)
    of weight 2 with sign 0 over Rational Field
    sage: M.new_subspace()
    Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 3 for Gamma_0(6) of weight 2 with sign 0 over Rational Field
    sage: M.old_subspace()
    Modular Symbols subspace of dimension 3 of Modular
    Symbols space of dimension 3 for Gamma_0(6) of weight
    2 with sign 0 over Rational Field
    sage: G = DirichletGroup(13)
    sage: G
    Group of Dirichlet characters modulo 13 with values in Cyclotomic Field of order 12 and degree 4
    sage: dimension_modular_forms(Gamma1(13),2)
    13
    sage: [dimension_modular_forms(e,2) for e in G]
    [1, 0, 3, 0, 2, 0, 2, 0, 2, 0, 3, 0]
    sage: G = DirichletGroup(100)
    sage: G
    Group of Dirichlet characters modulo 100 with values in Cyclotomic Field of order 20 and degree 8
    sage: dimension_modular_forms(Gamma1(100),2)
    370
    sage: v = [dimension_modular_forms(e,2) for e in G]; v
    [24, 0, 0, 17, 18, 0, 0, 17, 18, 0, 0, 21, 18, 0, 0, 17,
      18, 0, 0, 17, 24, 0, 0, 17, 18, 0, 0, 17, 18, 0, 0, 21,
      18, 0, 0, 17, 18, 0, 0, 17]
    sage: sum(v)
    370
    sage: S = CuspForms(Gamma0(45), 2, prec=14); S
    Cuspidal subspace of dimension 3 of Modular Forms space
    of dimension 10 for Congruence Subgroup Gamma0(45) of
    weight 2 over Rational Field
    sage: S.basis()
    [
    q - q^4 - q^10 - 2*q^13 + O(q^14),
    q^2 - q^5 - 3*q^8 + 4*q^11 + O(q^14),
    q^3 - q^6 - q^9 - q^12 + O(q^14)
    ]
    sage: S.new_subspace().basis()
    [
    q + q^2 - q^4 - q^5 - 3*q^8 - q^10 + 4*q^11 - 2*q^13 + O(q^14)
    ]
    sage: CuspForms(Gamma0(9),2)
    Cuspidal subspace of dimension 0 of Modular Forms space
    of dimension 3 for Congruence Subgroup Gamma0(9) of
    weight 2 over Rational Field
    sage: CuspForms(Gamma0(15),2, prec=10).basis()
    [
    q - q^2 - q^3 - q^4 + q^5 + q^6 + 3*q^8 + q^9 + O(q^10)
    ]
"""
