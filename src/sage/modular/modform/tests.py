# sage_setup: distribution = sagemath-schemes
# sage.doctest: needs sage.libs.pari
"""
TESTS:

Test consistency between modular forms and modular symbols Eisenstein
series charpolys, which are computed in completely separate ways. ::

    sage: for N in range(25,33):
    ....:    m = ModularForms(N)
    ....:    e = m.eisenstein_subspace()
    ....:    f = e.hecke_polynomial(2)
    ....:    g = ModularSymbols(N).eisenstein_submodule().hecke_polynomial(2)
    ....:    print("{} {}".format(N, f == g))
    25 True
    26 True
    27 True
    28 True
    29 True
    30 True
    31 True
    32 True

Another similar consistency check::

    sage: for N in range(1,26):
    ....:    eps = (N > 2 and DirichletGroup(N,QQ).0) or N
    ....:    m = ModularForms(eps)
    ....:    e = m.eisenstein_subspace()
    ....:    f = e.hecke_polynomial(2)
    ....:    g = ModularSymbols(eps).eisenstein_submodule().hecke_polynomial(2)
    ....:    print("{} {}".format(N, f == g))
    1 True
    2 True
    3 True
    4 True
    5 True
    6 True
    7 True
    8 True
    9 True
    10 True
    11 True
    12 True
    13 True
    14 True
    15 True
    16 True
    17 True
    18 True
    19 True
    20 True
    21 True
    22 True
    23 True
    24 True
    25 True

We check that bug :issue:`8541` (traced to a linear algebra problem) is fixed::

    sage: f = CuspForms(DirichletGroup(5).0,5).0
    sage: f[15]
    30*zeta4 - 210
"""
