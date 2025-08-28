# sage.doctest: needs sage.combinat sage.modules
r"""
Tests For Combinatorial Hopf Algebras

We check that the composition of `Sym \to FSym^* \to QSym`
agrees with `Sym \to QSym`::

    sage: s = SymmetricFunctions(QQ).s()
    sage: M = QuasiSymmetricFunctions(QQ).M()
    sage: F = algebras.FSym(QQ).dual().F()
    sage: all(F(s[p]).to_quasisymmetric_function() == M(s[p])
    ....:     for n in range(5)
    ....:     for p in Partitions(n))
    True

There are various morphisms between the Hopf algebras below.
We test that the diagram of morphisms is commutative::

    sage: NSym = NonCommutativeSymmetricFunctions(QQ)
    sage: QSym = QuasiSymmetricFunctions(QQ)
    sage: FQSym = algebras.FQSym(QQ)
    sage: FSym = algebras.FSym(QQ)
    sage: FSymDual = FSym.dual()
    sage: Sym = SymmetricFunctions(QQ)

    sage: def go(composition):
    ....:     x = NSym.a_realization()[composition]
    ....:     if QSym(FQSym(x)) != QSym(Sym(x)): return False
    ....:     if Sym(FSym(x)) != Sym(x): return False
    ....:     if FQSym(FSym(x)) != FQSym(x): return False
    ....:     if FSymDual(Sym(x)) != FSymDual(FQSym(x)): return False
    ....:     return True

    sage: go([2,1,2])  # not tested (needs more morphisms)
    True
    sage: all(all(go(comp) for comp in Compositions(n)) for n in range(5))  # not tested (needs more morphisms)
    True

    sage: def go2(n):
    ....:     for sigma in Permutations(n):
    ....:          x = FQSym.F()[sigma]
    ....:          if QSym(FSymDual(x)) != QSym(x): return False
    ....:     s = Sym.s()
    ....:     for mu in Partitions(n):
    ....:          x = s[mu]
    ....:          if QSym(FSymDual(x)) != QSym(x): return False
    ....:     return True

    sage: all(go2(n) for n in range(6))  # not tested (needs more morphisms)
    True
"""
