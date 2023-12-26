# sage_setup: distribution = sagemath-combinat
r"""
Combinatorics quickref

Integer Sequences::

    sage: s = oeis([1,3,19,211]); s                         # optional - internet
    0: A000275: Coefficients of a Bessel function (reciprocal of J_0(z));
                also pairs of permutations with rise/rise forbidden.
    sage: s[0].programs()                                   # optional - internet
    [('maple', ...),
     ('mathematica', ...),
     ('pari',
      0: {a(n) = if( n<0, 0, n!^2 * 4^n * polcoeff( 1 / besselj(0, x + x * O(x^(2*n))), 2*n))}; /* _Michael Somos_, May 17 2004 */)]

Combinatorial objects::

    sage: S = Subsets([1,2,3,4]); S.list(); S.<tab>                       # not tested
    sage: P = Partitions(10000); P.cardinality()                                        # needs sage.libs.flint
    3616...315650422081868605887952568754066420592310556052906916435144
    sage: Combinations([1,3,7]).random_element()            # random
    sage: Compositions(5, max_part=3).unrank(3)
    [2, 2, 1]

    sage: DyckWord([1,0,1,0,1,1,0,0]).to_binary_tree()                                  # needs sage.graphs
    [., [., [[., .], .]]]
    sage: Permutation([3,1,4,2]).robinson_schensted()
    [[[1, 2], [3, 4]], [[1, 3], [2, 4]]]
    sage: StandardTableau([[1, 4], [2, 5], [3]]).schuetzenberger_involution()
    [[1, 3], [2, 4], [5]]

Constructions and Species::

    sage: for (p, s) in cartesian_product([P,S]): print((p, s)) # not tested
    sage: def IV_3(n):
    ....:     return IntegerVectors(n, 3)
    sage: DisjointUnionEnumeratedSets(Family(IV_3, NonNegativeIntegers))  # not tested

Words::

    sage: Words('abc', 4).list()
    [word: aaaa, ..., word: cccc]

    sage: Word('aabcacbaa').is_palindrome()
    True
    sage: WordMorphism('a->ab,b->a').fixed_point('a')
    word: abaababaabaababaababaabaababaabaababaaba...

Polytopes::

    sage: points = random_matrix(ZZ, 6, 3, x=7).rows()                                  # needs sage.modules
    sage: L = LatticePolytope(points)                                                   # needs sage.geometry.polyhedron sage.modules
    sage: L.npoints(); L.plot3d()                           # random                    # needs sage.geometry.polyhedron sage.modules sage.plot

:ref:`Root systems, Coxeter and Weyl groups <sage.combinat.root_system.all>`::

    sage: WeylGroup(["B",3]).bruhat_poset()                                             # needs sage.graphs sage.modules
    Finite poset containing 48 elements
    sage: RootSystem(["A",2,1]).weight_lattice().plot()         # not tested            # needs sage.graphs sage.modules sage.plot

:ref:`Crystals <sage.combinat.crystals.all>`::

    sage: CrystalOfTableaux(["A",3], shape=[3,2]).some_flashy_feature()   # not tested

:mod:`Symmetric functions and combinatorial Hopf algebras <sage.combinat.algebraic_combinatorics>`::

    sage: Sym = SymmetricFunctions(QQ); Sym.inject_shorthands(verbose=False)            # needs sage.sage.modules
    sage: m( ( h[2,1] * (1 + 3 * p[2,1]) ) + s[2](s[3]) )                               # needs sage.sage.modules
    3*m[1, 1, 1] + ... + 10*m[5, 1] + 4*m[6]

:ref:`Discrete groups, Permutation groups <sage.groups.groups_catalog>`::

    sage: S = SymmetricGroup(4)                                                         # needs sage.groups
    sage: M = PolynomialRing(QQ, 'x0,x1,x2,x3')
    sage: M.an_element() * S.an_element()                                               # needs sage.groups
    x0

Graph theory, posets, lattices (:ref:`sage.graphs`, :ref:`sage.combinat.posets.all`)::

    sage: Poset({1: [2,3], 2: [4], 3: [4]}).linear_extensions().cardinality()           # needs sage.graphs sage.modules
    2
"""
