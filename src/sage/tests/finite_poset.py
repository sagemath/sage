# sage.doctest: needs sage.graphs sage.modules
"""
This file contains test functions that can be used to search
bugs by testing random finite posets and lattices.

As an examples: if a lattice is distributive, then it must be also
modular, and if a poset is ranked, then the dual poset must also
be ranked.
"""

from sage.misc.prandom import randint
from sage.misc.call import attrcall
from functools import reduce

implications = {
 'doubling_convex': ['doubling_any'],
 'doubling_interval': ['doubling_lower', 'doubling_upper'],
 'doubling_lower': ['doubling_convex', 'meet_semidistributive'],
 'doubling_upper': ['doubling_convex', 'join_semidistributive'],
 'cosectionally_complemented': ['complemented', 'coatomic', 'regular'],
 'distributive': ['modular', 'semidistributive', 'join_distributive', 'meet_distributive', 'subdirectly_reducible', 'doubling_interval', 'extremal'],
 'geometric': ['upper_semimodular', 'relatively_complemented'],
 'isoform': ['uniform'],
 'join_distributive': ['meet_semidistributive', 'upper_semimodular'],
 'join_semidistributive': ['join_pseudocomplemented', 'interval_dismantlable'],
 'lower_semimodular': ['graded'],
 'meet_distributive': ['join_semidistributive', 'lower_semimodular'],
 'meet_semidistributive': ['pseudocomplemented', 'interval_dismantlable'],
 'modular': ['upper_semimodular', 'lower_semimodular', 'supersolvable'],
 'orthocomplemented': ['self_dual', 'complemented'],
 'planar': ['dismantlable'],
 'dismantlable': ['sublattice_dismantlable'],
 'interval_dismantlable': ['sublattice_dismantlable'],
 'relatively_complemented': ['sectionally_complemented', 'cosectionally_complemented', 'isoform'],
 'sectionally_complemented': ['complemented', 'atomic', 'regular'],
 'semidistributive': ['join_semidistributive', 'meet_semidistributive'],
 'simple': ['isoform'],
 'supersolvable': ['graded'],
 'uniform': ['regular'],
 'uniq_orthocomplemented': ['orthocomplemented'],
 'upper_semimodular': ['graded'],
 'vertically_decomposable': ['subdirectly_reducible'],
}

dual_properties = [
 ['atomic', 'coatomic'],
 ['upper_semimodular', 'lower_semimodular'],
 ['sectionally_complemented', 'cosectionally_complemented'],
 ['join_distributive', 'meet_distributive'],
 ['join_semidistributive', 'meet_semidistributive'],
 ['pseudocomplemented', 'join_pseudocomplemented'],
 ['doubling_lower', 'doubling_upper'],
]

selfdual_properties = ['distributive', 'modular', 'semidistributive', 'complemented',
 'relatively_complemented', 'orthocomplemented', 'uniq_orthocomplemented', 'supersolvable', 'planar',
 'dismantlable', 'vertically_decomposable', 'simple', 'isoform', 'uniform', 'regular',
 'subdirectly_reducible', 'doubling_any', 'doubling_convex', 'doubling_interval',
 'interval_dismantlable', 'interval_dismantlable']

dual_elements = [
 ['atoms', 'coatoms'],
 ['meet_irreducibles', 'join_irreducibles'],
 ['meet_primes', 'join_primes']
]

two_to_one = [ ['distributive', 'dismantlable', 'planar'],
               ['upper_semimodular', 'lower_semimodular', 'modular'],
               ['meet_distributive', 'join_distributive', 'distributive'],
               ['meet_semidistributive', 'join_semidistributive', 'semidistributive'],
               ['lower_semimodular', 'meet_semidistributive', 'distributive'],
               ['upper_semimodular', 'join_semidistributive', 'distributive'],
               ['complemented', 'modular', 'relatively_complemented'],
             ]

mutually_exclusive = [
 ['doubling_any', 'simple'],
 ['vertically_decomposable', 'atomic'],
 ['vertically_decomposable', 'coatomic'],
 ['vertically_decomposable', 'regular'],
]

set_inclusions = [
 ['atoms', 'join_irreducibles'],
 ['coatoms', 'meet_irreducibles'],
 ['double_irreducibles', 'join_irreducibles'],
 ['double_irreducibles', 'meet_irreducibles'],
 ['meet_primes', 'meet_irreducibles'],
 ['join_primes', 'join_irreducibles'],
]

sublattice_closed = ['distributive', 'modular', 'semidistributive', 'join_semidistributive', 'meet_semidistributive']


def test_attrcall(name, L):
    """
    Return a function by name.

    This is a helper function for test_finite_lattice(). This
    will unify all Boolean-valued functions to a function without
    parameters.

    EXAMPLES::

        sage: from sage.tests.finite_poset import test_attrcall
        sage: N5 = posets.PentagonPoset()
        sage: N5.is_modular() == test_attrcall('is_modular', N5)
        True
        sage: N5.is_constructible_by_doublings('convex') == test_attrcall('is_doubling_convex', N5)                     # needs sage.combinat
        True
    """
    if name == 'is_doubling_any':
        return L.is_constructible_by_doublings('any')
    if name == 'is_doubling_lower':
        return L.is_constructible_by_doublings('upper')
    if name == 'is_doubling_upper':
        return L.is_constructible_by_doublings('lower')
    if name == 'is_doubling_convex':
        return L.is_constructible_by_doublings('convex')
    if name == 'is_doubling_interval':
        return L.is_constructible_by_doublings('interval')
    if name == 'is_uniq_orthocomplemented':
        return L.is_orthocomplemented(unique=True)
    return attrcall(name)(L)


def test_finite_lattice(L):
    """
    Test several functions on a given finite lattice.

    The function contains tests of different kinds:

    - Implications of Boolean properties. Examples: a distributive lattice is modular,
      a dismantlable and distributive lattice is planar, a simple lattice can not be
      constructible by Day's doublings.
    - Dual and self-dual properties. Examples: Dual of a modular lattice is modular,
      dual of an atomic lattice is co-atomic.
    - Certificate tests. Example: certificate for a non-complemented lattice must be
      an element without a complement.
    - Verification of some property by known property or by a random test.
      Examples: A lattice is distributive iff join-primes are exactly
      join-irreducibles and an interval of a relatively complemented
      lattice is complemented.
    - Set inclusions. Example: Every co-atom must be meet-irreducible.
    - And several other tests. Example: The skeleton of a pseudocomplemented
      lattice must be Boolean.

    EXAMPLES::

        sage: from sage.tests.finite_poset import test_finite_lattice
        sage: L = posets.RandomLattice(10, 0.98)
        sage: test_finite_lattice(L) is None  # Long time
        True
    """
    from sage.combinat.posets.lattices import LatticePoset

    from sage.sets.set import Set
    from sage.combinat.subset import Subsets

    from sage.misc.prandom import randint
    from sage.misc.flatten import flatten
    from sage.misc.call import attrcall

    from sage.misc.sageinspect import sage_getargspec

    if L.cardinality() < 4:
        # Special cases should be tested in specific TESTS-sections.
        return None

    all_props = set(list(implications) + flatten(implications.values()))
    P = {x: test_attrcall('is_' + x, L) for x in all_props}

    ### Relations between boolean-valued properties ###

    # Direct one-property implications
    for prop1 in implications:
        if P[prop1]:
            for prop2 in implications[prop1]:
                if not P[prop2]:
                    raise ValueError("error: %s should implicate %s" % (prop1, prop2))

    # Impossible combinations
    for p1, p2 in mutually_exclusive:
        if P[p1] and P[p2]:
            raise ValueError("error: %s and %s should be impossible combination" % (p1, p2))

    # Two-property implications
    for p1, p2, p3 in two_to_one:
        if P[p1] and P[p2] and not P[p3]:
            raise ValueError("error: %s and %s, so should be %s" % (p1, p2, p3))

    Ldual = L.dual()
    # Selfdual properties
    for p in selfdual_properties:
        if P[p] != test_attrcall('is_'+p, Ldual):
            raise ValueError("selfdual property %s error" % p)
    # Dual properties and elements
    for p1, p2 in dual_properties:
        if P[p1] != test_attrcall('is_'+p2, Ldual):
            raise ValueError("dual properties error %s" % p1)
    for e1, e2 in dual_elements:
        if set(attrcall(e1)(L)) != set(attrcall(e2)(Ldual)):
            raise ValueError("dual elements error %s" % e1)

    ### Certificates ###

    # Return value must be a pair with correct result as first element.
    for p_ in all_props:
        # Dirty fix first
        if p_[:9] == 'doubling_' or p_[:5] == 'uniq_':
            continue
        p = "is_"+p_
        if 'certificate' in sage_getargspec(getattr(L, p)).args:
            res = attrcall(p, certificate=True)(L)
            if not isinstance(res, tuple) or len(res) != 2:
                raise ValueError("certificate-option does not return a pair in %s" % p)
            if P[p_] != res[0]:
                raise ValueError("certificate-option changes result in %s" % p)

    # Test for "yes"-certificates
    if P['supersolvable']:
        a = L.is_supersolvable(certificate=True)[1]
        S = Subsets(L).random_element()
        if L.is_chain_of_poset(S):
            if not L.sublattice(a+list(S)).is_distributive():
                raise ValueError("certificate error in is_supersolvable")
    if P['dismantlable']:
        elms = L.is_dismantlable(certificate=True)[1]
        if len(elms) != L.cardinality():
            raise ValueError("certificate error 1 in is_dismantlable")
        elms = elms[:randint(0, len(elms)-1)]
        L_ = L.sublattice([x for x in L if x not in elms])
        if L_.cardinality() != L.cardinality() - len(elms):
            raise ValueError("certificate error 2 in is_dismantlable")
    if P['vertically_decomposable']:
        c = L.is_vertically_decomposable(certificate=True)[1]
        if c == L.bottom() or c == L.top():
            raise ValueError("certificate error 1 in is_vertically_decomposable")
        e = L.random_element()
        if L.compare_elements(c, e) is None:
            raise ValueError("certificate error 2 in is_vertically_decomposable")

    # Test for "no"-certificates
    if not P['atomic']:
        a = L.is_atomic(certificate=True)[1]
        if a in L.atoms() or a not in L.join_irreducibles():
            raise ValueError("certificate error in is_atomic")
    if not P['coatomic']:
        a = L.is_coatomic(certificate=True)[1]
        if a in L.coatoms() or a not in L.meet_irreducibles():
            raise ValueError("certificate error in is_coatomic")

    if not P['complemented']:
        a = L.is_complemented(certificate=True)[1]
        if L.complements(a):
            raise ValueError("compl. error 1")
    if not P['sectionally_complemented']:
        a, b = L.is_sectionally_complemented(certificate=True)[1]
        L_ = L.sublattice(L.interval(L.bottom(), a))
        if L_.is_complemented():
            raise ValueError("sec. compl. error 1")
        if len(L_.complements(b)) > 0:
            raise ValueError("sec. compl. error 2")
    if not P['cosectionally_complemented']:
        a, b = L.is_cosectionally_complemented(certificate=True)[1]
        L_ = L.sublattice(L.interval(a, L.top()))
        if L_.is_complemented():
            raise ValueError("cosec. compl. error 1")
        if L_.complements(b):
            raise ValueError("cosec. compl. error 2")
    if not P['relatively_complemented']:
        a, b, c = L.is_relatively_complemented(certificate=True)[1]
        I = L.interval(a, c)
        if len(I) != 3 or b not in I:
            raise ValueError("rel. compl. error 1")

    if not P['upper_semimodular']:
        a, b = L.is_upper_semimodular(certificate=True)[1]
        if not set(L.lower_covers(a)).intersection(set(L.lower_covers(b))) or set(L.upper_covers(a)).intersection(set(L.upper_covers(b))):
            raise ValueError("certificate error in is_upper_semimodular")
    if not P['lower_semimodular']:
        a, b = L.is_lower_semimodular(certificate=True)[1]
        if set(L.lower_covers(a)).intersection(set(L.lower_covers(b))) or not set(L.upper_covers(a)).intersection(set(L.upper_covers(b))):
            raise ValueError("certificate error in is_lower_semimodular")

    if not P['distributive']:
        x, y, z = L.is_distributive(certificate=True)[1]
        if L.meet(x, L.join(y, z)) == L.join(L.meet(x, y), L.meet(x, z)):
            raise ValueError("certificate error in is_distributive")
    if not P['modular']:
        x, a, b = L.is_modular(certificate=True)[1]
        if not L.is_less_than(x, b) or L.join(x, L.meet(a, b)) == L.meet(L.join(x, a), b):
            raise ValueError("certificate error in is_modular")

    if not P['pseudocomplemented']:
        a = L.is_pseudocomplemented(certificate=True)[1]
        L_ = L.subposet([e for e in L if L.meet(e, a) == L.bottom()])
        if L_.has_top():
            raise ValueError("certificate error in is_pseudocomplemented")
    if not P['join_pseudocomplemented']:
        a = L.is_join_pseudocomplemented(certificate=True)[1]
        L_ = L.subposet([e for e in L if L.join(e, a) == L.top()])
        if L_.has_bottom():
            raise ValueError("certificate error in is_join_pseudocomplemented")

    if not P['join_semidistributive']:
        e, x, y = L.is_join_semidistributive(certificate=True)[1]
        if L.join(e, x) != L.join(e, y) or L.join(e, x) == L.join(e, L.meet(x, y)):
            raise ValueError("certificate error in is_join_semidistributive")
    if not P['meet_semidistributive']:
        e, x, y = L.is_meet_semidistributive(certificate=True)[1]
        if L.meet(e, x) != L.meet(e, y) or L.meet(e, x) == L.meet(e, L.join(x, y)):
            raise ValueError("certificate error in is_meet_semidistributive")

    if not P['simple']:
        c = L.is_simple(certificate=True)[1]
        if len(L.congruence([c[randint(0, len(c)-1)]])) == 1:
            raise ValueError("certificate error in is_simple")
    if not P['isoform']:
        c = L.is_isoform(certificate=True)[1]
        if len(c) == 1:
            raise ValueError("certificate error in is_isoform")
        if all(L.subposet(c[i]).is_isomorphic(L.subposet(c[i+1])) for i in range(len(c)-1)):
            raise ValueError("certificate error in is_isoform")
    if not P['uniform']:
        c = L.is_uniform(certificate=True)[1]
        if len(c) == 1:
            raise ValueError("certificate error in is_uniform")
        if all(len(c[i]) == len(c[i+1]) for i in range(len(c)-1)):
            raise ValueError("certificate error in is_uniform")
    if not P['regular']:
        c = L.is_regular(certificate=True)[1]
        if len(c[0]) == 1:
            raise ValueError("certificate error 1 in is_regular")
        if Set(c[1]) not in c[0]:
            raise ValueError("certificate error 2 in is_regular")
        if L.congruence([c[1]]) == c[0]:
            raise ValueError("certificate error 3 in is_regular")

    if not P['subdirectly_reducible']:
        x, y = L.is_subdirectly_reducible(certificate=True)[1]
        a = L.random_element()
        b = L.random_element()
        c = L.congruence([[a, b]])
        if len(c) != L.cardinality():
            for c_ in c:
                if x in c_:
                    if y not in c_:
                        raise ValueError("certificate error 1 in is_subdirectly_reducible")
                    break
            else:
                raise ValueError("certificate error 2 in is_subdirectly_reducible")

    if not P['join_distributive']:
        a = L.is_join_distributive(certificate=True)[1]
        L_ = L.sublattice(L.interval(a, L.join(L.upper_covers(a))))
        if L_.is_distributive():
            raise ValueError("certificate error in is_join_distributive")
    if not P['meet_distributive']:
        a = L.is_meet_distributive(certificate=True)[1]
        L_ = L.sublattice(L.interval(L.meet(L.lower_covers(a)), a))
        if L_.is_distributive():
            raise ValueError("certificate error in is_meet_distributive")

    ### Other ###

    # Other ways to recognize some boolean property
    if P['distributive'] != (set(L.join_primes()) == set(L.join_irreducibles())):
        raise ValueError("every join-irreducible of a distributive lattice should be join-prime")
    if P['distributive'] != (set(L.meet_primes()) == set(L.meet_irreducibles())):
        raise ValueError("every meet-irreducible of a distributive lattice should be meet-prime")
    if P['join_semidistributive'] != all(L.canonical_joinands(e) is not None for e in L):
        raise ValueError("every element of join-semidistributive lattice should have canonical joinands")
    if P['meet_semidistributive'] != all(L.canonical_meetands(e) is not None for e in L):
        raise ValueError("every element of meet-semidistributive lattice should have canonical meetands")

    # Random verification of a Boolean property
    if P['relatively_complemented']:
        a = L.random_element()
        b = L.random_element()
        if not L.sublattice(L.interval(a, b)).is_complemented():
            raise ValueError("rel. compl. error 3")
    if P['sectionally_complemented']:
        a = L.random_element()
        if not L.sublattice(L.interval(L.bottom(), a)).is_complemented():
            raise ValueError("sec. compl. error 3")
    if P['cosectionally_complemented']:
        a = L.random_element()
        if not L.sublattice(L.interval(a, L.top())).is_complemented():
            raise ValueError("cosec. compl. error 2")

    # Element set inclusions
    for s1, s2 in set_inclusions:
        if not set(attrcall(s1)(L)).issubset(set(attrcall(s2)(L))):
            raise ValueError("%s should be a subset of %s" % (s1, s2))

    # Sublattice-closed properties
    L_ = L.sublattice(Subsets(L).random_element())
    for p in sublattice_closed:
        if P[p] and not test_attrcall('is_'+p, L_):
            raise ValueError("property %s should apply to sublattices" % p)

    # Some sublattices
    L_ = L.center()  # Center is a Boolean lattice
    if not L_.is_atomic() or not L_.is_distributive():
        raise ValueError("error in center")
    if P['pseudocomplemented']:
        L_ = L.skeleton()  # Skeleton is a Boolean lattice
        if not L_.is_atomic() or not L_.is_distributive():
            raise ValueError("error in skeleton")
    L_ = L.frattini_sublattice()
    S = Subsets(L).random_element()
    if L.sublattice(S) == L and L.sublattice([e for e in S if e not in L_]) != L:
        raise ValueError("error in Frattini sublattice")
    L_ = L.maximal_sublattices()
    L_ = L_[randint(0, len(L_)-1)]
    e = L.random_element()
    if e not in L_ and L.sublattice(list(L_)+[e]) != L:
        raise ValueError("error in maximal_sublattices")

    # Reverse functions: vertical composition and decomposition
    L_ = reduce(lambda a, b: a.vertical_composition(b), L.vertical_decomposition(), LatticePoset())
    if not L.is_isomorphic(L_):
        raise ValueError("error in vertical [de]composition")

    # Meet and join
    a = L.random_element()
    b = L.random_element()
    m = L.meet(a, b)
    j = L.join(a, b)
    m_ = L.subposet([e for e in L.principal_lower_set(a) if e in L.principal_lower_set(b)]).top()
    j_ = L.subposet([e for e in L.principal_upper_set(a) if e in L.principal_upper_set(b)]).bottom()
    if m != m_ or m != Ldual.join(a, b):
        raise ValueError("error in meet")
    if j != j_ or j != Ldual.meet(a, b):
        raise ValueError("error in join")

    # Misc misc
    e = L.neutral_elements()
    e = e[randint(0, len(e)-1)]
    a = L.random_element()
    b = L.random_element()
    if not L.sublattice([e, a, b]).is_distributive():
        raise ValueError("error in neutral_elements")


def test_finite_poset(P):
    """
    Test several functions on a given finite poset.

    The function contains tests of different kinds, for example

    - Numerical properties jump number, dimension etc. can't be a bigger
      in a subposet with one element less.
    - "Dual tests", for example the dual of meet-semilattice must be a join-semilattice.
    - Random tries: for example if the dimension of a poset is `k`, then it can't be the
      intersection of `k-1` random linear extensions.

    EXAMPLES::

        sage: from sage.tests.finite_poset import test_finite_poset
        sage: P = posets.RandomPoset(10, 0.15)
        sage: test_finite_poset(P) is None  # Long time
        True
    """
    from sage.combinat.posets.posets import Poset
    from sage.combinat.subset import Subsets
    from sage.misc.prandom import shuffle

    from sage.misc.call import attrcall

    e = P.random_element()
    P_one_less = P.subposet([x for x in P if x != e])

    # Cardinality
    if len(P) != P.cardinality():
        raise ValueError("error 1 in cardinality")
    if P.cardinality()-1 != P_one_less.cardinality():
        raise ValueError("error 5 in cardinality")

    # Height
    h1 = P.height()
    h2, chain = P.height(certificate=True)
    if h1 != h2:
        raise ValueError("error 1 in height")
    if h1 != len(chain):
        raise ValueError("error 2 in height")
    if not P.is_chain_of_poset(chain):
        raise ValueError("error 3 in height")
    if len(P.random_maximal_chain()) > h1:
        raise ValueError("error 4 in height")
    if h1-P_one_less.height() not in [0, 1]:
        raise ValueError("error 5 in height")

    # Width
    w1 = P.width()
    w2, antichain = P.width(certificate=True)
    if w1 != w2:
        raise ValueError("error 1 in width")
    if w1 != len(antichain):
        raise ValueError("error 2 in width")
    if not P.is_antichain_of_poset(antichain):
        raise ValueError("error 3 in width")
    if len(P.random_maximal_antichain()) > w1:
        raise ValueError("error 4 in width")
    if w1-P_one_less.width() not in [0, 1]:
        raise ValueError("error 5 in width")

    # Dimension
    dim1 = P.dimension()
    dim2, linexts = P.dimension(certificate=True)
    if dim1 != dim2:
        raise ValueError("error 1 in dimension")
    if dim1 != len(linexts):
        raise ValueError("error 2 in dimension")
    P_ = Poset( (P.list(), lambda a, b: all(linext.index(a) < linext.index(b) for linext in linexts)) )
    if P_ != Poset(P.hasse_diagram()):
        raise ValueError("error 3 in dimension")
    x = [P.random_linear_extension() for _ in range(dim1-1)]
    P_ = Poset( (P.list(), lambda a, b: all(linext.index(a) < linext.index(b) for linext in x)) )
    if P_ == Poset(P.hasse_diagram()):
        raise ValueError("error 4 in dimension")
    if dim1-P_one_less.dimension() < 0:
        raise ValueError("error 5 in dimension")

    # Jump number
    j1 = P.jump_number()
    j2, linext = P.jump_number(certificate=True)
    if j1 != j2:
        raise ValueError("error 1 in jump number")
    if P.linear_extension(linext).jump_count() != j1:
        raise ValueError("error 2 in jump number")
    if not P.is_linear_extension(linext):
        raise ValueError("error 3 in jump number")
    if P.linear_extension(P.random_linear_extension()).jump_count() < j1:
        raise ValueError("error 4 in jump number")
    if j1-P_one_less.jump_number() not in [0, 1]:
        raise ValueError("error 5 in jump number")

    P_dual = P.dual()
    selfdual_properties = ['chain', 'bounded', 'connected', 'graded', 'ranked', 'series_parallel', 'slender', 'lattice']
    for prop in selfdual_properties:
        f = attrcall('is_' + prop)
        if f(P) != f(P_dual):
            raise ValueError(f"error in self-dual property {prop}")
    if P.is_graded():
        if P.is_bounded():
            if P.is_eulerian() != P_dual.is_eulerian():
                raise ValueError("error in self-dual property eulerian")
            if P.is_eulerian():
                P_ = P.star_product(P)
                if not P_.is_eulerian():
                    raise ValueError("error in star product / eulerian")
        chain1 = P.random_maximal_chain()
        if len(chain1) != h1:
            raise ValueError("error in is_graded")
        if not P.is_ranked():
            raise ValueError("error in is_ranked / is_graded")

    if P.is_meet_semilattice() != P_dual.is_join_semilattice():
        raise ValueError("error in meet/join semilattice")

    if set(P.minimal_elements()) != set(P_dual.maximal_elements()):
        raise ValueError("error in min/max elements")
    if P.top() != P_dual.bottom():
        raise ValueError("error in top/bottom element")

    parts = P.connected_components()
    P_ = Poset()
    for part in parts:
        P_ = P_.disjoint_union(part)
    if not P.is_isomorphic(P_):
        raise ValueError("error in connected components / disjoint union")
    parts = P.ordinal_summands()
    P_ = Poset()
    for part in parts:
        P_ = P_.ordinal_sum(part)
    if not P.is_isomorphic(P_):
        raise ValueError("error in ordinal summands / ordinal sum")

    P_ = P.with_bounds().without_bounds()
    if not P.is_isomorphic(P_):
        raise ValueError("error in with bounds / without bounds")

    P_ = P.completion_by_cuts().irreducibles_poset()
    if not P.has_isomorphic_subposet(P_):
        raise ValueError("error in completion by cuts / irreducibles poset")

    P_ = P.subposet(Subsets(P).random_element())
    if not P_.is_induced_subposet(P):
        raise ValueError("error in subposet / is induced subposet")

    if not P.is_linear_extension(P.random_linear_extension()):
        raise ValueError("error in is linear extension")

    x = list(P)
    shuffle(x)
    if not P.is_linear_extension(P.sorted(x)):
        raise ValueError("error in sorted")

    dil = P.dilworth_decomposition()
    chain = dil[randint(0, len(dil)-1)]
    if not P.is_chain_of_poset(chain):
        raise ValueError("error in Dilworth decomposition")
    lev = P.level_sets()
    level = lev[randint(0, len(lev)-1)]
    if not P.is_antichain_of_poset(level):
        raise ValueError("error in level sets")

    # certificate=True must return a pair
    bool_with_cert = ['eulerian', 'greedy', 'join_semilattice',
                      'jump_critical', 'meet_semilattice', 'slender']
    for p in bool_with_cert:
        try:  # some properties are not always defined for all posets
            res1 = attrcall('is_' + p)(P)
        except ValueError:
            continue
        res2 = attrcall('is_' + p, certificate=True)(P)
        if not isinstance(res2, tuple) or len(res2) != 2:
            raise ValueError("certificate-option does not return a pair in %s" % p)
        if res1 != res2[0]:
            raise ValueError("certificate-option changes result in %s" % p)
