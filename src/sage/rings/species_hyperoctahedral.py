r"""
Polynomial `r`-species

This module is the analogue of :mod:`sage.rings.species` for polynomial
`r`-species in the sense of Henderson [Henderson2004]_, i.e. functors on
the category `\mathbf B_r` of finite sets with a free action of the
cyclic group `C_r` of order `r`.  The case `r = 2` is the
*hyperoctahedral* case.

The combinatorial model of the canonical permutation representation of
the wreath product `W(r,n) = C_r \wr S_n` lives in
:mod:`sage.groups.perm_gps.hyperoctahedral_group`; the cycle-index
algebra `\Lambda(r)` lives in :mod:`sage.rings.sf_hyperoctahedral`.

A *molecular* `r`-species of degree `n` is a transitive `W(r,n)`-set,
represented by a subgroup `H \leq W(r,n)` up to conjugacy *inside*
`W(r,n)`.  It is important that conjugacy is never taken in the larger
symmetric group `S_{rn}`.

EXAMPLES:

Atomic `r`-species are represented by subgroups of `W(r,n)` up to
`W(r,n)`-conjugacy.  For `r = 2` and `n = 1` there are two of them,
the trivial species `X` and the sign species `X°`::

    sage: from sage.rings.species_hyperoctahedral import (
    ....:     AtomicHyperoctahedralSpecies, MolecularHyperoctahedralSpecies)
    sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
    sage: A = AtomicHyperoctahedralSpecies(2)
    sage: W = _wreath_group(2, 1)
    sage: A(W.subgroup([]))
    X
    sage: A(W)
    X°

For larger `r` there is one degree-one species for every divisor `d`
of `r`, displayed as ``X@d``, where `d` is the order of the stabilizer::

    sage: for r in [1, 2, 4, 6]:
    ....:     A = AtomicHyperoctahedralSpecies(r)
    ....:     print(sorted(A.graded_component(1), key=str))
    [X]
    [X, X°]
    [X, X@2, X°]
    [X, X@2, X@3, X°]

Molecular `r`-species are the products of atomic ones.  They enumerate
the conjugacy classes of subgroups of `W(r,n)`::

    sage: M = MolecularHyperoctahedralSpecies(2)
    sage: W = _wreath_group(2, 2)
    sage: M(W.subgroup([]))
    X^2
    sage: M(W.subgroup([(1, 2), (3, 4)]))
    X°^2
    sage: M(W.subgroup([(1, 2)]))
    X°*X

Conjugacy is taken inside `W(r,n)`: the two sign subgroups of `W(2,2)`
are conjugate, while they would not be conjugate in `S_4`::

    sage: M(W.subgroup([(1, 2)])) == M(W.subgroup([(3, 4)]))
    True

Multisort `r`-species of grade `(n_1,\ldots,n_k)` are transitive sets
for the wreath Young subgroup `W(r;n_1,\ldots,n_k)`.  For example, the
diagonal `C_2` acting simultaneously on a block of each sort is an
atomic species of grade `(1,1)`, because it is directly
indecomposable::

    sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_young_subgroup
    sage: A = AtomicHyperoctahedralSpecies(2, "X, Y")
    sage: W = _wreath_young_subgroup(2, [1, 1])
    sage: a = A(W.subgroup([[(1, 2), (3, 4)]]), {0: [1, 2], 1: [3, 4]})
    sage: a.grade()
    [1, 1]
    sage: a.is_atomic()
    True

By contrast, the full product `C_2 \times C_2` of the sign changes on
each sort decomposes into two molecular factors::

    sage: M = MolecularHyperoctahedralSpecies(2, "X, Y")
    sage: M(W.subgroup([(1, 2), (3, 4)]), {0: [1, 2], 1: [3, 4]})
    X°*Y°
    sage: M(W.subgroup([]), {0: [1, 2], 1: [3, 4]})
    X*Y

REFERENCES:

.. [Henderson2004] Anthony Henderson.
   *Representations of wreath products on cohomology of De Concini-Procesi
   compactifications*.
   International Mathematics Research Notices 2004, no. 20, 981-1021.

AUTHORS:

- Martin Rubey (2026): initial version
"""


# ****************************************************************************
#       Copyright (C) 2026 Martin Rubey <martin.rubey@tuwien.ac.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from itertools import chain, product

from sage.arith.misc import divisors
from sage.categories.graded_algebras_with_basis import GradedAlgebrasWithBasis
from sage.categories.modules_with_basis import Modules
from sage.categories.monoids import Monoids
from sage.categories.tensor import tensor
from sage.categories.sets_cat import cartesian_product
from sage.categories.sets_with_grading import SetsWithGrading
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.integer_vector import IntegerVectors
from sage.combinat.partition import Partition, Partitions
from sage.combinat.partition_tuple import PartitionTuples_level
from sage.combinat.set_partition_ordered import OrderedSetPartitions
from sage.combinat.sf.sf import SymmetricFunctions
from sage.groups.perm_gps.constructor import PermutationGroupElement
from sage.groups.perm_gps.hyperoctahedral_group import (
    _check_grade,
    _wreath_blocks,
    _wreath_group,
    _wreath_young_subgroup,
    _wreath_young_subgroup_on_domain,
    _hyperoctahedral_disjoint_direct_product_decomposition,
)
from sage.groups.perm_gps.permgroup import PermutationGroup, PermutationGroup_generic
from sage.libs.gap.libgap import libgap
from sage.misc.cachefunc import cached_function, cached_method
from sage.misc.fast_methods import WithEqualityById
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.monoids.indexed_free_monoid import (IndexedFreeAbelianMonoid,
                                              IndexedFreeAbelianMonoidElement)
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.sf_hyperoctahedral import HyperoctahedralSymmetricFunctions
from sage.sets.set import Set
from sage.structure.category_object import normalize_names
from sage.structure.element import Element, parent
from sage.structure.parent import Parent
from sage.structure.unique_representation import (UniqueRepresentation,
                                                  WithPicklingByInitArgs)

GAP_FAIL = libgap.eval('fail')

# for each key (currently r, grade, order and orbit sizes) a list of
# canonical representatives of directly indecomposable subgroups of
# the wreath Young subgroup W(r; grade)
_dis_cache = dict()


@cached_function
def _wreath_young_subgroup_classes(r, grade):
    r"""
    Return representatives of the conjugacy classes of subgroups of
    `W(r;` ``grade`` `)`.

    INPUT:

    - ``r`` -- positive integer; the order of the cyclic group
    - ``grade`` -- a tuple of nonnegative integers, the number of
      `C_r`-blocks in each sort

    EXAMPLES::

        sage: from sage.rings.species_hyperoctahedral import _wreath_young_subgroup_classes
        sage: len(_wreath_young_subgroup_classes(2, (2,)))
        8
        sage: all(G.degree() == 4 for G in _wreath_young_subgroup_classes(2, (2,)))
        True
        sage: len(_wreath_young_subgroup_classes(2, (1, 1)))
        5
    """
    grade = _check_grade(grade)
    return _wreath_young_subgroup(r, grade).conjugacy_classes_subgroups()


def _canonical_dompart(r, grade):
    r"""
    Return the canonical partition of the standard domain into sorts.

    The blocks of sort `s` are the consecutive blocks with numbers
    `n_1+\cdots+n_s+1, \ldots, n_1+\cdots+n_{s+1}`.

    INPUT:

    - ``r`` -- positive integer; the order of the cyclic group
    - ``grade`` -- a tuple of nonnegative integers, the number of
      `C_r`-blocks in each sort

    OUTPUT:

    A tuple of frozensets, the `s`-th entry consisting of the points
    of sort `s`.

    EXAMPLES::

        sage: from sage.rings.species_hyperoctahedral import _canonical_dompart
        sage: _canonical_dompart(2, (2, 1))
        (frozenset({1, 2, 3, 4}), frozenset({5, 6}))
        sage: _canonical_dompart(2, (0, 1))
        (frozenset(), frozenset({1, 2}))
    """
    return tuple(frozenset(range(r * sum(grade[:s]) + 1,
                                    r * sum(grade[:s + 1]) + 1))
                 for s in range(len(grade)))


def _standardize_component(H, comp, r):
    r"""
    Return the component ``comp`` of ``H`` as a subgroup of `W(r,k)`.

    The set ``comp`` must be a union of complete consecutive `C_r`-blocks
    of the standard domain of `H`.  Its points are relabelled to
    `\{1, \ldots, rk\}` in order, which is a block-respecting bijection.

    INPUT:

    - ``H`` -- a permutation group on a union of complete `C_r`-blocks
    - ``comp`` -- an iterable; the subset of the domain to restrict to
    - ``r`` -- positive integer; the order of the cyclic group

    OUTPUT:

    A permutation group on `\{1, \ldots, rk\}`.

    EXAMPLES::

        sage: from sage.rings.species_hyperoctahedral import _standardize_component
        sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
        sage: H = _wreath_group(2, 3).subgroup([(3, 4)])
        sage: K = _standardize_component(H, {3, 4}, 2)
        sage: K.degree()
        2
        sage: K == _wreath_group(2, 1)
        True
    """
    comp = sorted(comp)

    gens = []
    for gen in H.gens():
        cycles = [cyc for cyc in gen.cycle_tuples() if cyc[0] in comp]
        if cycles:
            gens.append(PermutationGroupElement(cycles))
    Hc = PermutationGroup(gens, domain=comp)

    relabel = {p: i + 1 for i, p in enumerate(comp)}
    new_gens = []
    for g in Hc.gens():
        perm = list(range(1, len(comp) + 1))
        for p in comp:
            perm[relabel[p] - 1] = relabel[g(p)]
        new_gens.append(PermutationGroupElement(perm))
    return PermutationGroup(new_gens, domain=range(1, len(comp) + 1))


def _wreath_stabilizers(X, a, side, pi, r, arity=1, check=True):
    r"""
    Return the stabilizers of an action of `W(r;n_1,\ldots,n_k)` on the
    set ``X``, together with the assignment of sorts.

    This is the `r`-species analogue of calling
    :func:`~sage.rings.species._stabilizer_subgroups` with a symmetric
    group: the acting group is the wreath Young subgroup
    `W(r;n_1,\ldots,n_k)` built on the domain specified by ``pi``, and
    the returned stabilizers are relabelled so that they act on the
    standard domain `\{1, \ldots, rn\}` with the standard block system.

    INPUT:

    - ``X`` -- the set of structures being acted on
    - ``a`` -- the action, cf.
      :func:`~sage.rings.species._stabilizer_subgroups`
    - ``side`` -- ``'left'`` or ``'right'``
    - ``pi`` -- a dictionary (or iterable) mapping sorts to domains; the
      domains must be unions of complete `C_r`-blocks
    - ``r`` -- positive integer; the order of the cyclic group
    - ``arity`` -- the number of sorts
    - ``check`` -- boolean (default: ``True``); whether to check that
      ``a`` is a group action

    OUTPUT:

    A list of pairs ``(H, pi_H)``, one for each orbit of the action,
    where ``H`` is a subgroup of `W(r,n)` on the standard domain and
    ``pi_H`` assigns the points of its domain to sorts.

    EXAMPLES::

        sage: from sage.rings.species_hyperoctahedral import _wreath_stabilizers
        sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
        sage: W = _wreath_group(2, 2)
        sage: X = list(W.domain())
        sage: a = lambda g, x: g(x)
        sage: H = _wreath_stabilizers(X, a, 'left', {0: list(W.domain())}, 2)
        sage: len(H)
        1
        sage: H[0][0].order()
        2
        sage: H[0][1]
        {0: [1, 2, 3, 4]}

    For several sorts the acting group is the wreath Young subgroup,
    which has two orbits on its own domain::

        sage: H = _wreath_stabilizers(X, a, 'left', {0: [1, 2], 1: [3, 4]}, 2, arity=2)
        sage: len(H)
        2
        sage: H[0][0].order()
        2
        sage: H[0][1]
        {0: [1, 2], 1: [3, 4]}
    """
    from sage.rings.species import _stabilizer_subgroups
    if pi is None:
        raise ValueError("the assignment of sorts to the domain elements must be provided")
    if not isinstance(pi, dict):
        pi = dict(enumerate(pi))
    dompart = [sorted(pi.get(s, [])) for s in range(arity)]
    domain = list(chain.from_iterable(dompart))
    W = _wreath_young_subgroup_on_domain(dompart, r)
    stabilizers = _stabilizer_subgroups(W, X, a, side=side, check=check)

    relabel = {p: i + 1 for i, p in enumerate(sorted(domain))}
    return [(_standardize_component(H, H.domain(), r),
             {s: [relabel[p] for p in b] for s, b in enumerate(dompart)})
            for H in stabilizers]


def _check_standard_domain(G, r):
    r"""
    Check that ``G`` acts on the standard domain `\{1, \ldots, rn\}`.

    INPUT:

    - ``G`` -- a permutation group
    - ``r`` -- positive integer; the order of the cyclic group

    EXAMPLES::

        sage: from sage.rings.species_hyperoctahedral import _check_standard_domain
        sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
        sage: _check_standard_domain(_wreath_group(2, 2), 2)
        sage: _check_standard_domain(_wreath_group(2, 2).subgroup([(1, 2)]), 2)

    A group on a non-standard domain is rejected::

        sage: from sage.groups.perm_gps.permgroup import PermutationGroup
        sage: G = PermutationGroup([(2, 3)], domain=[2, 3, 4, 5])
        sage: _check_standard_domain(G, 2)
        Traceback (most recent call last):
        ...
        ValueError: Permutation Group with generators [(2,3)] must act on the standard domain {1, ..., 4}
    """
    degree = G.degree()
    if degree % r:
        raise ValueError(f"the degree {degree} is not divisible by r = {r}")
    W = _wreath_group(r, degree // r)
    if set(G.domain()) != set(W.domain()):
        raise ValueError(f"{G} must act on the standard domain {{1, ..., {degree}}}")


def _sorted_orbit_list(block):
    r"""
    Return a sorted copy of a list of ``C_r``-orbits.

    INPUT:

    - ``block`` -- an iterable of orbits, each an iterable of labels

    EXAMPLES::

        sage: from sage.rings.species_hyperoctahedral import _sorted_orbit_list
        sage: _sorted_orbit_list([(3, 4), (1, 2)])
        [(1, 2), (3, 4)]
        sage: _sorted_orbit_list([("a", 1), (1, "a")])
        [('a', 1), (1, 'a')]
    """
    try:
        return sorted(block)
    except TypeError:
        return sorted(block, key=str)


def _orbit_label_sets(arity, r, labels):
    r"""
    Return labels as a list of per-sort tuples of orbits.

    INPUT:

    - ``arity`` -- the arity (number of sorts)

    - ``r`` -- positive integer; the order of the cyclic group

    - ``labels`` -- an iterable of ``arity`` iterables of orbits, each
      orbit an iterable of exactly ``r`` labels in cyclic phase order

    The orbits of each sort are returned as a sorted tuple of ``r``-tuples,
    so that the first ``r`` labels of each sort form the first
    `C_r`-orbit, etc.  The labels within an orbit keep their cyclic order.

    EXAMPLES::

        sage: from sage.rings.species_hyperoctahedral import _orbit_label_sets
        sage: _orbit_label_sets(1, 2, [[[1, 2], [3, 4]]])
        [((1, 2), (3, 4))]
        sage: _orbit_label_sets(1, 2, [[[3, 4], [1, 2]]])
        [((1, 2), (3, 4))]

        The order within an orbit is the phase order and is kept::

        sage: _orbit_label_sets(1, 2, [[[2, 1], [3, 4]]])
        [((2, 1), (3, 4))]

        TESTS::

        sage: _orbit_label_sets(1, 2, [[[1, 2]], [[3, 4]]])
        Traceback (most recent call last):
        ...
        ValueError: number of args must match arity of self
        sage: _orbit_label_sets(1, 2, [[[1, 2, 3]]])
        Traceback (most recent call last):
        ...
        ValueError: each orbit of labels must consist of exactly 2 labels, but [(1, 2, 3)] does not
        sage: _orbit_label_sets(1, 2, [[[1, 1]]])
        Traceback (most recent call last):
        ...
        ValueError: the labels of each sort must be distinct, but [(1, 1)] has duplicates
    """
    if len(labels) != arity:
        raise ValueError("number of args must match arity of self")
    result = []
    for orbits in labels:
        orbits = [tuple(orbit) for orbit in orbits]
        if any(len(orbit) != r for orbit in orbits):
            raise ValueError(f"each orbit of labels must consist of exactly "
                             f"{r} labels, but {orbits} does not")
        flat = [x for orbit in orbits for x in orbit]
        if len(set(flat)) != len(flat):
            raise ValueError(f"the labels of each sort must be distinct, "
                             f"but {orbits} has duplicates")
        result.append(tuple(_sorted_orbit_list(orbits)))
    return result


class AtomicHyperoctahedralSpecies(UniqueRepresentation, Parent):
    r"""
    The set of atomic `r`-species.

    An atomic `r`-species of grade `(n_1,\ldots,n_k)` is a directly
    indecomposable transitive `W(r;n_1,\ldots,n_k)`-set, represented up
    to conjugacy inside `W(r;n_1,\ldots,n_k)`.

    INPUT:

    - ``r`` -- positive integer; the order of the cyclic group
    - ``names`` -- an iterable of strings for the sorts of the species
      (default: ``"X"``)

    EXAMPLES::

        sage: from sage.rings.species_hyperoctahedral import AtomicHyperoctahedralSpecies
        sage: A = AtomicHyperoctahedralSpecies(2)
        sage: A
        Atomic 2-species in X
        sage: A.grading_set()
        Integer vectors of length 1

        sage: A = AtomicHyperoctahedralSpecies(2, "X, Y")
        sage: A
        Atomic 2-species in X, Y

    TESTS::

        sage: AtomicHyperoctahedralSpecies(2) is AtomicHyperoctahedralSpecies(ZZ(2))
        True
        sage: AtomicHyperoctahedralSpecies(2) is AtomicHyperoctahedralSpecies(2, "X")
        True
        sage: AtomicHyperoctahedralSpecies(0)
        Traceback (most recent call last):
        ...
        ValueError: r must be a positive integer
    """
    @staticmethod
    def __classcall__(cls, r, names="X"):
        r"""
        Normalize the arguments for unique representation.

        TESTS::

            sage: from sage.rings.species_hyperoctahedral import AtomicHyperoctahedralSpecies
            sage: AtomicHyperoctahedralSpecies(2) is AtomicHyperoctahedralSpecies(2)
            True
        """
        r = ZZ(r)
        if r < 1:
            raise ValueError("r must be a positive integer")
        names = normalize_names(-1, names)
        return super().__classcall__(cls, r, names)

    def __init__(self, r, names):
        r"""
        Initialize the class of atomic `r`-species.

        EXAMPLES::

            sage: from sage.rings.species_hyperoctahedral import AtomicHyperoctahedralSpecies
            sage: A = AtomicHyperoctahedralSpecies(2)
            sage: A._r
            2
            sage: A._arity
            1

        TESTS:

        We have to exclude ``_test_graded_components``, because
        :meth:`~sage.combinat.integer_vector.IntegerVectors.some_elements`
        yields degrees that are too large::

            sage: TestSuite(AtomicHyperoctahedralSpecies(2)).run(skip="_test_graded_components")
            sage: TestSuite(AtomicHyperoctahedralSpecies(2, "X, Y")).run(skip="_test_graded_components")
        """
        category = SetsWithGrading().Infinite()
        Parent.__init__(self, names=names, category=category)
        self._r = ZZ(r)
        self._arity = len(names)
        self._cache = dict()
        # the grades whose standard species have already been renamed
        self._renamed = set()

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.rings.species_hyperoctahedral import AtomicHyperoctahedralSpecies
            sage: AtomicHyperoctahedralSpecies(3)
            Atomic 3-species in X
            sage: AtomicHyperoctahedralSpecies(3, "X, Y")
            Atomic 3-species in X, Y
        """
        return f"Atomic {self._r}-species in {', '.join(self._names)}"

    def _an_element_(self):
        r"""
        Return an element of ``self``.

        This is the diagonal cyclic `r`-species: it acts as a full
        `r`-cycle on one block per sort.

        TESTS::

            sage: from sage.rings.species_hyperoctahedral import AtomicHyperoctahedralSpecies
            sage: AtomicHyperoctahedralSpecies(2).an_element()
            X°
            sage: A = AtomicHyperoctahedralSpecies(2, "X, Y")
            sage: A.an_element()
            {((1,2)(3,4),): ({1, 2}, {3, 4})}

            sage: A = AtomicHyperoctahedralSpecies(1, "X, Y")
            sage: A.an_element()
            {((1,2)(3,4),): ({1, 2}, {3, 4})}
        """
        r = self._r
        k = self._arity
        if r == 1:
            # swap two singletons per sort
            gens = [[(2*s + 1, 2*s + 2) for s in range(k)]]
            pi = {s: [2*s + 1, 2*s + 2] for s in range(k)}
        else:
            # rotate one block per sort
            gens = [[tuple(range(r*s + 1, r*s + r + 1)) for s in range(k)]]
            pi = {s: list(range(r*s + 1, r*s + r + 1)) for s in range(k)}
        G = PermutationGroup(gens)
        return self._element_constructor_(G, pi)

    def _rename(self, grade):
        r"""
        Give the standard atomic `r`-species of the given grade their
        names.

        For a unit grade the atomic `r`-species correspond to the
        divisors `d \mid r`: the stabilizer of the species `C_r / C_d`
        is the cyclic group `C_d` of order `d`.  Following the
        convention of [Henderson2004]_, we display it as ``X@d``, with
        the extremal cases ``X = X@1`` and ``X° = X@r``, and with
        ``X`` for the unique species when `r = 1`.  Every sort receives
        the names obtained from its own name in this way.

        If the grade has a single nonzero entry `n \geq 2` in sort
        `s`, then the compositions `F(X@d)` of the classical species
        `F` of degree `n`, cf.
        :meth:`~sage.rings.species.AtomicSpecies._rename`, with the
        atomic species `X@d` of unit grade in sort `s` are the atomic
        `r`-species `C_d^n \rtimes H`, where `H` is the group of `F`;
        they are named `F_n(X@d)` accordingly.

        Moreover, the type 2 substitutions `X@d(F)` of the atomic
        species `X@d` of unit grade in sort `s` with the classical
        species `F` of degree `n`, cf.
        :meth:`~sage.rings.species_hyperoctahedral.MolecularHyperoctahedralSpecies.Element._type2_substitute_molecular`,
        are atomic `r`-species; they are named `X@d(F_n)` accordingly.
        The case `X(F_n)` is omitted, because `X(F_n) = F_n(X)` has
        already been named above.  In the unlikely event that a type
        2 substitution coincides with an already named species, it is
        left unnamed.

        EXAMPLES::

            sage: from sage.rings.species_hyperoctahedral import AtomicHyperoctahedralSpecies
            sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
            sage: A = AtomicHyperoctahedralSpecies(2)
            sage: A(_wreath_group(2, 1))
            X°
            sage: A = AtomicHyperoctahedralSpecies(4, "X, Y")
            sage: sorted(A.graded_component([0, 1]), key=str)
            [Y, Y@2, Y°]

            sage: A = AtomicHyperoctahedralSpecies(2)
            sage: sorted(A.graded_component(2), key=str)
            [E_2(X), E_2(X°), X°(E_2), {((1,2)(3,4),)}, {((1,4,2,3),)}]
            sage: [a for a in sorted(A.graded_component(3), key=str) if a.get_custom_name()]
            [C_3(X), C_3(X°), E_3(X), E_3(X°), X°(C_3), X°(E_3)]

            sage: A = AtomicHyperoctahedralSpecies(4, "X, Y")
            sage: [a for a in sorted(A.graded_component([0, 2]), key=str) if a.get_custom_name()]
            [E_2(Y), E_2(Y@2), E_2(Y°), Y@2(E_2), Y°(E_2)]
        """
        if sum(grade) == 0:
            return
        nonzero = [(s, n) for s, n in enumerate(grade) if n]
        if len(nonzero) != 1:
            return
        s, n = nonzero[0]

        # the atomic r-species X@d of unit grade in sort s
        W = _wreath_group(self._r, 1)
        cycle = W.gens()[0]
        unit_grade = [0] * self._arity
        unit_grade[s] = 1
        dompart = _canonical_dompart(self._r, unit_grade)
        atoms = {}
        atom_names = {}
        for d in divisors(self._r):
            G = W.subgroup([cycle ** (self._r // d)])
            if d == 1:
                name = self._names[s]
            elif d == self._r:
                name = self._names[s] + "°"
            else:
                name = f"{self._names[s]}@{d}"
            atoms[d] = self(G, dompart, check=False)
            atoms[d].rename(name)
            atom_names[d] = name

        if n == 1:
            return

        # the compositions F(X@d) of the classical species of degree n
        # with X@d, computed with the type 1 substitution
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        from sage.rings.species import MolecularSpecies, _classical_species_groups
        M = MolecularHyperoctahedralSpecies(self._r, self._names)
        molecules = {d: M({atoms[d]: ZZ.one()}) for d in divisors(self._r)}
        # the ordinary molecular species of degree n in sort s and the
        # singleton species of the remaining sorts, which act as
        # placeholders for the type 2 substitution
        O = MolecularSpecies(self._names)
        singletons = [O(SymmetricGroup(1), {i: [1]})
                      for i in range(self._arity)]
        for name, H in _classical_species_groups(n):
            # the classical species as an ordinary molecular species of
            # a single sort for the type 1 substitution and in sort s of
            # all sorts for the type 2 substitution
            F = MolecularSpecies(self._names[s])(H, {0: range(1, n + 1)})
            Fs = O(H, {s: range(1, n + 1)})
            for d in divisors(self._r):
                molecule, = M._type1_substitute_molecular(F, [molecules[d]])._monomial
                molecule.rename(f"{name}({atom_names[d]})")
            # the type 2 substitutions X@d(F), where X(F) = F(X) has
            # been renamed above already
            for d in divisors(self._r):
                if d == 1:
                    continue
                args = list(singletons)
                args[s] = Fs
                molecule, = molecules[d]._type2_substitute_molecular(args)._monomial
                if not molecule.get_custom_name():
                    molecule.rename(f"{atom_names[d]}({name})")

    def _element_constructor_(self, G, pi=None, check=True):
        r"""
        Construct the atomic `r`-species with the given data.

        INPUT:

        - ``G`` -- element of ``self`` (in this case ``pi`` must be
          ``None``) or a permutation group which is a subgroup of
          `W(r,n)` on the standard domain
        - ``pi`` -- a dictionary (or iterable) mapping sorts to
          iterables whose union is the domain of ``G``; if the arity is
          one, ``pi`` can be omitted
        - ``check`` -- boolean (default: ``True``); whether to check
          that ``G`` is directly indecomposable

        EXAMPLES::

            sage: from sage.rings.species_hyperoctahedral import AtomicHyperoctahedralSpecies
            sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
            sage: A = AtomicHyperoctahedralSpecies(2)
            sage: G = _wreath_group(2, 2)
            sage: A(G.subgroup([[(1, 3), (2, 4)]]))
            E_2(X)
            sage: A(G.subgroup([(1, 2), (3, 4)]))
            Traceback (most recent call last):
            ...
            ValueError: ((1,2), (3,4)) is not directly indecomposable

        TESTS::

            The full wreath product `W(2,2)` is `E_2(X^\circ)`::

            sage: A(G)
            E_2(X°)
        """
        if parent(G) is self:
            if pi is None:
                return G
            raise ValueError("cannot reassign sorts to an atomic species")
        if not isinstance(G, PermutationGroup_generic):
            raise ValueError(f"{G} must be a permutation group")
        if check:
            _check_standard_domain(G, self._r)
            if len(_hyperoctahedral_disjoint_direct_product_decomposition(G, self._r)) != 1:
                raise ValueError(f"{G.gens()} is not directly indecomposable")
        if pi is None:
            if self._arity == 1:
                pi = {0: G.domain()}
            else:
                raise ValueError("the assignment of sorts to the domain elements must be provided")
        elif not isinstance(pi, dict):
            pi = dict(enumerate(pi))
        if check:
            if not set(pi).issubset(range(self._arity)):
                raise ValueError(f"keys of pi (={pi.keys()}) must be in range({self._arity})")
            if (sum(len(p) for p in pi.values()) != len(G.domain())
                    or set(chain.from_iterable(pi.values())) != set(G.domain())):
                raise ValueError(f"values of pi (={pi.values()}) must partition the domain of G (={G.domain()})")
        dompart = [sorted(pi.get(s, [])) for s in range(self._arity)]
        return self.element_class(self, G, dompart)

    def grading_set(self):
        r"""
        Return the grading set of ``self``.

        This is the set of integer vectors whose length is the arity of
        ``self``.

        EXAMPLES::

            sage: from sage.rings.species_hyperoctahedral import AtomicHyperoctahedralSpecies
            sage: AtomicHyperoctahedralSpecies(2).grading_set()
            Integer vectors of length 1
            sage: AtomicHyperoctahedralSpecies(2, "X, Y").grading_set()
            Integer vectors of length 2
        """
        return IntegerVectors(length=self._arity)

    def subset(self, size):
        r"""
        Return the set of atomic `r`-species of total degree ``size``.

        EXAMPLES::

            sage: from sage.rings.species_hyperoctahedral import AtomicHyperoctahedralSpecies
            sage: A = AtomicHyperoctahedralSpecies(2)
            sage: sorted(a.degree() for a in A.subset(2))
            [2, 2, 2, 2, 2]

            sage: A = AtomicHyperoctahedralSpecies(2, "X, Y")
            sage: sorted(A.subset(1), key=str)
            [X, X°, Y, Y°]
        """
        result = Set()
        for grade in IntegerVectors(size, length=self._arity):
            result = result.union(self.graded_component(grade))
        return result

    def graded_component(self, grade):
        r"""
        Return the set of atomic `r`-species of grade ``grade``.

        EXAMPLES::

            sage: from sage.rings.species_hyperoctahedral import AtomicHyperoctahedralSpecies
            sage: A = AtomicHyperoctahedralSpecies(2)
            sage: len(A.graded_component(2))
            5

            sage: A = AtomicHyperoctahedralSpecies(2, "X, Y")
            sage: len(A.graded_component([1, 0]))
            2
            sage: len(A.graded_component([1, 1]))
            1

        The unique atom of grade `(1,1)` for `r = 2` is the diagonal
        `C_2`::

            sage: a = next(iter(A.graded_component([1, 1]))); a
            {((1,2)(3,4),): ({1, 2}, {3, 4})}

        For `r = 1` the number of atoms agrees with the number of
        atomic species::

            sage: from sage.rings.species import AtomicSpecies
            sage: A1 = AtomicHyperoctahedralSpecies(1, "X, Y")
            sage: B = AtomicSpecies("X, Y")
            sage: all(len(A1.graded_component(list(mc))) == len(B.graded_component(list(mc)))
            ....:     for mc in IntegerVectors(3, length=2))
            True
        """
        if not hasattr(grade, '__len__'):
            grade = (grade,)
        if len(grade) != self._arity:
            raise ValueError("invalid degree")
        grade = _check_grade(grade)
        return Set([self(G, _canonical_dompart(self._r, grade), check=False)
                    for G in _wreath_young_subgroup_classes(self._r, grade)
                    if len(_hyperoctahedral_disjoint_direct_product_decomposition(G, self._r)) == 1])

    def __contains__(self, x):
        r"""
        Return whether ``x`` is an atomic `r`-species, or a subgroup
        together with an assignment of sorts defining one.

        EXAMPLES::

            sage: from sage.rings.species_hyperoctahedral import AtomicHyperoctahedralSpecies
            sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
            sage: A = AtomicHyperoctahedralSpecies(1)
            sage: _wreath_group(1, 1).subgroup([]) in A
            True

        TESTS::

            sage: A = AtomicHyperoctahedralSpecies(2, "X, Y")
            sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_young_subgroup
            sage: W = _wreath_young_subgroup(2, [1, 1])
            sage: (W.subgroup([[(1, 2), (3, 4)]]), {0: [1, 2], 1: [3, 4]}) in A
            True
            sage: (W.subgroup([(1, 2), (3, 4)]), {0: [1, 2], 1: [3, 4]}) in A
            False
            sage: (W.subgroup([[(1, 2), (3, 4)]]), {0: [1, 3], 1: [2, 4]}) in A
            False
        """
        if parent(x) is self:
            return True
        if isinstance(x, PermutationGroup_generic):
            if self._arity == 1:
                G = x
                pi = {0: G.domain()}
            else:
                return False
        else:
            G, pi = x
            if not isinstance(G, PermutationGroup_generic):
                return False
            if not isinstance(pi, dict):
                pi = dict(enumerate(pi))
        if not set(pi).issubset(range(self._arity)):
            return False
        if (sum(len(p) for p in pi.values()) != len(G.domain())
                or set(chain.from_iterable(pi.values())) != set(G.domain())):
            return False
        degree = G.degree()
        if degree % self._r:
            return False
        W = _wreath_group(self._r, degree // self._r)
        if set(G.domain()) != set(W.domain()):
            return False
        dompart = [sorted(pi.get(s, [])) for s in range(self._arity)]
        # each sort must be a union of complete C_r-blocks
        for s_points in dompart:
            s_set = set(s_points)
            for b in _wreath_blocks(self._r, degree // self._r):
                if len(s_set.intersection(b)) not in (0, len(b)):
                    return False
        for orbit in G.orbits():
            if not any(set(orbit).issubset(s) for s in dompart):
                return False
        if libgap.IsSubgroup(W.gap(), G.gap()) != True:
            return False
        return len(_hyperoctahedral_disjoint_direct_product_decomposition(G, self._r)) == 1

    class Element(WithEqualityById,
                  Element,
                  WithPicklingByInitArgs,
                  metaclass=InheritComparisonClasscallMetaclass):
        r"""
        An atomic `r`-species.
        """
        @staticmethod
        def __classcall__(cls, parent, G, dompart):
            r"""
            Normalize the input for unique representation.

            The `C_r`-blocks of the domain are relabelled so that the
            blocks of sort `0` come first, then the blocks of sort
            `1`, and so on.  Afterwards the group is replaced by the
            canonical representative of its conjugacy class in the
            ambient wreath Young subgroup and the sort partition is the
            canonical one.

            INPUT:

            - ``G`` -- a directly indecomposable permutation group
            - ``dompart`` -- an iterable of `k` iterables, where `k` is
              the arity, assigning each element of the domain of ``G``
              to a sort

            .. WARNING::

                We do not check whether ``G`` is indeed directly
                indecomposable.

            TESTS::

                sage: from sage.rings.species_hyperoctahedral import AtomicHyperoctahedralSpecies
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
                sage: A = AtomicHyperoctahedralSpecies(2)
                sage: G = _wreath_group(2, 2)
                sage: A(G.subgroup([[(1, 3), (2, 4)]])) is A(G.subgroup([[(1, 3), (2, 4)]]))
                True

            The relabelling of the sorts identifies species which
            differ only by the ordering of the blocks.  Here the
            diagonal `C_2` is the unique atom of grade `(1,1)`::

                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_young_subgroup
                sage: A = AtomicHyperoctahedralSpecies(2, "X, Y")
                sage: W = _wreath_young_subgroup(2, [1, 1])
                sage: a = A(W.subgroup([[(1, 2), (3, 4)]]), {0: [1, 2], 1: [3, 4]})
                sage: a is A(W.subgroup([[(1, 2), (3, 4)]]), {0: [3, 4], 1: [1, 2]})
                True
                sage: a is A(W.subgroup([[(1, 2), (3, 4)]]), {0: [2, 1], 1: [4, 3]})
                True

            The subgroup acting on one sort only is not atomic::

                sage: a is A(W.subgroup([(1, 2)]), {0: [1, 2], 1: [3, 4]})
                Traceback (most recent call last):
                ...
                ValueError: ((1,2),) is not directly indecomposable

            A group which is not a subgroup of the wreath Young
            subgroup is rejected.  Here the cyclic group of order
            `4` does not respect the `C_2`-blocks::

                sage: A = AtomicHyperoctahedralSpecies(2)
                sage: A(CyclicPermutationGroup(4))
                Traceback (most recent call last):
                ...
                ValueError: Cyclic group of order 4 as a permutation group is not a subgroup of Permutation Group with generators [(3,4), (1,2), (1,3)(2,4)]
            """
            r = parent._r
            _check_standard_domain(G, r)

            dompart = [sorted(b) for b in dompart]

            # each sort must be a union of complete C_r-blocks
            point_sort = {}
            for s, b in enumerate(dompart):
                point_sort.update({p: s for p in b})
            blocks = _wreath_blocks(r, G.degree() // r)
            for blk in blocks:
                s = point_sort.get(blk[0])
                for p in blk:
                    if point_sort.get(p) != s:
                        raise ValueError(f"the assignment of sorts {dompart} "
                                         f"must be a union of C_{r}-blocks")

            # every orbit of G must be contained in a single sort; the
            # orbits are grouped by sort, in decreasing order of size
            orbits_by_sort = [[] for _ in range(len(dompart))]
            for orbit in G.orbits():
                s = point_sort.get(orbit[0])
                if s is None or not all(point_sort.get(p) == s for p in orbit):
                    raise ValueError(f"all elements of orbit {list(orbit)} "
                                     f"must have the same sort")
                orbits_by_sort[s].append(orbit)
            for orbits in orbits_by_sort:
                orbits.sort(key=lambda o: (-len(o), min(o)))

            grade = tuple(len(b) // r for b in dompart)

            # relabel the blocks so that the blocks of sort 0 come
            # first, then the blocks of sort 1, and so on
            ordered = list(chain.from_iterable(dompart))
            if ordered != list(range(1, G.degree() + 1)):
                relabel = {p: i + 1 for i, p in enumerate(ordered)}
                gens = []
                for gen in G.gens():
                    cycles = [tuple(relabel[p] for p in cyc)
                              for cyc in gen.cycle_tuples()]
                    cycles = [cyc for cyc in cycles if cyc]
                    if cycles:
                        gens.append(PermutationGroupElement(cycles))
                G = PermutationGroup(gens, domain=range(1, G.degree() + 1))
                orbits_by_sort = [[tuple(relabel[p] for p in o)
                                   for o in orbits]
                                  for orbits in orbits_by_sort]

            G_gap = G.gap()
            W = _wreath_young_subgroup(r, grade)
            W_gap = W.gap()
            if libgap.IsSubgroup(W_gap, G_gap) != True:
                raise ValueError(f"{G} is not a subgroup of {W}")

            def new_dis():
                """
                Return a representative of the conjugacy class of
                ``G`` in the wreath Young subgroup.

                Its `C_r`-blocks are ordered in such a way that blocks
                touched by larger orbits have smaller numbers and it
                has a small generating set.
                """
                # the image of each point: the C_r-blocks of each sort
                # are ordered in such a way that blocks touched by
                # larger orbits come first
                pos = {}
                for s, orbits in enumerate(orbits_by_sort):
                    base = sum(grade[:s])
                    placed = set()
                    m = 0
                    for o in orbits:
                        for j in range(base, base + grade[s]):
                            if j not in placed and any(p in o for p in blocks[j]):
                                placed.add(j)
                                for k, p in enumerate(blocks[j]):
                                    pos[p] = (base + m) * r + k + 1
                                m += 1
                gens = []
                for gen in G_gap.SmallGeneratingSet().sage():
                    cycles = [tuple(pos[p] for p in cyc)
                              for cyc in gen.cycle_tuples()]
                    cycles = [cyc for cyc in cycles if len(cyc) > 1]
                    if cycles:
                        gens.append(PermutationGroupElement(cycles))
                return PermutationGroup(gens, domain=range(1, G.degree() + 1))

            # find the canonical representative of the conjugacy class
            # of G in the wreath Young subgroup, or create it
            key = (r, grade, G.order(),
                   tuple(tuple(len(o) for o in orbits)
                         for orbits in orbits_by_sort))
            lookup_dis = _dis_cache.get(key)
            if lookup_dis is None:
                dis = new_dis()
                _dis_cache[key] = [dis]
            else:
                for dis in lookup_dis:
                    if libgap.RepresentativeAction(W_gap, G_gap,
                                                   dis.gap()) != GAP_FAIL:
                        break
                else:
                    dis = new_dis()
                    lookup_dis.append(dis)

            key = (grade, dis)
            if key in parent._cache:
                return parent._cache[key]
            elm = WithPicklingByInitArgs.__classcall__(
                cls, parent, dis, _canonical_dompart(r, grade))
            parent._cache[key] = elm
            return elm

        def __init__(self, parent, dis, domain_partition):
            r"""
            Initialize an atomic `r`-species.

            EXAMPLES::

                sage: from sage.rings.species_hyperoctahedral import AtomicHyperoctahedralSpecies
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
                sage: A = AtomicHyperoctahedralSpecies(2)
                sage: a = A(_wreath_group(2, 2).subgroup([[(1, 3), (2, 4)]]))
                sage: TestSuite(a).run()
                sage: loads(dumps(a)) is a
                True

                sage: A = AtomicHyperoctahedralSpecies(2, "X, Y")
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_young_subgroup
                sage: W = _wreath_young_subgroup(2, [1, 1])
                sage: a = A(W.subgroup([[(1, 2), (3, 4)]]), {0: [1, 2], 1: [3, 4]})
                sage: TestSuite(a).run()
                sage: loads(dumps(a)) is a
                True
            """
            Element.__init__(self, parent)
            self._dis = dis
            self._dompart = domain_partition
            self._mc = tuple(len(v) // parent._r for v in self._dompart)
            self._tc = sum(self._mc)

        def _repr_(self):
            r"""
            Return a string representation of ``self``.

            The first time a grade is displayed, the standard names for
            that grade are installed by :meth:`_rename`.

            EXAMPLES::

                sage: from sage.rings.species_hyperoctahedral import AtomicHyperoctahedralSpecies
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
                sage: A = AtomicHyperoctahedralSpecies(2)
                sage: A(_wreath_group(2, 2).subgroup([[(1, 3), (2, 4)]]))
                E_2(X)

                sage: A = AtomicHyperoctahedralSpecies(2, "X, Y")
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_young_subgroup
                sage: W = _wreath_young_subgroup(2, [1, 1])
                sage: A(W.subgroup([[(1, 2), (3, 4)]]), {0: [1, 2], 1: [3, 4]})
                {((1,2)(3,4),): ({1, 2}, {3, 4})}
            """
            P = self.parent()
            grade = tuple(self._mc)
            if grade not in P._renamed:
                P._renamed.add(grade)
                P._rename(grade)
                return repr(self)
            if P._arity == 1:
                return "{" + f"{self._dis.gens()}" + "}"
            dompart = ', '.join("{" + repr(sorted(b))[1:-1] + "}"
                                for b in self._dompart)
            return "{" + f"{self._dis.gens()}: ({dompart})" + "}"

        def permutation_group(self):
            r"""
            Return the permutation group representing ``self``, together
            with the partition of its domain into sorts.

            The group acts on the standard domain and is the canonical
            representative of its conjugacy class in the ambient wreath
            Young subgroup.  The domain partition consists of complete
            `C_r`-blocks.

            EXAMPLES::

                sage: from sage.rings.species_hyperoctahedral import AtomicHyperoctahedralSpecies
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
                sage: A = AtomicHyperoctahedralSpecies(2)
                sage: A(_wreath_group(2, 1)).permutation_group()
                (Permutation Group with generators [(1,2)], (frozenset({1, 2}),))

                sage: A = AtomicHyperoctahedralSpecies(2, "X, Y")
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_young_subgroup
                sage: W = _wreath_young_subgroup(2, [1, 1])
                sage: A(W.subgroup([[(1, 2), (3, 4)]]), {0: [1, 2], 1: [3, 4]}).permutation_group()
                (Permutation Group with generators [(1,2)(3,4)],
                 (frozenset({1, 2}), frozenset({3, 4})))
            """
            return self._dis, self._dompart

        def degree(self):
            r"""
            Return the degree of ``self``.

            This is the total number of `C_r`-blocks in its domain.

            EXAMPLES::

                sage: from sage.rings.species_hyperoctahedral import AtomicHyperoctahedralSpecies
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
                sage: A = AtomicHyperoctahedralSpecies(2)
                sage: A(_wreath_group(2, 2).subgroup([[(1, 3), (2, 4)]])).degree()
                2
            """
            return self._tc

        def grade(self):
            r"""
            Return the grade of ``self``.

            This is the number of `C_r`-blocks in each sort.

            EXAMPLES::

                sage: from sage.rings.species_hyperoctahedral import AtomicHyperoctahedralSpecies
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
                sage: A = AtomicHyperoctahedralSpecies(2)
                sage: A(_wreath_group(2, 2).subgroup([[(1, 3), (2, 4)]])).grade()
                [2]

                sage: A = AtomicHyperoctahedralSpecies(2, "X, Y")
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_young_subgroup
                sage: W = _wreath_young_subgroup(2, [1, 1])
                sage: A(W.subgroup([[(1, 2), (3, 4)]]), {0: [1, 2], 1: [3, 4]}).grade()
                [1, 1]
            """
            return self.parent().grading_set()(list(self._mc))

        def is_atomic(self):
            r"""
            Return whether ``self`` is atomic, which is always the case.

            EXAMPLES::

                sage: from sage.rings.species_hyperoctahedral import AtomicHyperoctahedralSpecies
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
                sage: A = AtomicHyperoctahedralSpecies(2)
                sage: A(_wreath_group(2, 2).subgroup([[(1, 3), (2, 4)]])).is_atomic()
                True
            """
            return True

        def structures(self, *labels):
            r"""
            Iterate over the structures on the given free `C_r`-set of labels.

            The labels are given as one list of `C_r`-orbits per sort,
            each orbit an iterable of exactly `r` labels in cyclic phase
            order.

            This yields flat tuples of all labels, one representative per
            coset of the stabilizer of ``self`` in the corresponding
            wreath product `W(r; n)`.  The labels are ordered such that
            the first `r` labels form the first `C_r`-orbit of the first
            sort, etc.

            EXAMPLES::

            The cyclic singleton `X^\circ` has a single structure on
            one `C_r`-orbit and none on more; the free singleton has one
            structure per phase::

                sage: from sage.rings.species_hyperoctahedral import AtomicHyperoctahedralSpecies
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
                sage: A = AtomicHyperoctahedralSpecies(2)
                sage: Xo = A(_wreath_group(2, 1))
                sage: Xo
                X°
                sage: list(Xo.structures([('a', 'b')]))
                [('a', 'b')]
                sage: list(Xo.structures([('a', 'b'), ('c', 'd')]))
                []
                sage: X_free = A(_wreath_group(2, 1).subgroup([]))
                sage: sorted(X_free.structures([('a', 'b')]))
                [('a', 'b'), ('b', 'a')]

            The number of structures of the set-like species `E_n(X^\circ)`
            is one on any free `C_r`-set of `n` orbits::

                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_young_subgroup
                sage: E2Xo = A(_wreath_young_subgroup(2, [2]))
                sage: E2Xo
                E_2(X°)
                sage: list(E2Xo.structures([('a', 'b'), ('c', 'd')]))
                [('a', 'b', 'c', 'd')]
                sage: list(A(_wreath_young_subgroup(2, [3])).structures([('a', 'b'), ('c', 'd'), ('e', 'f')]))
                [('a', 'b', 'c', 'd', 'e', 'f')]

            The number of structures is `|W(r; n)| / |L|`::

                sage: a = A(_wreath_group(2, 2).subgroup([[(1, 2), (3, 4)]]))
                sage: a
                {((1,2)(3,4),)}
                sage: sorted(a.structures([('a', 'b'), ('c', 'd')]))
                [('a', 'b', 'c', 'd'), ('a', 'b', 'd', 'c'), ('c', 'd', 'a', 'b'), ('c', 'd', 'b', 'a')]
                sage: a.permutation_group()[0].order()
                2

            Multisort species take one argument per sort, each a list of
            orbits of that sort; the following atom acts as a simultaneous
            flip on one orbit of each sort, so its two structures differ by
            a flip of the labels of the second sort only::

                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_young_subgroup
                sage: AXY = AtomicHyperoctahedralSpecies(2, "X, Y")
                sage: W = _wreath_young_subgroup(2, [1, 1])
                sage: a = AXY(W.subgroup([[(1, 2), (3, 4)]]), {0: [1, 2], 1: [3, 4]})
                sage: a
                {((1,2)(3,4),): ({1, 2}, {3, 4})}
                sage: sorted(a.structures([('a', 'b')], [('c', 'd')]))
                [('a', 'b', 'c', 'd'), ('a', 'b', 'd', 'c')]

            TESTS::

                sage: list(Xo.structures([('a', 'b')], [('c', 'd')]))
                Traceback (most recent call last):
                ...
                ValueError: number of args must match arity of self
                sage: list(Xo.structures([('a', 'b', 'c')]))
                Traceback (most recent call last):
                ...
                ValueError: each orbit of labels must consist of exactly 2 labels, but [('a', 'b', 'c')] does not
                sage: list(Xo.structures([('a', 'a')]))
                Traceback (most recent call last):
                ...
                ValueError: the labels of each sort must be distinct, but [('a', 'a')] has duplicates
            """
            P = self.parent()
            labels = _orbit_label_sets(P._arity, P._r, labels)
            n = tuple(len(s) for s in labels)
            if self._mc != n:
                # wrong number of orbits
                return
            S = _wreath_young_subgroup(P._r, list(n))
            l = [x for s in labels for o in s for x in o]
            for rep in libgap.RightTransversal(S, self._dis):
                yield tuple(S(rep)._act_on_list_on_position(l))


class MolecularHyperoctahedralSpecies(IndexedFreeAbelianMonoid):
    r"""
    The monoid of molecular `r`-species.

    This is the commutative free abelian monoid generated by the atomic
    `r`-species, exactly as
    :class:`~sage.rings.species.MolecularSpecies` is generated by
    :class:`~sage.rings.species.AtomicSpecies`.

    INPUT:

    - ``r`` -- positive integer; the order of the cyclic group
    - ``names`` -- an iterable of strings for the sorts of the species
      (default: ``"X"``)

    EXAMPLES::

        sage: from sage.rings.species_hyperoctahedral import MolecularHyperoctahedralSpecies
        sage: M = MolecularHyperoctahedralSpecies(2)
        sage: M
        Molecular 2-species in X

        sage: M = MolecularHyperoctahedralSpecies(2, "X, Y")
        sage: M
        Molecular 2-species in X, Y

    TESTS::

        sage: from sage.rings.species_hyperoctahedral import MolecularHyperoctahedralSpecies
        sage: MolecularHyperoctahedralSpecies(2) is MolecularHyperoctahedralSpecies(ZZ(2))
        True
        sage: MolecularHyperoctahedralSpecies(2) is MolecularHyperoctahedralSpecies(2, "X")
        True
    """
    @staticmethod
    def __classcall__(cls, r, names="X"):
        r"""
        Normalize the arguments for unique representation.

        TESTS::

            sage: from sage.rings.species_hyperoctahedral import MolecularHyperoctahedralSpecies
            sage: MolecularHyperoctahedralSpecies(2) is MolecularHyperoctahedralSpecies(2)
            True
        """
        r = ZZ(r)
        if r < 1:
            raise ValueError("r must be a positive integer")
        names = normalize_names(-1, names)
        return UniqueRepresentation.__classcall__(cls, r, names)

    def __init__(self, r, names):
        r"""
        Initialize the monoid of molecular `r`-species.

        EXAMPLES::

            sage: from sage.rings.species_hyperoctahedral import MolecularHyperoctahedralSpecies
            sage: M = MolecularHyperoctahedralSpecies(2)
            sage: M._r
            2
            sage: M._arity
            1

        TESTS:

        We have to exclude ``_test_graded_components``, because
        :meth:`~sage.combinat.integer_vector.IntegerVectors.some_elements`
        yields degrees that are too large::

            sage: TestSuite(MolecularHyperoctahedralSpecies(2)).run(skip="_test_graded_components")
            sage: TestSuite(MolecularHyperoctahedralSpecies(2, "X, Y")).run(skip="_test_graded_components")
        """
        indices = AtomicHyperoctahedralSpecies(r, names)
        category = Monoids().Commutative() & SetsWithGrading().Infinite()
        IndexedFreeAbelianMonoid.__init__(self, indices, prefix='',
                                          bracket=False, category=category)
        self._r = ZZ(r)
        self._arity = indices._arity

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.rings.species_hyperoctahedral import MolecularHyperoctahedralSpecies
            sage: MolecularHyperoctahedralSpecies(3)
            Molecular 3-species in X
            sage: MolecularHyperoctahedralSpecies(3, "X, Y")
            Molecular 3-species in X, Y
        """
        return (f"Molecular {self._r}-species in "
                f"{', '.join(self._indices._names)}")

    def _first_ngens(self, n):
        r"""
        Used by the preparser for ``F.<x> = ...``.

        We do not use the generic implementation of
        :class:`sage.monoids.indexed_free_monoid.IndexedFreeAbelianMonoid`,
        because the atomic species cannot be enumerated.

        EXAMPLES::

            sage: from sage.rings.species_hyperoctahedral import MolecularHyperoctahedralSpecies
            sage: M = MolecularHyperoctahedralSpecies(2, "X, Y")
            sage: M._first_ngens(1)
            (X,)
            sage: M._first_ngens(2)
            (X, Y)
        """
        singletons = [sorted(self._indices.graded_component(grade), key=str)[0]
                      for grade in IntegerVectors(1, length=self._arity)]
        return tuple(self.gen(a) for a in singletons[:n])

    def _element_constructor_(self, G, pi=None, check=True):
        r"""
        Construct the molecular `r`-species with the given data.

        INPUT:

        - ``G`` -- one of the following:

          - an element of ``self``
          - a permutation group which is a subgroup of `W(r,n)` on the
            standard domain
          - a dictionary from atoms to exponents
          - a triple ``(X, a, side)`` consisting of a finite set, an
            action and a string ``'left'`` or ``'right'``; the side can
            be omitted, it is then assumed to be ``'right'``

        - ``pi`` -- a dictionary (or iterable) mapping sorts to iterables
          whose union is the domain of ``G``, resp. the domain of the
          acting wreath product if ``G`` is an action; if the arity is
          one and ``G`` is a permutation group, ``pi`` can be omitted
        - ``check`` -- boolean (default: ``True``); whether to check
          the dictionary input and the group action

        EXAMPLES:

        The trivial subgroup of `W(2,2)` is the square of the trivial
        degree-one `r`-species::

            sage: from sage.rings.species_hyperoctahedral import MolecularHyperoctahedralSpecies
            sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
            sage: M = MolecularHyperoctahedralSpecies(2)
            sage: M(_wreath_group(2, 2).subgroup([]))
            X^2

        The subgroup generated by the two sign changes is the product of
        two copies of the sign species::

            sage: M(_wreath_group(2, 2).subgroup([(1, 2), (3, 4)]))
            X°^2

        The subgroup generated by the block interchange is atomic::

            sage: M(_wreath_group(2, 2).subgroup([[(1, 3), (2, 4)]]))
            E_2(X)

        A single sign change acts on one block and fixes the other::

            sage: M(_wreath_group(2, 2).subgroup([(1, 2)]))
            X°*X

        The stabilizer of a point under the natural action of
        `W(2,2)`::

            sage: W = _wreath_group(2, 2)
            sage: X = list(W.domain())
            sage: a = lambda g, x: g(x)
            sage: M((X, a, 'left'), {0: list(W.domain())})
            X*X°

        For several sorts, each `C_r`-block of the domain is assigned to
        a sort.  The trivial subgroup of `W(2;1,1)` decomposes into the
        product of the trivial species in each sort::

            sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_young_subgroup
            sage: M = MolecularHyperoctahedralSpecies(2, "X, Y")
            sage: W = _wreath_young_subgroup(2, [1, 1])
            sage: M(W.subgroup([]), {0: [1, 2], 1: [3, 4]})
            X*Y
            sage: M(W.subgroup([(1, 2), (3, 4)]), {0: [1, 2], 1: [3, 4]})
            X°*Y°

        The diagonal `C_2` is directly indecomposable, whence atomic::

            sage: M(W.subgroup([[(1, 2), (3, 4)]]), {0: [1, 2], 1: [3, 4]})
            {((1,2)(3,4),): ({1, 2}, {3, 4})}

        TESTS::

            sage: TestSuite(M(W.subgroup([(1, 2), (3, 4)]), {0: [1, 2], 1: [3, 4]})).run()
        """
        if parent(G) is self:
            raise ValueError("cannot reassign data to a molecular species")
        if isinstance(G, dict):
            if check:
                if not all(A.parent() is self._indices for A in G):
                    raise ValueError(f"all keys of the dict {G} must be {self._indices}")
                if not all(isinstance(e, Integer) for e in G.values()):
                    raise ValueError(f"all values of the dict {G} must be Integers")
            return self.element_class(self, G)
        if isinstance(G, tuple):
            if len(G) == 2:
                X, a = G
                side = 'right'
            else:
                X, a, side = G
                if side not in ['left', 'right']:
                    raise ValueError(f"the side must be 'right' or 'left', but is {side}")
            stabilizers = _wreath_stabilizers(X, a, side, pi, self._r,
                                              self._arity, check=check)
            if len(stabilizers) > 1:
                raise ValueError("action is not transitive")
            G, pi = stabilizers[0]
        if not isinstance(G, PermutationGroup_generic):
            raise ValueError(f"{G} must be a permutation group")
        _check_standard_domain(G, self._r)
        if pi is None:
            if self._arity == 1:
                pi = {0: G.domain()}
            elif G.degree() == 0:
                pi = {}
            else:
                raise ValueError("the assignment of sorts to the domain elements must be provided")
        elif not isinstance(pi, dict):
            pi = dict(enumerate(pi))
        if not set(pi).issubset(range(self._arity)):
            raise ValueError(f"keys of pi (={pi.keys()}) must be in range({self._arity})")
        if (sum(len(p) for p in pi.values()) != len(G.domain())
                or set(chain.from_iterable(pi.values())) != set(G.domain())):
            raise ValueError(f"values of pi (={pi.values()}) must partition the domain of G (={G.domain()})")
        dompart = [sorted(pi.get(s, [])) for s in range(self._arity)]

        decomposition = _hyperoctahedral_disjoint_direct_product_decomposition(G, self._r)
        result = self.one()
        for comp in decomposition:
            comp = sorted(comp)
            comp_set = set(comp)
            H = _standardize_component(G, comp, self._r)
            relabel = {p: i + 1 for i, p in enumerate(comp)}
            pi_H = {}
            for s in range(self._arity):
                pts = [relabel[p] for p in dompart[s] if p in comp_set]
                if pts:
                    pi_H[s] = pts
            result *= self.gen(self._indices(H, pi_H))
        return result

    def grading_set(self):
        r"""
        Return the grading set of ``self``.

        This is the set of integer vectors whose length is the arity of
        ``self``.

        EXAMPLES::

            sage: from sage.rings.species_hyperoctahedral import MolecularHyperoctahedralSpecies
            sage: MolecularHyperoctahedralSpecies(2).grading_set()
            Integer vectors of length 1
            sage: MolecularHyperoctahedralSpecies(2, "X, Y").grading_set()
            Integer vectors of length 2
        """
        return IntegerVectors(length=self._arity)

    def subset(self, size):
        r"""
        Return the set of molecular `r`-species of total degree ``size``.

        EXAMPLES::

            sage: from sage.rings.species_hyperoctahedral import MolecularHyperoctahedralSpecies
            sage: M = MolecularHyperoctahedralSpecies(2)
            sage: len(M.subset(2))
            8

            sage: M = MolecularHyperoctahedralSpecies(2, "X, Y")
            sage: len(M.subset(2))
            21
        """
        result = Set()
        for grade in IntegerVectors(size, length=self._arity):
            result = result.union(self.graded_component(grade))
        return result

    def graded_component(self, grade):
        r"""
        Return the set of molecular `r`-species of grade ``grade``.

        EXAMPLES::

            sage: from sage.rings.species_hyperoctahedral import MolecularHyperoctahedralSpecies
            sage: M = MolecularHyperoctahedralSpecies(2)
            sage: len(M.graded_component(2))
            8

        The number of molecular `r`-species of degree `n` is the
        number of conjugacy classes of subgroups of `W(r,n)`::

            sage: M = MolecularHyperoctahedralSpecies(2)
            sage: [len(M.graded_component(n)) for n in range(1, 5)]
            [2, 8, 33, 193]

        For several sorts the ambient group is the wreath Young
        subgroup `W(r;n_1,\ldots,n_k)`::

            sage: M = MolecularHyperoctahedralSpecies(2, "X, Y")
            sage: len(M.graded_component([1, 1]))
            5

        Since the molecular species form the free commutative monoid on
        the atomic species, the atomic numbers can be recovered from the
        molecular numbers by the Euler transform::

            sage: from sage.arith.misc import divisors, moebius
            sage: from sage.rings.power_series_ring import PowerSeriesRing
            sage: from sage.rings.species_hyperoctahedral import AtomicHyperoctahedralSpecies
            sage: M = MolecularHyperoctahedralSpecies(2)
            sage: m = [ZZ(1)] + [len(M.graded_component(n)) for n in range(1, 5)]
            sage: R.<x> = PowerSeriesRing(QQ, default_prec=5)
            sage: f = R(m)
            sage: b = [ZZ(0)] + [n * f.log()[n] for n in range(1, 5)]
            sage: a = [sum(moebius(d) * b[n // d] for d in divisors(n)) / n
            ....:      for n in range(1, 5)]
            sage: a
            [2, 5, 19, 120]
            sage: A = AtomicHyperoctahedralSpecies(2)
            sage: [len(A.graded_component(n)) for n in range(1, 5)] == a
            True
        """
        if not hasattr(grade, '__len__'):
            grade = (grade,)
        if len(grade) != self._arity:
            raise ValueError("invalid degree")
        grade = _check_grade(grade)
        return Set([self(rep, _canonical_dompart(self._r, grade))
                    for rep in _wreath_young_subgroup_classes(self._r, grade)])

    def _type1_substitute_molecular(self, M, molecules):
        r"""
        Substitute `r`-molecular species into an ordinary molecular species.

        This is the final, group-theoretic step of the first Henderson
        substitution: the ordinary molecular species ``M`` is an
        abstract coloured species whose sorts correspond to the
        entries of ``molecules``, and every ordinary point is
        replaced by the corresponding substituted `r`-molecule.

        INPUT:

        - ``M`` -- a molecular species of an ordinary
          :class:`~sage.rings.species.PolynomialSpecies`;
          its ``permutation_group()`` supplies the ordinary group
          `H` and the partition of its points into sorts
        - ``molecules`` -- a list of molecular `r`-species, one for
          each sort of ``M``

        OUTPUT:

        An element of ``self``.

        EXAMPLES::

            sage: from sage.rings.species import PolynomialSpecies
            sage: from sage.rings.species_hyperoctahedral import MolecularHyperoctahedralSpecies
            sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
            sage: P = PolynomialSpecies(QQ, "A, B")
            sage: M = MolecularHyperoctahedralSpecies(2)
            sage: X = M(_wreath_group(2, 1).subgroup([]))
            sage: Xo = M(_wreath_group(2, 1))
            sage: G = PermutationGroup([(1, 2)], domain=[1, 2, 3])
            sage: E2AB = P(G, {0: [1, 2], 1: [3]}).support()[0]
            sage: M._type1_substitute_molecular(E2AB, [Xo, X])
            X*E_2(X°)

        The target may be multisort; the sorts of the chunks are
        determined by the substituted molecules::

            sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_young_subgroup
            sage: M = MolecularHyperoctahedralSpecies(2, "X, Y")
            sage: G = _wreath_young_subgroup(2, [1, 1]).subgroup([[(1, 2), (3, 4)]])
            sage: XY = M(G, {0: [1, 2], 1: [3, 4]})
            sage: C3 = PolynomialSpecies(QQ, "X")(CyclicPermutationGroup(3)).support()[0]
            sage: M._type1_substitute_molecular(C3, [XY]).grade()
            [3, 3]

        TESTS:

        The group of `E_2(X^\circ)` is the full wreath product::

            sage: P = PolynomialSpecies(QQ, "X")
            sage: M = MolecularHyperoctahedralSpecies(2)
            sage: Xo = M(_wreath_group(2, 1))
            sage: E2 = P(SymmetricGroup(2)).support()[0]
            sage: G = M._type1_substitute_molecular(E2, [Xo]).permutation_group()[0]
            sage: G.order()
            8
        """
        r = self._r
        H, dompart = M.permutation_group()
        n = H.degree()

        # the sort of each ordinary point
        sort_of = {}
        for i, block in enumerate(dompart):
            for j in block:
                sort_of[j] = i

        # sizes (in C_r-blocks) and offsets of the chunks
        sizes = {}
        offsets = {}
        offset = 0
        for j in range(1, n + 1):
            sizes[j] = molecules[sort_of[j]].degree()
            offsets[j] = offset
            offset += sizes[j]
        N = offset

        gens = []
        # Lift the generators of the ordinary outer group to the chunks.
        for h in H.gens():
            perm = list(range(1, N * r + 1))
            for j in range(1, n + 1):
                jp = h(j)
                for p in range(sizes[j] * r):
                    perm[offsets[j] * r + p] = offsets[jp] * r + p + 1
            gens.append(PermutationGroupElement(perm))

        # Embed the generators of each substituted molecule into
        # every occurrence of its sort.
        for i, block in enumerate(dompart):
            K, _ = molecules[i].permutation_group()
            for j in block:
                shift = offsets[j] * r
                for gen in K.gens():
                    cycles = [tuple(shift + p for p in cyc)
                              for cyc in gen.cycle_tuples()]
                    cycles = [cyc for cyc in cycles if cyc]
                    if cycles:
                        gens.append(PermutationGroupElement(cycles))

        W = PermutationGroup(gens, domain=range(1, N * r + 1))

        # the sort assignment: each chunk is filled with the
        # corresponding molecule, whose points are assigned to sorts by
        # its canonical layout
        pi = {s: [] for s in range(self._arity)}
        for j in range(1, n + 1):
            K_dompart = molecules[sort_of[j]].permutation_group()[1]
            shift = offsets[j] * r
            for s in range(self._arity):
                pi[s].extend(shift + p for p in K_dompart[s])

        return self(W, pi)

    class Element(IndexedFreeAbelianMonoidElement):
        r"""
        A molecular `r`-species.
        """
        @cached_method
        def grade(self):
            r"""
            Return the grade of ``self``.

            This is the vector counting the `C_r`-blocks in each
            sort.

            EXAMPLES::

                sage: from sage.rings.species_hyperoctahedral import MolecularHyperoctahedralSpecies
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
                sage: M = MolecularHyperoctahedralSpecies(2)
                sage: M(_wreath_group(2, 2).subgroup([(1, 2), (3, 4)])).grade()
                [2]

                sage: M = MolecularHyperoctahedralSpecies(2, "X, Y")
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_young_subgroup
                sage: W = _wreath_young_subgroup(2, [1, 1])
                sage: M(W.subgroup([]), {0: [1, 2], 1: [3, 4]}).grade()
                [1, 1]
            """
            P = self.parent()
            S = P.grading_set()
            mons = self._monomial
            if not mons:
                return S([0] * P._arity)
            mc = [sum(n * a._mc[s] for a, n in mons.items())
                  for s in range(P._arity)]
            return S(mc)

        def degree(self):
            r"""
            Return the degree of ``self``.

            This is the total number of `C_r`-blocks in its domain.

            EXAMPLES::

                sage: from sage.rings.species_hyperoctahedral import MolecularHyperoctahedralSpecies
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
                sage: M = MolecularHyperoctahedralSpecies(2)
                sage: M(_wreath_group(2, 2).subgroup([(1, 2), (3, 4)])).degree()
                2
            """
            return Integer(sum(n * a.degree() for a, n in self._monomial.items()))

        def is_molecular(self):
            r"""
            Return whether ``self`` is molecular, which is always the case.

            EXAMPLES::

                sage: from sage.rings.species_hyperoctahedral import MolecularHyperoctahedralSpecies
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
                sage: M = MolecularHyperoctahedralSpecies(2)
                sage: M(_wreath_group(2, 2).subgroup([[(1, 3), (2, 4)]])).is_molecular()
                True
            """
            return True

        def is_atomic(self):
            r"""
            Return whether ``self`` consists of a single atomic species.

            EXAMPLES::

                sage: from sage.rings.species_hyperoctahedral import MolecularHyperoctahedralSpecies
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
                sage: M = MolecularHyperoctahedralSpecies(2)
                sage: M(_wreath_group(2, 2).subgroup([[(1, 3), (2, 4)]])).is_atomic()
                True
                sage: M(_wreath_group(2, 2).subgroup([(1, 2), (3, 4)])).is_atomic()
                False
            """
            return len(self._monomial) == 1 and next(iter(self._monomial.values())) == 1

        @cached_method
        def permutation_group(self):
            r"""
            Return a permutation group representing ``self``, together
            with the partition of its domain into sorts.

            The result acts on the standard domain `\{1, \ldots, rn\}`
            and is the direct product of the permutation groups of the
            atomic factors of ``self``.  The points of each sort of
            each factor are appended after the points of the same sort
            contributed by previous factors.

            EXAMPLES::

                sage: from sage.rings.species_hyperoctahedral import MolecularHyperoctahedralSpecies
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
                sage: M = MolecularHyperoctahedralSpecies(2)
                sage: M(_wreath_group(2, 2).subgroup([(1, 2), (3, 4)])).permutation_group()
                (Permutation Group with generators [(3,4), (1,2)], (frozenset({1, 2, 3, 4}),))
                sage: M.one().permutation_group()
                (Permutation Group with generators [()], (frozenset(),))

            For several sorts the domain partition records which points
            belong to which sort::

                sage: M = MolecularHyperoctahedralSpecies(2, "X, Y")
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_young_subgroup
                sage: W = _wreath_young_subgroup(2, [1, 1])
                sage: G, dompart = M(W.subgroup([(1, 2)]), {0: [1, 2], 1: [3, 4]}).permutation_group()
                sage: G.order()
                2
                sage: dompart
                (frozenset({1, 2}), frozenset({3, 4}))

            An atom whose group acts nontrivially on both sorts, and a
            product of such atoms, are laid out sort by sort::

                sage: G, dompart = M(W.subgroup([[(1, 2), (3, 4)]]), {0: [1, 2], 1: [3, 4]}).permutation_group()
                sage: G, dompart
                (Permutation Group with generators [(1,2)(3,4)], (frozenset({1, 2}), frozenset({3, 4})))
                sage: F = M(W.subgroup([[(1, 2), (3, 4)]]), {0: [1, 2], 1: [3, 4]}) ^ 2
                sage: G, dompart = F.permutation_group()
                sage: sorted(str(gen) for gen in G.gens())
                ['(1,2)(5,6)', '(3,4)(7,8)']
                sage: sorted(map(sorted, dompart))
                [[1, 2, 3, 4], [5, 6, 7, 8]]
            """
            P = self.parent()
            r = P._r
            arity = P._arity

            # the grade of self and the starting point of the points
            # of each sort in the standard domain
            grade = [0] * arity
            for A, e in self._monomial.items():
                for s in range(arity):
                    grade[s] += e * A._mc[s]
            base = [0] * arity
            for s in range(1, arity):
                base[s] = base[s - 1] + r * grade[s - 1]

            offset = [0] * arity
            gens = []
            for A, e in self._monomial.items():
                H = A.permutation_group()[0]
                for _ in range(e):
                    relabel = {}
                    for s in range(arity):
                        points = sorted(A._dompart[s])
                        for i, p in enumerate(points):
                            relabel[p] = base[s] + offset[s] + i + 1
                    for gen in H.gens():
                        cycles = [tuple(relabel[p] for p in cyc)
                                  for cyc in gen.cycle_tuples()]
                        cycles = [cyc for cyc in cycles if cyc]
                        if cycles:
                            gens.append(PermutationGroupElement(cycles))
                    for s in range(arity):
                        offset[s] += len(A._dompart[s])
            dompart = _canonical_dompart(r, tuple(grade))
            return PermutationGroup(gens, domain=range(1, sum(offset) + 1)), dompart

        def cycle_index(self, parent=None):
            r"""
            Return the cycle index of ``self``.

            This is the average over the conjugacy classes of the
            permutation group of ``self`` of the power sums indexed by
            the coloured cycle type, see Henderson [Henderson2004].

            For unisort species, the result is an element of the power
            sum basis of the hyperoctahedral symmetric functions of
            level `r`; for multisort species, a tensor product of such.

            EXAMPLES::

                sage: from sage.rings.species_hyperoctahedral import MolecularHyperoctahedralSpecies
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group, _wreath_young_subgroup
                sage: M = MolecularHyperoctahedralSpecies(2)
                sage: Xo = M(_wreath_group(2, 1))
                sage: Xo.cycle_index()
                1/2*p_1(ζ^1) + 1/2*p_1(ζ^0)

                sage: Xf = M(_wreath_group(2, 1).subgroup([]))
                sage: Xf.cycle_index()
                p_1(ζ^0)

            The cycle index of the atomic species `E_2(X^\circ)` has
            the five terms of Henderson's decomposition of `W(2,2)`::

                sage: E2Xo = M(_wreath_young_subgroup(2, [2]))
                sage: E2Xo.cycle_index()
                1/8*p_{1,1}(ζ^1) + 1/4*p_2(ζ^1) + 1/4*p_1(ζ^0)*p_1(ζ^1) + 1/8*p_{1,1}(ζ^0) + 1/4*p_2(ζ^0)

            The cycle index is multiplicative::

                sage: (Xo * Xo).cycle_index() == Xo.cycle_index()^2
                True
                sage: (Xo * Xf).cycle_index()
                1/2*p_1(ζ^0)*p_1(ζ^1) + 1/2*p_{1,1}(ζ^0)

                sage: M.one().cycle_index()
                1

            For multisort species, the result is a tensor product of
            hyperoctahedral symmetric functions::

                sage: from sage.rings.species_hyperoctahedral import AtomicHyperoctahedralSpecies
                sage: MXY = MolecularHyperoctahedralSpecies(2, "X, Y")
                sage: AXY = AtomicHyperoctahedralSpecies(2, "X, Y")
                sage: d0 = AXY(_wreath_group(2, 1), {0: [1, 2]})
                sage: d1 = AXY(_wreath_group(2, 1), {1: [1, 2]})
                sage: MXY({d0: 1, d1: 1}).cycle_index()
                1/4*p_1(ζ^1) # p_1(ζ^1) + 1/4*p_1(ζ^1) # p_1(ζ^0) + 1/4*p_1(ζ^0) # p_1(ζ^1) + 1/4*p_1(ζ^0) # p_1(ζ^0)

            An atom acting on both sorts::

                sage: W = _wreath_young_subgroup(2, [1, 1])
                sage: d = AXY(W.subgroup([[(1, 2), (3, 4)]]), {0: [1, 2], 1: [3, 4]})
                sage: MXY({d: 1}).cycle_index()
                1/2*p_1(ζ^1) # p_1(ζ^1) + 1/2*p_1(ζ^0) # p_1(ζ^0)

            For `r = 1`, the values agree with the ordinary cycle index
            of species.py, with keys `\lambda` instead of `(\lambda,)`::

                sage: M1 = MolecularHyperoctahedralSpecies(1)
                sage: M1(SymmetricGroup(3), {0: [1, 2, 3]}).cycle_index()
                1/6*p_{1,1,1}(ζ^0) + 1/2*p_{2,1}(ζ^0) + 1/3*p_3(ζ^0)

            TESTS::

            Check that we support different parents::

                sage: from sage.combinat.partition_tuple import PartitionTuples_level
                sage: F = CombinatorialFreeModule(QQ, PartitionTuples_level(2))
                sage: P = Xo.cycle_index(parent=F)
                sage: P
                1/2*B[([], [1])] + 1/2*B[([1], [])]
                sage: P.parent() is F
                True

            This parent should be a module with basis indexed by
            partition tuples::

                sage: Xo.cycle_index(parent=QQ)
                Traceback (most recent call last):
                  ...
                ValueError: `parent` should be a module with basis indexed by partition tuples
            """
            P = self.parent()
            r = P._r
            k = P._arity
            if parent is None:
                p = HyperoctahedralSymmetricFunctions(r,
                                                   SymmetricFunctions(QQ).powersum())
                if k == 1:
                    parent = p
                else:
                    parent = tensor([p] * k)
            elif parent not in Modules.WithBasis:
                raise ValueError("`parent` should be a module with basis indexed "
                                 "by partition tuples")
            base_ring = parent.base_ring()
            Pt = PartitionTuples_level(r)
            G, dompart = self.permutation_group()
            base = [min(s) if s else 1 for s in dompart]

            def cycle_type(g):
                types = []
                for s in range(k):
                    n_s = len(dompart[s]) // r
                    sigma = []
                    signs = []
                    for i in range(n_s):
                        image = g(base[s] + i * r)
                        sigma.append((image - base[s]) // r)
                        signs.append((image - base[s]) % r)
                    parts = [[] for _ in range(r)]
                    seen = [False] * n_s
                    for start in range(n_s):
                        if seen[start]:
                            continue
                        current = start
                        length = 0
                        typ = 0
                        while not seen[current]:
                            seen[current] = True
                            length += 1
                            typ = (typ + signs[current]) % r
                            current = sigma[current]
                        parts[typ].append(length)
                    types.append([Partition(sorted(p, reverse=True))
                                  for p in parts])
                if k == 1:
                    return Pt(types[0])
                return tuple(Pt(t) for t in types)

            return (parent.sum_of_terms([cycle_type(C.an_element()),
                                         base_ring(C.cardinality())]
                                        for C in G.conjugacy_classes())
                    / G.cardinality())

        def structures(self, *labels):
            r"""
            Iterate over the structures on the given free `C_r`-set of labels.

            The labels are given as one list of `C_r`-orbits per sort,
            each orbit an iterable of exactly `r` labels in cyclic phase
            order, see :meth:`AtomicHyperoctahedralSpecies.Element.structures
            <structures>`.

            This yields tuples with one flat structure per atom, in the
            order of the atoms in ``self``; the relabelling is such that
            the first `C_r`-orbits correspond to the first factor in the
            atomic decomposition, etc.

            EXAMPLES::

                sage: from sage.rings.species_hyperoctahedral import MolecularHyperoctahedralSpecies, AtomicHyperoctahedralSpecies
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group, _wreath_young_subgroup
                sage: A = AtomicHyperoctahedralSpecies(2)
                sage: M = MolecularHyperoctahedralSpecies(2)
                sage: Xo = A(_wreath_group(2, 1))
                sage: X_free = A(_wreath_group(2, 1).subgroup([]))

            The product of two free singletons has eight structures::

                sage: sorted(M({X_free: 2}).structures([('a', 'b'), ('c', 'd')]))
                [(('a', 'b'), ('c', 'd')),
                 (('a', 'b'), ('d', 'c')),
                 (('b', 'a'), ('c', 'd')),
                 (('b', 'a'), ('d', 'c')),
                 (('c', 'd'), ('a', 'b')),
                 (('c', 'd'), ('b', 'a')),
                 (('d', 'c'), ('a', 'b')),
                 (('d', 'c'), ('b', 'a'))]

            The number of structures of a molecule is `|W(r; n)| / |L|`,
            where `L` is the group of the molecule; `X^\circ \cdot
            E_2(X^\circ)` has three structures on three orbits::

                sage: E2Xo = A(_wreath_young_subgroup(2, [2]))
                sage: m = M({Xo: 1, E2Xo: 1})
                sage: m
                X°*E_2(X°)
                sage: sorted(m.structures([('a', 'b'), ('c', 'd'), ('e', 'f')]))
                [(('a', 'b'), ('c', 'd', 'e', 'f')),
                 (('c', 'd'), ('a', 'b', 'e', 'f')),
                 (('e', 'f'), ('a', 'b', 'c', 'd'))]
                sage: m.permutation_group()[0].order()
                16

            A molecule has no structures on the wrong number of orbits,
            and a multisort molecule takes one list of orbits per sort::

                sage: list(M({Xo: 2}).structures([('a', 'b'), ('c', 'd'), ('e', 'f')]))
                []
                sage: AXY = AtomicHyperoctahedralSpecies(2, "X, Y")
                sage: MXY = MolecularHyperoctahedralSpecies(2, "X, Y")
                sage: W = _wreath_young_subgroup(2, [1, 1])
                sage: d = AXY(W.subgroup([[(1, 2), (3, 4)]]), {0: [1, 2], 1: [3, 4]})
                sage: sorted(MXY({d: 1}).structures([('a', 'b')], [('c', 'd')]))
                [(('a', 'b', 'c', 'd'),), (('a', 'b', 'd', 'c'),)]
                sage: list(MXY({d: 1}).structures([('a', 'b')]))
                Traceback (most recent call last):
                ...
                ValueError: number of args must match arity of self
            """
            P = self.parent()
            labels = _orbit_label_sets(P._arity, P._r, labels)
            atoms = [a for a, n in self._monomial.items() for _ in range(n)]
            sizes = [a._mc for a in atoms]
            try:
                # raises if the sizes do not match the number of labels
                dissections = [OrderedSetPartitions(list(l), [mc[i] for mc in sizes])
                               for i, l in enumerate(labels)]
            except ValueError:
                return
            for d in product(*dissections):
                # d[i][j] is the set of orbits of sort i given to atom j;
                # sort each block, since the orbit order is positional data
                yield from product(*[
                    a.structures(*[_sorted_orbit_list(d[i][j])
                                   for i in range(P._arity)])
                    for j, a in enumerate(atoms)])

        def _type2_substitute_molecular(self, molecules):
            r"""
            Substitute ordinary molecular species into the sorts of
            ``self``.

            This is the group-theoretic core of the second Henderson
            substitution of [Henderson2004]_, (4.8), in the case where
            every argument is molecular: every `C_r`-orbit of the
            `i`-th sort of ``self`` is replaced by `m_i` new `C_r`-orbits
            corresponding to the points of the molecular structure of
            ``molecules[i]``, transported along the phases.  The result
            is the molecular `r`-species `X_r^{\mathbf N}/H` with the
            sorts of the ``molecules``, in which the sort of a new
            `C_r`-orbit is the sort of the corresponding point of the
            molecule.  The stabilizer `H` of the substituted structure
            on the standard `C_r`-set of `N = \sum_i d_i m_i` orbits is
            generated by:

            - the diagonal lift of the group `L` of ``self``, which maps
              the new orbit `(j, a)` to `(\sigma(j), a)`, where
              `\sigma` is the permutation of the outer orbits induced
              by `L`, shifting the phases of all `a` simultaneously, and

            - the diagonal copies of the groups of the molecules, each
              acting simultaneously on all `r` phases of the new
              orbits of one outer `C_r`-orbit.

            In particular, `|H| = |L| \prod_i |K_i|^{d_i}`, where `K_i`
            is the group of ``molecules[i]``.

            INPUT:

            - ``molecules`` -- a list of ordinary molecular species
              with the same parent, one for each sort of ``self``; each
              must have positive degree

            OUTPUT: an element of the
            :class:`MolecularHyperoctahedralSpecies` with the same `r`
            as ``self`` and the sorts of the ``molecules``.

            EXAMPLES::

                sage: from sage.rings.species import MolecularSpecies
                sage: from sage.rings.species_hyperoctahedral import MolecularHyperoctahedralSpecies
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group, _wreath_young_subgroup
                sage: M = MolecularHyperoctahedralSpecies(2)
                sage: X = M(_wreath_group(2, 1).subgroup([]))
                sage: Xo = M(_wreath_group(2, 1))
                sage: E2 = MolecularSpecies("X")(SymmetricGroup(2))
                sage: X._type2_substitute_molecular([E2])
                E_2(X)

            The two new `C_r`-orbits are exchanged by the diagonal copy
            of `S_2`; the cyclic singleton additionally shifts their
            phases::

                sage: Xo._type2_substitute_molecular([E2])
                X°(E_2)
                sage: _.permutation_group()[0].order()
                4

            Substituting the singleton species leaves every molecular
            `r`-species unchanged, and substituting `E_2` into every
            orbit of `E_2(X^\circ)` gives a transitive species whose
            group has order `8 \cdot 2^2 = 32`::

                sage: E2Xo = M(_wreath_young_subgroup(2, [2]))
                sage: E2Xo
                E_2(X°)
                sage: E1 = MolecularSpecies("X")(SymmetricGroup(1))
                sage: E2Xo._type2_substitute_molecular([E1])
                E_2(X°)
                sage: E2Xo._type2_substitute_molecular([E2]).permutation_group()[0].order()
                32

            The sorts of ``self`` are substituted independently, and the
            sorts of the result are those of the molecules: substituting
            `E_2(X)` into the sort `X` and `E_1(Y)` into the sort `Y` of
            `X^\circ Y^\circ` yields a structure with two orbits of sort
            `X` and one orbit of sort `Y`::

                sage: Mxy = MolecularHyperoctahedralSpecies(2, "X, Y")
                sage: G = _wreath_young_subgroup(2, [1, 1]).subgroup([(1, 2), (3, 4)])
                sage: Q = Mxy(G, {0: [1, 2], 1: [3, 4]})
                sage: Q
                X°*Y°
                sage: Moxy = MolecularSpecies("X, Y")
                sage: E2X = Moxy(SymmetricGroup(2), {0: [1, 2]})
                sage: E1Y = Moxy(SymmetricGroup(1), {1: [1]})
                sage: Q._type2_substitute_molecular([E2X, E1Y])
                X°(E_2)*Y°

            The molecules may be multisort; the sort of a new `C_r`-orbit
            is then the sort of the corresponding point of the molecule.
            Substituting `E_2(Y)` into `X^\circ` gives a structure with
            two orbits, both of sort `Y`::

                sage: Y2 = Moxy(SymmetricGroup(2), {0: [], 1: [1, 2]})
                sage: Xo._type2_substitute_molecular([Y2]).grade()
                [0, 2]
                sage: Xo._type2_substitute_molecular([Y2]).permutation_group()[0].order()
                4

            TESTS:

            At `r=1` the substitution agrees with the ordinary
            composition of molecular species [BLL1998]_::

                sage: E2o = MolecularSpecies("X")(SymmetricGroup(2))
                sage: E2h = MolecularHyperoctahedralSpecies(1)(SymmetricGroup(2))
                sage: G = E2o(E2o).permutation_group()[0]
                sage: E2h._type2_substitute_molecular([E2o]) == MolecularHyperoctahedralSpecies(1)(G)
                True

            Check the number of molecules, their type and their degree.
            Constant terms, that is, molecules of degree zero, are not
            supported, since the substituted structures live on the
            nonempty blocks of a partition::

                sage: Xo._type2_substitute_molecular([E2, E2])
                Traceback (most recent call last):
                ...
                ValueError: the number of molecules must match the arity of self
                sage: Xo._type2_substitute_molecular([X])
                Traceback (most recent call last):
                ...
                ValueError: all molecules must be ordinary molecular species
                sage: Xo._type2_substitute_molecular([MolecularSpecies("X").one()])
                Traceback (most recent call last):
                ...
                ValueError: all molecules must have positive degree
                sage: Moxy = MolecularSpecies("X, Y")
                sage: Y2 = Moxy(SymmetricGroup(2), {0: [], 1: [1, 2]})
                sage: Q._type2_substitute_molecular([E2, Y2])
                Traceback (most recent call last):
                ...
                ValueError: all molecules must have the same parent
            """
            from sage.rings.species import MolecularSpecies

            if len(molecules) != self.parent()._arity:
                raise ValueError("the number of molecules must match "
                                 "the arity of self")
            if not all(isinstance(M, MolecularSpecies.Element)
                       for M in molecules):
                raise ValueError("all molecules must be ordinary "
                                 "molecular species")
            if any(sum(M.grade()) == 0 for M in molecules):
                raise ValueError("all molecules must have positive degree")
            if len(set(M.parent() for M in molecules)) > 1:
                raise ValueError("all molecules must have the same parent")

            r = self.parent()._r
            if not molecules:
                return self.parent().one()
            L, _ = self.permutation_group()
            groups = [M.permutation_group() for M in molecules]
            grade = list(self.grade())
            n = sum(grade)
            arity = molecules[0].parent()._arity

            # the color of each outer C_r-orbit, the degree of the
            # molecule substituted into it, and the offset of its new
            # C_r-orbits
            color = [s for s, d in enumerate(grade) for _ in range(d)]
            m = [sum(molecules[color[j]].grade()) for j in range(n)]
            offsets = [0] * (n + 1)
            for j in range(n):
                offsets[j + 1] = offsets[j] + m[j]
            N = offsets[n]

            gens = []

            # lift the generators of the outer group: the outer orbit
            # j is mapped to the orbit jp, preserving the new orbits
            # and shifting their phases
            for gen in L.gens():
                perm = list(range(1, r * N + 1))
                for j in range(n):
                    # the image of the first point of orbit j determines
                    # the target orbit and the phase shift
                    p = gen(r * j + 1)
                    jp = (p - 1) // r
                    shift = (p - 1) % r
                    for a in range(1, m[j] + 1):
                        for xi in range(r):
                            target = r * (offsets[jp] + a - 1) + (xi + shift) % r + 1
                            perm[r * (offsets[j] + a - 1) + xi] = target
                gens.append(PermutationGroupElement(perm))

            # insert a diagonal copy of the generators of each inner
            # molecule into every outer orbit of its color, acting
            # simultaneously on all r phases of the new orbits
            for j in range(n):
                K = groups[color[j]][0]
                base = offsets[j]
                for gen in K.gens():
                    cycles = []
                    for cyc in gen.cycle_tuples():
                        for xi in range(r):
                            cycles.append(tuple(r * (base + a - 1) + xi + 1
                                               for a in cyc))
                    cycles = [cyc for cyc in cycles if len(cyc) > 1]
                    if cycles:
                        gens.append(PermutationGroupElement(cycles))

            H = PermutationGroup(gens, domain=range(1, r * N + 1))

            # the sort of the new C_r-orbit (j, a) is the sort of the
            # point a in the molecular structure substituted into the
            # outer orbit j, so the orbits of one molecule are assigned
            # to the sorts of its points
            pi = {t: [] for t in range(arity)}
            for j in range(n):
                K_dompart = groups[color[j]][1]
                for t, v in enumerate(K_dompart):
                    for a in v:
                        pi[t].extend(range(r * (offsets[j] + a - 1) + 1,
                                           r * (offsets[j] + a - 1) + r + 1))
            return MolecularHyperoctahedralSpecies(
                r, molecules[0].parent()._indices._names)(H, pi)


class PolynomialHyperoctahedralSpecies(CombinatorialFreeModule):
    r"""
    The ring of polynomial `r`-species.

    This is the commutative graded algebra over ``base_ring`` whose
    basis is given by the molecular `r`-species, exactly as
    :class:`~sage.rings.species.PolynomialSpecies` is the graded algebra
    whose basis is given by the molecular species.

    INPUT:

    - ``base_ring`` -- a ring
    - ``r`` -- positive integer; the order of the cyclic group
    - ``names`` -- an iterable of strings for the sorts of the
      species (default: ``"X"``)

    EXAMPLES::

        sage: from sage.rings.species_hyperoctahedral import PolynomialHyperoctahedralSpecies
        sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
        sage: P = PolynomialHyperoctahedralSpecies(QQ, 2)
        sage: P
        Polynomial 2-species in X over Rational Field
        sage: W = _wreath_group(2, 1)
        sage: P(W.subgroup([])) * P(W)
        X*X°

    The product is associative, commutative and distributive::

        sage: P = PolynomialHyperoctahedralSpecies(QQ, 2)
        sage: W1 = _wreath_group(2, 1)
        sage: W2 = _wreath_group(2, 2)
        sage: a = P(W1)
        sage: b = P(W2.subgroup([]))
        sage: c = P(W2.subgroup([[(1, 3), (2, 4)]]))
        sage: a * b == b * a
        True
        sage: (a * b) * c == a * (b * c)
        True
        sage: a * (b + c) == a * b + a * c
        True
        sage: (a * b).degree()
        3

    For several sorts the generators are the degree one species in
    each sort::

        sage: P.<X, Y> = PolynomialHyperoctahedralSpecies(QQ, 2)
        sage: W = _wreath_group(2, 1)
        sage: P(W.subgroup([]), {0: [1, 2]}) * P(W, {1: [1, 2]})
        X*Y°

    TESTS::

        sage: from sage.rings.species_hyperoctahedral import PolynomialHyperoctahedralSpecies
        sage: PolynomialHyperoctahedralSpecies(QQ, 2) is PolynomialHyperoctahedralSpecies(QQ, ZZ(2))
        True
        sage: PolynomialHyperoctahedralSpecies(QQ, 2) is PolynomialHyperoctahedralSpecies(QQ, ZZ(2), "X")
        True
        sage: PolynomialHyperoctahedralSpecies(QQ, 0)
        Traceback (most recent call last):
        ...
        ValueError: r must be a positive integer
    """
    @staticmethod
    def __classcall__(cls, base_ring, r, names="X"):
        r"""
        Normalize the arguments for unique representation.

        TESTS::

            sage: from sage.rings.species_hyperoctahedral import PolynomialHyperoctahedralSpecies
            sage: PolynomialHyperoctahedralSpecies(QQ, 2) is PolynomialHyperoctahedralSpecies(QQ, 2)
            True
        """
        r = ZZ(r)
        if r < 1:
            raise ValueError("r must be a positive integer")
        names = normalize_names(-1, names)
        return super().__classcall__(cls, base_ring, r, names)

    def __init__(self, base_ring, r, names):
        r"""
        Initialize the ring of polynomial `r`-species.

        TESTS::

            sage: from sage.rings.species_hyperoctahedral import PolynomialHyperoctahedralSpecies
            sage: P = PolynomialHyperoctahedralSpecies(ZZ, 2)
            sage: TestSuite(P).run()
            sage: P2 = PolynomialHyperoctahedralSpecies(ZZ, 2, "X, Y")
            sage: TestSuite(P2).run()
        """
        self._r = ZZ(r)
        self._arity = len(names)
        category = GradedAlgebrasWithBasis(base_ring).Commutative()
        CombinatorialFreeModule.__init__(self, base_ring,
                                         basis_keys=MolecularHyperoctahedralSpecies(r, names),
                                         category=category,
                                         prefix='', bracket=False)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.rings.species_hyperoctahedral import PolynomialHyperoctahedralSpecies
            sage: PolynomialHyperoctahedralSpecies(ZZ, 3)
            Polynomial 3-species in X over Integer Ring
            sage: PolynomialHyperoctahedralSpecies(ZZ, 3, "X, Y")
            Polynomial 3-species in X, Y over Integer Ring
        """
        names = self._indices._indices._names
        return (f"Polynomial {self._r}-species in {', '.join(names)} "
                f"over {self.base_ring()}")

    def _first_ngens(self, n):
        r"""
        Used by the preparser for ``F.<x> = ...``.

        We do not use the generic implementation of
        :class:`sage.combinat.CombinatorialFreeModule`, because we do
        not want to implement `gens`.

        EXAMPLES::

            sage: from sage.rings.species_hyperoctahedral import PolynomialHyperoctahedralSpecies
            sage: P.<X, Y> = PolynomialHyperoctahedralSpecies(QQ, 2)  # indirect doctest
            sage: X.degree()
            1
            sage: X + 2*Y
            X + 2*Y

        Only the first ``n`` degree one species are returned::

            sage: P = PolynomialHyperoctahedralSpecies(QQ, 2, "X, Y")
            sage: P._first_ngens(1)
            (X,)
            sage: P._first_ngens(2)
            (X, Y)
        """
        singletons = [sorted(self._indices._indices.graded_component(grade), key=str)[0]
                      for grade in IntegerVectors(1, length=self._arity)]
        return tuple(self(a) for a in singletons[:n])

    def _element_constructor_(self, G, pi=None, check=True):
        r"""
        Construct the polynomial `r`-species with the given data.

        INPUT:

        - ``G`` -- one of the following:

          - an element of ``self``
          - an atomic or molecular `r`-species
          - a permutation group which is a subgroup of `W(r,n)` on the
            standard domain
          - a dictionary from molecular `r`-species to elements of the
            base ring
          - a triple ``(X, a, side)`` consisting of a finite set, an
            action and a string ``'left'`` or ``'right'``; the side can
            be omitted, it is then assumed to be ``'right'``

        - ``pi`` -- a dictionary (or iterable) mapping sorts to
          iterables whose union is the domain of ``G``, resp. the
          domain of the acting wreath product if ``G`` is an action;
          if the arity is one and ``G`` is a permutation group, ``pi``
          can be omitted
        - ``check`` -- boolean (default: ``True``); whether to check
          the dictionary input and the group action

        EXAMPLES::

            sage: from sage.rings.species_hyperoctahedral import (
            ....:     PolynomialHyperoctahedralSpecies, MolecularHyperoctahedralSpecies)
            sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
            sage: P = PolynomialHyperoctahedralSpecies(QQ, 2)
            sage: M = MolecularHyperoctahedralSpecies(2)
            sage: W = _wreath_group(2, 2)
            sage: P(M(W.subgroup([])))
            X^2
            sage: P(W.subgroup([(1, 2), (3, 4)]))
            X°^2

        A subgroup which is not directly indecomposable is decomposed
        into molecular factors::

            sage: P(W.subgroup([(1, 2)]))
            X°*X

        The stabilizer of a point under the natural action of
        `W(2,2)`::

            sage: X = list(W.domain())
            sage: a = lambda g, x: g(x)
            sage: P((X, a, 'left'), {0: list(W.domain())})
            X*X°

        For several sorts the assignment of the `C_r`-blocks to sorts
        must be provided::

            sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_young_subgroup
            sage: P = PolynomialHyperoctahedralSpecies(QQ, 2, "X, Y")
            sage: W = _wreath_young_subgroup(2, [1, 1])
            sage: P(W.subgroup([[(1, 2), (3, 4)]]), {0: [1, 2], 1: [3, 4]})
            {((1,2)(3,4),): ({1, 2}, {3, 4})}
            sage: P(W.subgroup([]), {0: [1, 2], 1: [3, 4]})
            X*Y

        TESTS::

            sage: P = PolynomialHyperoctahedralSpecies(QQ, 2)
            sage: f = P.one()
            sage: P(f) is f
            True
        """
        if parent(G) is self:
            raise ValueError("cannot reassign data to a polynomial species")

        if isinstance(G, AtomicHyperoctahedralSpecies.Element):
            G = self._indices({G: ZZ.one()})

        if isinstance(G, MolecularHyperoctahedralSpecies.Element):
            if check and G.parent() is not self._indices:
                raise ValueError(f"{G} must be a {self._indices}")
            return self._from_dict({G: self.base_ring().one()})

        if isinstance(G, dict):
            if check:
                if not all(M.parent() is self._indices for M in G):
                    raise ValueError(f"all keys of the dict {G} must be {self._indices}")
                if not all(e in self.base_ring() for e in G.values()):
                    raise ValueError(f"all values of the dict {G} must be in {self.base_ring()}")
            return self._from_dict(G)

        if isinstance(G, tuple):
            if len(G) == 2:
                X, a = G
                side = 'right'
            else:
                X, a, side = G
                if side not in ['left', 'right']:
                    raise ValueError(f"the side must be 'right' or 'left', but is {side}")
            stabilizers = _wreath_stabilizers(X, a, side, pi, self._r,
                                              self._arity, check=check)
            result = self.zero()
            for H, pi_H in stabilizers:
                result += self.monomial(self._indices(H, pi_H))
            return result

        if isinstance(G, PermutationGroup_generic):
            M = self._indices(G, pi)
            return self._from_dict({M: ZZ.one()})

        raise ValueError(f"{G} must be an element of the base ring, a permutation "
                         "group or a molecular species")

    def change_ring(self, R):
        r"""
        Return the base change of ``self`` to `R`.

        EXAMPLES::

            sage: from sage.rings.species_hyperoctahedral import PolynomialHyperoctahedralSpecies
            sage: P = PolynomialHyperoctahedralSpecies(ZZ, 2)
            sage: P.change_ring(QQ)
            Polynomial 2-species in X over Rational Field
            sage: P.change_ring(ZZ) is P
            True
            sage: P2 = PolynomialHyperoctahedralSpecies(ZZ, 2, "X, Y")
            sage: P2.change_ring(QQ)
            Polynomial 2-species in X, Y over Rational Field
        """
        if R is self.base_ring():
            return self
        return PolynomialHyperoctahedralSpecies(R, self._r,
                                                self._indices._indices._names)

    def degree_on_basis(self, m):
        r"""
        Return the degree of the molecular `r`-species ``m``.

        EXAMPLES::

            sage: from sage.rings.species_hyperoctahedral import PolynomialHyperoctahedralSpecies
            sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
            sage: P = PolynomialHyperoctahedralSpecies(QQ, 2)
            sage: P.degree_on_basis(P(_wreath_group(2, 3)).support()[0])
            3
        """
        return m.degree()

    @cached_method
    def one_basis(self):
        r"""
        Return the index of the multiplicative unit.

        EXAMPLES::

            sage: from sage.rings.species_hyperoctahedral import PolynomialHyperoctahedralSpecies
            sage: PolynomialHyperoctahedralSpecies(QQ, 2).one_basis()
            1
        """
        return self._indices.one()

    def product_on_basis(self, H, K):
        r"""
        Return the product of the basis elements indexed by ``H`` and ``K``.

        The product of two molecular `r`-species is their direct product
        with disjoint supports, which is the product in the free
        commutative monoid of molecular `r`-species.

        EXAMPLES::

            sage: from sage.rings.species_hyperoctahedral import PolynomialHyperoctahedralSpecies
            sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
            sage: P = PolynomialHyperoctahedralSpecies(QQ, 2)
            sage: M = P._indices
            sage: H = M(_wreath_group(2, 1).subgroup([]))
            sage: K = M(_wreath_group(2, 1))
            sage: P.product_on_basis(H, K)
            X*X°
        """
        return self.element_class(self, {H * K: ZZ.one()})

    def _exponential(self, multiplicities, degrees):
        r"""
        Return `E_r(\sum_i m_i X_i)` in the specified degrees.

        The weighted `r`-exponential is expanded into molecular species
        using the `\lambda`-ring coefficients `m_\lambda(c)`, which are
        obtained from the monomial symmetric functions by substituting
        the Adams operations `p_k \mapsto \psi_k(c)`:

        .. MATH::

            E_r(cX)
            = \sum_{\lambda} m_\lambda(c)
              \frac{X_r^{|\lambda|}}{W(r;\lambda)},

        where `W(r;\lambda)` is the wreath Young subgroup
        `W(r,\lambda_1)\times W(r,\lambda_2)\times\cdots`, so that
        `X_r^d/W(r;\lambda)` is the molecular `r`-species
        `E(r)_{\lambda_1}E(r)_{\lambda_2}\cdots`, the inflation of the
        ordinary molecular species `E_\lambda`.

        INPUT:

        - ``multiplicities`` -- a list of weights, elements of the base
          ring, of length the arity of ``self``

        - ``degrees`` -- a list of nonnegative integers, of the same
          length, such that the result is homogeneous of degree
          ``degrees[i]`` in sort ``i``

        EXAMPLES::

            sage: from sage.rings.species_hyperoctahedral import PolynomialHyperoctahedralSpecies
            sage: P = PolynomialHyperoctahedralSpecies(QQ, 2)
            sage: P._exponential([1], [3])
            E_3(X°)
            sage: P._exponential([3/2], [3])
            3/2*E_3(X°) + 3/4*E_2(X°)*X° - 1/16*X°^3

        The weights use the usual lambda-ring structure of the base
        ring.  In particular, if `q` has lambda degree one, then
        `\psi_n(q) = q^n`::

            sage: R.<q> = QQ[]
            sage: P = PolynomialHyperoctahedralSpecies(R, 2)
            sage: P._exponential([1+q], [2])
            (q^2+1)*E_2(X°) + q*X°^2
            sage: P._exponential([1+q], [3])
            (q^3+1)*E_3(X°) + (q^2+q)*E_2(X°)*X°
            sage: P._exponential([1-q], [2])
            (-q^2+1)*E_2(X°) + (q^2-q)*X°^2

        Each sort contributes a separate factor::

            sage: P = PolynomialHyperoctahedralSpecies(R, 2, "X, Y")
            sage: P._exponential([1, q], [2, 2])
            q^2*E_2(X°)*E_2(Y°)

        TESTS::

            sage: P = PolynomialHyperoctahedralSpecies(QQ, 2)
            sage: P._exponential([1], [0])
            1
            sage: P._exponential([1], [0]).parent() is P
            True
        """
        B = self.base_ring()
        Sym = SymmetricFunctions(B)
        p = Sym.p()
        m = Sym.m()

        def stretch(c, k):
            r"""
            Substitute in ``c`` all variables appearing in the
            base ring with their ``k``-th power.
            """
            if callable(c):
                return c(*[g ** k for g in B.gens() if g != B.one()])
            return c

        def monomial(c, lam):
            r"""
            Return `m_\lambda(c)`, the monomial symmetric function
            evaluated at the virtual alphabet `c`, by expanding
            `m_\lambda` into powersum symmetric functions and
            substituting `p_k \mapsto \psi_k(c)`.
            """
            total = B.zero()
            for mu, coeff in p(m[lam]).monomial_coefficients().items():
                prod = B.one()
                for k, mult in mu.to_exp_dict().items():
                    prod *= stretch(c, k) ** mult
                total += coeff * prod
            return total

        r = self._r

        def factor(s, c, d):
            r"""
            Return `E_r(c X_s)_d` as a molecular `r`-species in
            sort ``s``.
            """
            grade = [0] * self._arity
            grade[s] = d
            dompart = _canonical_dompart(r, grade)
            return self.sum(monomial(c, lam)
                            * self(_wreath_young_subgroup(r, list(lam)), dompart)
                            for lam in Partitions(d))

        return self.prod(factor(s, multiplicities[s], degrees[s])
                         for s in range(self._arity))

    class Element(CombinatorialFreeModule.Element):
        r"""
        A (virtual) polynomial `r`-species.
        """
        def is_constant(self):
            r"""
            Return whether ``self`` is a constant polynomial `r`-species.

            EXAMPLES::

                sage: from sage.rings.species_hyperoctahedral import PolynomialHyperoctahedralSpecies
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
                sage: P = PolynomialHyperoctahedralSpecies(QQ, 2)
                sage: X = P(_wreath_group(2, 1).subgroup([]))
                sage: X.is_constant()
                False
                sage: (3*P.one()).is_constant()
                True
                sage: P(0).is_constant()
                True
                sage: (1 + X).is_constant()
                False
            """
            return self.is_zero() or not self.maximal_degree()

        def is_virtual(self):
            r"""
            Return whether ``self`` is a virtual polynomial `r`-species.

            EXAMPLES::

                sage: from sage.rings.species_hyperoctahedral import PolynomialHyperoctahedralSpecies
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
                sage: P = PolynomialHyperoctahedralSpecies(QQ, 2)
                sage: X = P(_wreath_group(2, 1).subgroup([]))
                sage: Y = P(_wreath_group(2, 1))
                sage: V = 2*X - 3*Y
                sage: V.is_virtual()
                True
                sage: (X*Y).is_virtual()
                False
            """
            return any(c < 0 for c in self.coefficients(sort=False))

        def is_molecular(self):
            r"""
            Return whether ``self`` is a molecular `r`-species.

            EXAMPLES::

                sage: from sage.rings.species_hyperoctahedral import PolynomialHyperoctahedralSpecies
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
                sage: P = PolynomialHyperoctahedralSpecies(QQ, 2)
                sage: X = P(_wreath_group(2, 1).subgroup([]))
                sage: Y = P(_wreath_group(2, 1))
                sage: (2*X).is_molecular()
                False
                sage: (X*Y).is_molecular()
                True
            """
            coefficients = self.coefficients(sort=False)
            return len(coefficients) == 1 and coefficients[0] == 1

        def is_atomic(self):
            r"""
            Return whether ``self`` is an atomic `r`-species.

            EXAMPLES::

                sage: from sage.rings.species_hyperoctahedral import PolynomialHyperoctahedralSpecies
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
                sage: P = PolynomialHyperoctahedralSpecies(QQ, 2)
                sage: X = P(_wreath_group(2, 1).subgroup([]))
                sage: Y = P(_wreath_group(2, 1))
                sage: (2*Y).is_atomic()
                False
                sage: (X*Y).is_atomic()
                False
                sage: Y.is_atomic()
                True
            """
            return self.is_molecular() and self.support()[0].is_atomic()

        def _compose_with_singletons(self, names, args):
            r"""
            Return the type 2 substitution of ``self`` with sums of
            singleton species.

            This is the analogue for `r`-species of Exercise 2.6.16 in
            [BLL1998]_, using the second kind of substitution of
            [Henderson2004]_, (4.8).  The homogeneous `k`-sort
            `r`-species ``self`` is substituted with the sums of
            singleton species `X_{1,1}+\cdots+X_{1,m_1},\ldots,
            X_{k,1}+\cdots+X_{k,m_k}` and the result is expanded into
            multisort molecular `r`-species.

            For a molecular term `X_1^{d_1}\cdots X_k^{d_k}/L` of
            ``self`` and compositions `e^{(i)} = (e_{i,1},\ldots,e_{i,m_i})`
            of the `d_i`, let `E = (e_{i,j})`, let `W(r;E)` be the
            corresponding wreath Young subgroup and let
            `\mathbf d = (d_1,\ldots,d_k)`.  The inner sum of the
            expansion is

            .. MATH::

                \sum_{\tau\in W(r;E)\backslash W(r;\mathbf d)/L}
                \frac{\prod_{i,j} X_{i,j}^{e_{i,j}}}
                     {W(r;E)\cap\tau L\tau^{-1}},

            exactly as in the ordinary case: the isomorphism classes of
            the substituted structures on the standard graded `C_r`-set
            of grade `E` are indexed by the double cosets, the stabilizer
            of the structure corresponding to `\tau` being
            `W(r;E)\cap\tau L\tau^{-1}`.

            INPUT:

            - ``names`` -- the (flat) list of the names of the sorts of
              the result, whose length is the total number of parts of
              the compositions in ``args``

            - ``args`` -- a sequence of `k` compositions, where `k` is
              the arity of ``self``.  The parts of the `i`-th
              composition sum to the grade of ``self`` in sort `i` and
              count `C_r`-orbits.

            OUTPUT: the polynomial `r`-species

            ``self(X_{1,1}+\cdots+X_{1,m_1},...,X_{k,1}+\cdots+X_{k,m_k})``

            in the sorts given by ``names``.

            EXAMPLES:

            Substituting `X+Y` into `E_2(X)`, whose group interchanges
            the two `C_2`-orbits, yields the free species `X\cdot Y`::

                sage: from sage.rings.species_hyperoctahedral import PolynomialHyperoctahedralSpecies
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
                sage: P = PolynomialHyperoctahedralSpecies(QQ, 2)
                sage: W = _wreath_group(2, 2)
                sage: E2X = P(W.subgroup([[(1, 3), (2, 4)]]))
                sage: E2X._compose_with_singletons(["X", "Y"], [[1, 1]])
                X*Y

            Substituting `X+Y` into `X^{\circ 2}` yields two copies of
            `X^\circ Y^\circ`: either of the two `C_2`-orbits may
            carry the `X`-structure::

                sage: P(W.subgroup([(1, 2), (3, 4)]))._compose_with_singletons(["X", "Y"], [[1, 1]])
                2*X°*Y°

            Compositions with empty parts are allowed::

                sage: E2X._compose_with_singletons(["X", "Y"], [[2, 0]])
                E_2(X)

            The molecular terms of a polynomial species are expanded
            separately::

                sage: (E2X + P(W.subgroup([(1, 2), (3, 4)])))._compose_with_singletons(["X", "Y"], [[1, 1]])
                X*Y + 2*X°*Y°

            Each sort of a multisort species is refined separately.
            Here the sort `X` of `E_2(X)\cdot Y^\circ` is split into
            `X_1` and `X_2`::

                sage: P2 = PolynomialHyperoctahedralSpecies(QQ, 2, "X, Y")
                sage: E2XYo = P2(W.subgroup([[(1, 3), (2, 4)]]), {0: [1, 2, 3, 4], 1: []})
                sage: E2XYo *= P2(_wreath_group(2, 1), {0: [], 1: [1, 2]})
                sage: E2XYo._compose_with_singletons(["X1", "X2", "Y"], [[1, 1], [1]])
                X1*X2*Y°

            For `r = 1` the result agrees with the ordinary
            substitution of Exercise 2.6.16 in [BLL1998]_::

                sage: P1 = PolynomialHyperoctahedralSpecies(QQ, 1)
                sage: C4 = P1(CyclicPermutationGroup(4))
                sage: F1 = C4._compose_with_singletons(["X", "Y"], [[2, 2]]); F1
                X^2*Y^2 + {((1,2)(3,4),): ({1, 2}, {3, 4})}

                sage: from sage.rings.species import PolynomialSpecies
                sage: Po = PolynomialSpecies(QQ, "X")
                sage: F = Po(CyclicPermutationGroup(4))._compose_with_singletons(["X", "Y"], [[2, 2]])
                sage: sorted(str(M.permutation_group()) for M in F1.support()) == sorted(str(M.permutation_group()) for M in F.support())
                True

            TESTS::

                sage: P.one()._compose_with_singletons(["X"], [[0]])
                1
                sage: P.zero()._compose_with_singletons(["X"], [[0]])
                0

                sage: E2X._compose_with_singletons(["X", "Y"], [[1], [1]])
                Traceback (most recent call last):
                ...
                ValueError: the number of compositions should be the arity of self
                sage: E2X._compose_with_singletons(["X"], [[1, 1]])
                Traceback (most recent call last):
                ...
                ValueError: the total length of the compositions must match the number of names
                sage: E2X._compose_with_singletons(["X", "Y"], [[1, 2]])
                Traceback (most recent call last):
                ...
                ValueError: the size of the i-th composition should be the grade of self in sort i
                sage: (E2X + P(_wreath_group(2, 1).subgroup([])))._compose_with_singletons(["X", "Y"], [[1, 1]])
                Traceback (most recent call last):
                ...
                ValueError: self should be homogeneous with respect to all sorts
                sage: P.zero()._compose_with_singletons(["X"], [[1]])
                Traceback (most recent call last):
                ...
                ValueError: the size of the i-th composition should be the grade of self in sort i, which is zero
            """
            P = self.parent()
            r = P._r
            names = normalize_names(-1, names)
            if len(args) != P._arity:
                raise ValueError("the number of compositions should be the arity of self")
            if sum(len(c) for c in args) != len(names):
                raise ValueError("the total length of the compositions must match the number of names")
            Q = PolynomialHyperoctahedralSpecies(P.base_ring(), r, names)

            grades = set(M.grade() for M in self.support())
            if len(grades) > 1:
                raise ValueError("self should be homogeneous with respect to all sorts")
            if not grades:
                if any(sum(c) for c in args):
                    raise ValueError("the size of the i-th composition should be "
                                     "the grade of self in sort i, which is zero")
                return Q.zero()
            mc = next(iter(grades))
            if not all(sum(c) == d for c, d in zip(args, mc)):
                raise ValueError("the size of the i-th composition should be the grade of self in sort i")

            comp = list(chain.from_iterable(args))
            # the group of colour preserving relabellings of the
            # standard domain, together with the assignment of its
            # C_r-orbits to the sorts of the result
            S_down = _wreath_young_subgroup(r, comp)
            dompart = _canonical_dompart(r, comp)
            domain = range(1, r * sum(comp) + 1)
            # the group of sort preserving relabellings of the
            # standard graded r-set of grade mc
            S_up = _wreath_young_subgroup(r, mc)

            result = Q.zero()
            for M, c in self:
                # the group of M acts on the standard graded r-set of
                # grade mc, its sorts being consecutive
                H, _ = M.permutation_group()
                taus = libgap.DoubleCosetRepsAndSizes(S_up, S_down, H)

                # sum over double coset representatives: the group of
                # the coloured structure with representative tau is the
                # intersection tau H tau^-1 n S_down
                summand = Q.zero()
                for tau, _ in taus:
                    K = libgap.Intersection(libgap.ConjugateGroup(H, tau.Inverse()),
                                           S_down)
                    summand += Q(PermutationGroup(gap_group=K, domain=domain),
                                 dompart)
                result += c * summand
            return result

        def hadamard_product(self, other):
            r"""
            Return the Hadamard product of ``self`` and ``other``.

            The Hadamard product `F\boxtimes G` of two `r`-species is
            given by `(F\boxtimes G)[U] = F[U]\times G[U]` on every
            graded `C_r`-set `U`.  On molecular terms it is computed
            with the double coset formula

            .. MATH::

                \frac{X^{\mathbf d}}{H}\boxtimes\frac{X^{\mathbf d}}{K}
                = \sum_{\tau\in H\backslash W(r;\mathbf d)/K}
                \frac{X^{\mathbf d}}{H\cap\tau K\tau^{-1}},

            where `W(r;\mathbf d)` is the group of sort preserving
            relabellings of the standard graded `C_r`-set of grade
            `\mathbf d`.  The pairs of structures are indexed by the
            double cosets, the stabilizer of the pair corresponding to
            `\tau` being the intersection of the stabilizers.

            EXAMPLES:

            Exercise 2.1.9 in [BLL1998]_::

                sage: from sage.rings.species_hyperoctahedral import PolynomialHyperoctahedralSpecies
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
                sage: P = PolynomialHyperoctahedralSpecies(QQ, 1)
                sage: W3 = _wreath_group(1, 3)
                sage: C3 = P(W3.subgroup([(1, 2, 3)]))
                sage: X3 = P(W3.subgroup([]))
                sage: XE2 = P(W3.subgroup([(2, 3)]))
                sage: E2 = P(_wreath_group(1, 2))
                sage: C3.hadamard_product(C3)
                2*C_3(X)
                sage: X3.hadamard_product(C3)
                2*X^3
                sage: XE2.hadamard_product(XE2)
                X*E_2(X) + X^3
                sage: C3.hadamard_product(E2**2)
                0

            The table of marks of the molecular `2`-species of degree
            two::

                sage: P = PolynomialHyperoctahedralSpecies(QQ, 2)
                sage: W = _wreath_group(2, 2)
                sage: C = [P(W.subgroup([])), P(W.subgroup([(1, 2), (3, 4)])), P(W)]
                sage: table([(b, [(a.hadamard_product(b)).coefficient(b.support()[0]) for a in C]) for b in C])
                  X^2      [8, 2, 1]
                  X°^2     [0, 2, 1]
                  E_2(X°)  [0, 0, 1]

            Structures with different grades do not interact, as shown
            above for `C_3(X)` and `E_2(X)^2`.

            A multisort example::

                sage: Q = PolynomialHyperoctahedralSpecies(QQ, 2, "X, Y")
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_young_subgroup
                sage: W11 = _wreath_young_subgroup(2, [1, 1])
                sage: XY = Q(W11.subgroup([]), {0: [1, 2], 1: [3, 4]})
                sage: XoYo = Q(W11, {0: [1, 2], 1: [3, 4]})
                sage: XY.hadamard_product(XY)
                4*X*Y
                sage: XY.hadamard_product(XoYo)
                X*Y

            TESTS::

                sage: C3.hadamard_product(-C3)
                -2*C_3(X)

                sage: PolynomialHyperoctahedralSpecies(QQ, 2).one().hadamard_product(PolynomialHyperoctahedralSpecies(QQ, 2, "Y").one())
                Traceback (most recent call last):
                ...
                ValueError: the factors of a Hadamard product must have the same parent
            """
            P = self.parent()
            if P is not other.parent():
                raise ValueError("the factors of a Hadamard product must have the same parent")
            r = P._r

            result = P.zero()
            for M, c in self:
                mc = M.grade()
                domain = range(1, r * sum(mc) + 1)
                # the group of sort preserving relabellings of the
                # standard graded r-set of grade mc, together with
                # the assignment of its C_r-orbits to the sorts
                S_up = _wreath_young_subgroup(r, list(mc))
                dompart = _canonical_dompart(r, tuple(mc))
                # the group of M acts on the standard graded r-set of
                # grade mc, its sorts being consecutive
                H, _ = M.permutation_group()

                for N, d in other:
                    if mc != N.grade():
                        continue
                    K, _ = N.permutation_group()

                    # sum over double coset representatives: the group
                    # of the pair of structures with representative tau
                    # is the intersection H n tau K tau^-1
                    summand = P.zero()
                    for tau, _ in libgap.DoubleCosetRepsAndSizes(S_up, H, K):
                        G = libgap.Intersection(libgap.ConjugateGroup(K, tau.Inverse()),
                                               H)
                        summand += P(PermutationGroup(gap_group=G, domain=domain),
                                     dompart)
                    result += c * d * summand
            return result

        def _compose_with_weighted_singletons(self, names, multiplicities, args):
            r"""
            Return the type 2 substitution of ``self`` with sums of
            weighted singleton species.

            This is the weighted version of
            :meth:`_compose_with_singletons`: the homogeneous `k`-sort
            `r`-species ``self`` is substituted with
            `(\sum_j m_{1,j} X_{1,j}, \ldots, \sum_j m_{k,j} X_{k,j})`.

            The entries of ``multiplicities`` are weights, elements of
            the base ring, not necessarily integer multiplicities.  The
            unweighted composition with the corresponding singleton
            species is computed first; the result is then combined by a
            Hadamard product with the product of the weighted
            exponentials

            .. MATH::

                \prod_i E_r\left(\sum_j m_{i,j} X_{i,j}\right),

            so that the usual lambda-ring structure of the base ring is
            used.  In particular, if `q` has lambda degree one, then
            `\psi_n(q) = q^n` and, writing `R = \QQ[q]`, the theory of
            [Henderson2004]_ yields

            .. MATH::

                E_r((1+q)X) = (1+q^2)E(r)_2 + q(X^\circ)^2.

            INPUT:

            - ``names`` -- the (flat) list of the names of the sorts of
              the result, whose length is the total number of parts of
              the compositions in ``args``

            - ``multiplicities`` -- a (flat) list of weights, elements
              of the base ring, of the same length as ``names``

            - ``args`` -- a sequence of `k` compositions, where `k` is
              the arity of ``self``.  The parts of the `i`-th
              composition sum to the grade of ``self`` in sort `i` and
              count `C_r`-orbits.

            EXAMPLES:

            The `r`-analogue of Equation (2.5.41) in [BLL1998]_, with
            `E(r)_2 = E_2(X^\circ)` and `X^\circ` the degree one
            species with stabilizer `C_r`::

                sage: from sage.rings.species_hyperoctahedral import PolynomialHyperoctahedralSpecies
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group
                sage: P = PolynomialHyperoctahedralSpecies(QQ, 2)
                sage: W = _wreath_group(2, 2)
                sage: P(W)._compose_with_weighted_singletons(["X"], [-1], [[2]])
                -E_2(X°) + X°^2

                sage: E2X = P(W.subgroup([[(1, 3), (2, 4)]]))
                sage: E2X._compose_with_weighted_singletons(["X"], [-1], [[2]])
                -E_2(X) + X^2

            The weighted exponential of [Henderson2004]_, (4.8), is recovered
            by substituting into the set-like species::

                sage: R.<q> = QQ[]
                sage: P = PolynomialHyperoctahedralSpecies(R, 2)
                sage: P(W)._compose_with_weighted_singletons(["X"], [1+q], [[2]])
                (q^2+1)*E_2(X°) + q*X°^2
                sage: P(_wreath_group(2, 3))._compose_with_weighted_singletons(["X"], [1+q], [[3]])
                (q^3+1)*E_3(X°) + (q^2+q)*E_2(X°)*X°

            The weights of the sorts are independent; here the two
            singleton structures of `E_2(X^\circ)` carry the weights
            `1` and `q`::

                sage: P(W)._compose_with_weighted_singletons(["X", "Y"], [1, q], [[1, 1]])
                q*X°*Y°
                sage: E2X = P(W.subgroup([[(1, 3), (2, 4)]]))
                sage: E2X._compose_with_weighted_singletons(["X", "Y"], [1, q], [[1, 1]])
                q*X*Y

            For `r = 1` the result agrees with the ordinary weighted
            composition, Equation (2.5.41) in [BLL1998]_::

                sage: P1 = PolynomialHyperoctahedralSpecies(QQ, 1)
                sage: E2 = P1(SymmetricGroup(2))
                sage: C4 = P1(CyclicPermutationGroup(4))
                sage: E2._compose_with_weighted_singletons(["X"], [-1], [[2]])
                -E_2(X) + X^2
                sage: C4._compose_with_weighted_singletons(["X", "Y"], [1, -1], [[2, 2]])
                2*X^2*Y^2 - {((1,2)(3,4),): ({1, 2}, {3, 4})}

                sage: (C4 + E2^2)._compose_with_weighted_singletons(["X"], [-1], [[4]])
                -C_4(X) + {((1,2)(3,4),)} + E_2(X)^2 - 2*E_2(X)*X^2 + X^4

            TESTS::

                sage: P.zero()._compose_with_weighted_singletons(["X"], [-1], [[0]])
                0

                sage: C4._compose_with_weighted_singletons(["X"], [-1, 0], [[4]])
                Traceback (most recent call last):
                ...
                ValueError: the number of names must match the number of multiplicities

                sage: C4._compose_with_weighted_singletons(["X"], [-1], [[2,2]])
                Traceback (most recent call last):
                ...
                ValueError: the total length of the compositions must match the number of names
            """
            if len(names) != len(multiplicities):
                raise ValueError("the number of names must match the number of multiplicities")
            if sum(len(c) for c in args) != len(names):
                raise ValueError("the total length of the compositions must match the number of names")
            left = self._compose_with_singletons(names, args)
            right = left.parent()._exponential(multiplicities,
                                               list(chain.from_iterable(args)))
            return left.hadamard_product(right)

        def __call__(self, *args):
            r"""
            Return the type `2` substitution of ``args`` into ``self``.

            This implements the composition `F \circ_2 (G_1, \ldots,
            G_k)` of [Henderson2004]_, Equation (4.8), where
            ``self`` is a `C_r`-equivariant species and each `G_i` is an
            ordinary species.

            The args may be multisort, in which case the composite is a
            `C_r`-equivariant species with the sorts of the args, as in
            the ordinary composition of species
            (:class:`~sage.rings.species.PolynomialSpecies`): the sort of
            a `C_r`-orbit of the composite is the sort of the
            corresponding point of the inner structure attached to the
            block.

            The result is a `C_r`-equivariant species with the sorts of
            the args over their base ring; the coefficients of ``self``
            have to coerce into it.

            EXAMPLES::

                sage: from sage.rings.species_hyperoctahedral import PolynomialHyperoctahedralSpecies, MolecularHyperoctahedralSpecies
                sage: from sage.rings.species import PolynomialSpecies
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group, _wreath_young_subgroup
                sage: PoX = PolynomialSpecies(QQ, "X")
                sage: X = PoX(SymmetricGroup(1))
                sage: E2 = PoX(SymmetricGroup(2))
                sage: E3 = PoX(SymmetricGroup(3))
                sage: P = PolynomialHyperoctahedralSpecies(QQ, 2)

            Substituting `E_2` into the cyclic singleton `X^\circ` gives
            the cyclic composition; substituting into the free singleton
            forgets the `C_r`-action of the inner orbits.  The atoms are
            only canonical up to the choice of generators, so we list
            the orders of the molecular groups where the display would
            depend on that choice::

                sage: Xo = P(_wreath_group(2, 1))
                sage: Xo(E2)
                X°(E_2)
                sage: [(M.permutation_group()[0].order(), c) for M, c in Xo(E3)]
                [(12, 1)]
                sage: X_free = P(_wreath_group(2, 1).subgroup([]))
                sage: X_free(E2)
                E_2(X)
                sage: X_free(X + X^2)
                X + X^2

            Substituting `2 E_2` into the set-like species `E_2(X^\circ)`::

                sage: E2Xo = P(_wreath_young_subgroup(2, [2]))
                sage: sorted((M.permutation_group()[0].order(), c) for M, c in E2Xo(2*E2))
                [(16, 1), (32, 2)]

            Weighted species are supported::

                sage: R.<q> = QQ[]
                sage: PoXq = PolynomialSpecies(R, "X")
                sage: Pq = PolynomialHyperoctahedralSpecies(R, 2)
                sage: E2Xo_q = Pq(_wreath_young_subgroup(2, [2]))
                sage: E2Xo_q((1+q)*PoXq(SymmetricGroup(1)))
                (q^2+1)*E_2(X°) + q*X°^2

            Substitution is linear in the outer species.  Each sort of
            a multisort species is substituted with the corresponding
            argument.  The two sorts of `E_2(X^\circ, Y^\circ)` below
            carry `E_2` and `E_3`-structures::

                sage: PXY = PolynomialHyperoctahedralSpecies(QQ, 2, "X, Y")
                sage: MXY = MolecularHyperoctahedralSpecies(2, "X, Y")
                sage: W22 = _wreath_young_subgroup(2, [2])
                sage: A = PXY(MXY(W22, {0: [1, 2, 3, 4]}))
                sage: B = PXY(MXY(W22, {1: [1, 2, 3, 4]}))
                sage: sorted((tuple(M.grade()), M.permutation_group()[0].order())
                ....:     for M, c in (A + B)(E2, E3))
                [((4,), 32), ((6,), 288)]

            For `r = 1` the result agrees with the ordinary composition of
            species::

                sage: P1 = PolynomialHyperoctahedralSpecies(QQ, 1)
                sage: P1(CyclicPermutationGroup(4))(E2)
                {((7,8), (1,3,5,7)(2,4,6,8))}
                sage: C4 = PoX(CyclicPermutationGroup(4))
                sage: C4(E2)
                {((7,8), (1,3,5,7)(2,4,6,8))}

            This holds also for multisort args: substituting the product of
            the singletons of two sorts into `E_2` yields the atomic
            species of perfect matchings between the two sorts::

                sage: PoXY = PolynomialSpecies(QQ, "X, Y")
                sage: X2 = PoXY(SymmetricGroup(1), {0: [1]})
                sage: Y2 = PoXY(SymmetricGroup(1), {1: [1]})
                sage: P1(SymmetricGroup(2))(X2 * Y2)
                {((1,2)(3,4),): ({1, 2}, {3, 4})}
                sage: E2(X2 * Y2)
                E_2(X*Y)

            The composite lives over the base ring of the args, and the
            coefficients of ``self`` have to coerce into it::

                sage: E2Xo(PoXq(SymmetricGroup(1))).parent().base_ring()
                Univariate Polynomial Ring in q over Rational Field
                sage: (q*E2Xo_q)(E2)
                Traceback (most recent call last):
                ...
                ValueError: the coefficients of self must coerce into the base ring of the args

            Substituting the singleton `X` forgets the sorts, and a zero
            argument annihilates all structures using it::

                sage: E2Xo(X)
                E_2(X°)
                sage: (Xo + E2Xo)(X)
                X° + E_2(X°)
                sage: (Xo + E2Xo)(PoX.zero())
                0

            Check the case of arity zero::

                sage: P0 = PolynomialHyperoctahedralSpecies(QQ, 2, [])
                sage: P0.one()()
                1
                sage: (5*P0.one())()
                5

            TESTS::

                sage: E2Xo()
                Traceback (most recent call last):
                ...
                ValueError: number of args must match arity of self
                sage: E2Xo(2)
                Traceback (most recent call last):
                ...
                ValueError: all args must be ordinary polynomial species
                sage: E2Xo(E2, E2)
                Traceback (most recent call last):
                ...
                ValueError: number of args must match arity of self
                sage: PoXY = PolynomialSpecies(QQ, "X, Y")
                sage: X2 = PoXY(SymmetricGroup(1), {0: [1]})
                sage: Y2 = PoXY(SymmetricGroup(1), {1: [1]})

            The args may be multisort; the composite then has their
            sorts.  Substituting the product of the singletons of two
            sorts into `E_2(X^\circ)` gives a structure with four
            orbits, two of each sort::

                sage: E2Xo(X2 * Y2)
                {((3,4)(7,8), (1,2)(5,6), (1,3)(2,4)(5,7)(6,8)): ({1, 2, 3, 4}, {5, 6, 7, 8})}
                sage: [(tuple(M.grade()), M.permutation_group()[0].order()) for M, c in E2Xo(X2 * Y2)]
                [((2, 2), 8)]
                sage: A(E2, X2 * Y2)
                Traceback (most recent call last):
                ...
                ValueError: all args must have the same parent
                sage: E2Xo(PoX.one())
                Traceback (most recent call last):
                ...
                ValueError: all args must have positive degree
            """
            from sage.rings.species import PolynomialSpecies

            P = self.parent()
            if len(args) != P._arity:
                raise ValueError("number of args must match arity of self")
            if not all(isinstance(arg, PolynomialSpecies.Element)
                       for arg in args):
                raise ValueError("all args must be ordinary polynomial species")
            if len(set(arg.parent() for arg in args)) > 1:
                raise ValueError("all args must have the same parent")
            if any(sum(M.grade()) == 0 for g in args for M, _ in g):
                raise ValueError("all args must have positive degree")

            r = P._r
            # the composite is a C_r-species with the sorts and over
            # the base ring of the args
            P0 = args[0].parent() if args else P
            R = P0.base_ring()
            if P.base_ring() is not R:
                # the composite carries the weights of self and of the
                # args, so the coefficients of self must coerce into the
                # base ring of the args
                try:
                    terms = [(M, R(c)) for M, c in self]
                except (TypeError, ValueError):
                    raise ValueError("the coefficients of self must coerce "
                                     "into the base ring of the args")
                P = PolynomialHyperoctahedralSpecies(
                    R, r, P._indices._indices._names)
                self = P.sum_of_terms(terms)
            H = PolynomialHyperoctahedralSpecies(
                R, r, P0._indices._indices._names)
            if not self.support():
                return H.zero()
            if P._arity == 0:
                return H.sum_of_terms((H._indices.one(), c) for _, c in self)

            arg_terms = [sorted(g, key=lambda x: x[0].grade()) for g in args]
            multiplicities = list(chain.from_iterable(
                [[c for _, c in g] for g in arg_terms]))
            substituted = list(chain.from_iterable(
                [[M for M, _ in g] for g in arg_terms]))
            F_degrees = sorted(set(M.grade() for M, _ in self))
            names = ["X%s" % i for i in range(len(substituted))]

            result = H.zero()
            for mc in F_degrees:
                F = P.sum_of_terms((M, c) for M, c in self if M.grade() == mc)
                for degrees in cartesian_product(
                        [IntegerVectors(d, length=len(arg))
                         for d, arg in zip(mc, arg_terms)]):
                    FX = F._compose_with_weighted_singletons(names,
                                                             multiplicities,
                                                             degrees)
                    for M, c in FX:
                        result += c * H.monomial(
                            M._type2_substitute_molecular(substituted))
            return result

        def structures(self, *labels):
            r"""
            Iterate over the structures on the given free `C_r`-set of labels.

            The labels are given as one list of `C_r`-orbits per sort,
            each orbit an iterable of exactly `r` labels in cyclic phase
            order, see :meth:`AtomicHyperoctahedralSpecies.Element.structures
            <structures>`.

            This yields pairs consisting of a molecular `r`-species and a
            structure for it, and triples with an additional index between
            `0` and the coefficient of the molecule, if that coefficient is
            larger than one.

            EXAMPLES::

                sage: from sage.rings.species_hyperoctahedral import PolynomialHyperoctahedralSpecies
                sage: from sage.groups.perm_gps.hyperoctahedral_group import _wreath_group, _wreath_young_subgroup
                sage: P = PolynomialHyperoctahedralSpecies(QQ, 2)
                sage: Xo = P(_wreath_group(2, 1))
                sage: X = P(_wreath_group(2, 1).subgroup([]))
                sage: f = Xo + 2*X
                sage: sorted(f.structures([('a', 'b')]), key=str)
                [(X, (('a', 'b'),), 0),
                 (X, (('a', 'b'),), 1),
                 (X, (('b', 'a'),), 0),
                 (X, (('b', 'a'),), 1),
                 (X°, (('a', 'b'),))]

                sage: E2Xo = P(_wreath_young_subgroup(2, [2]))
                sage: list((Xo + E2Xo).structures([('a', 'b'), ('c', 'd')]))
                [(E_2(X°), (('a', 'b', 'c', 'd'),))]

            TESTS::

                sage: list((-Xo).structures([('a', 'b')]))
                Traceback (most recent call last):
                ...
                NotImplementedError: only implemented for proper non-virtual species
            """
            labels = _orbit_label_sets(self.parent()._arity,
                                       self.parent()._r, labels)
            for M, c in self.monomial_coefficients().items():
                if c not in ZZ or c < 0:
                    raise NotImplementedError("only implemented for proper "
                                              "non-virtual species")
                if c == 1:
                    for s in M.structures(*labels):
                        yield M, s
                else:
                    for e, s in cartesian_product([range(c),
                                                    M.structures(*labels)]):
                        yield M, s, e
