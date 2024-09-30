# sage_setup: distribution = sagemath-objects
"""
The C3 algorithm, under control of a total order

Abstract
========

Python handles multiple inheritance by computing, for each class,
a linear extension of the poset of all its super classes (the Method
Resolution Order, MRO). The MRO is calculated recursively from local
information (the *ordered* list of the direct super classes), with
the so-called ``C3`` algorithm. This algorithm can fail if the local
information is not consistent; worst, there exist hierarchies of
classes with provably no consistent local information.

For large hierarchy of classes, like those derived from categories in
Sage, maintaining consistent local information by hand does not scale
and leads to unpredictable ``C3`` failures (the dreaded "could not
find a consistent method resolution order"); a maintenance nightmare.

This module implements a final solution to this problem. Namely, it
allows for building automatically the local information from the bare
class hierarchy in such a way that guarantees that the ``C3``
algorithm will never fail.

Err, but you said that this was provably impossible? Well, not if one
relaxes a bit the hypotheses; but that's not something one would want
to do by hand :-)

The problem
===========

Consider the following hierarchy of classes::

    sage: class A1(): pass
    sage: class A2():
    ....:     def foo(self): return 2
    sage: class A3(): pass
    sage: class A4():
    ....:     def foo(self): return 4
    sage: class A5(A2, A1):
    ....:     def foo(self): return 5
    sage: class A6(A4, A3): pass
    sage: class A7(A6, A5): pass

If ``a`` is an instance of ``A7``, then Python needs to choose which
implementation to use upon calling ``a.foo()``: that of ``A4`` or
``A5``, but obviously not that of ``A2``. In Python, like in many
other dynamic object oriented languages, this is achieved by
calculating once for all a specific linear extension of the hierarchy
of the super classes of each class, called its Method Resolution Order
(MRO)::

    sage: [cls.__name__ for cls in A7.mro()]
    ['A7', 'A6', 'A4', 'A3', 'A5', 'A2', 'A1', 'object']

Thus, in our example, the implementation in ``A4`` is chosen::

    sage: a = A7()
    sage: a.foo()
    4

Specifically, the MRO is calculated using the so-called ``C3``
algorithm which guarantees that the MRO respects not only inheritance,
but also the order in which the bases (direct super classes) are given
for each class.

However, for large hierarchies of classes with lots of multiple
inheritance, like those derived from categories in Sage, this
algorithm easily fails if the order of the bases is not chosen
consistently (here for ``A2`` w.r.t. ``A1``)::

    sage: class B6(A1,A2): pass
    sage: class B7(B6,A5): pass
    Traceback (most recent call last):
    ...
    TypeError: Cannot create a consistent method resolution
    order (MRO) for bases A1, A2

There actually exist hierarchies of classes for which ``C3`` fails
whatever order of the bases is chosen; the smallest such example,
admittedly artificial, has ten classes (see below). Still, this
highlights that this problem has to be tackled in a systematic way.

Fortunately, one can trick ``C3``, without changing the inheritance
semantic, by adding some super classes of ``A`` to the bases of
``A``. In the following example, we completely force a given MRO by
specifying *all* the super classes of ``A`` as bases::

    sage: class A7(A6, A5, A4, A3, A2, A1): pass
    sage: [cls.__name__ for cls in A7.mro()]
    ['A7', 'A6', 'A5', 'A4', 'A3', 'A2', 'A1', 'object']

Luckily this can be optimized; here it is sufficient to add a single
base to enforce the same MRO::

    sage: class A7(A6, A5, A4): pass
    sage: [cls.__name__ for cls in A7.mro()]
    ['A7', 'A6', 'A5', 'A4', 'A3', 'A2', 'A1', 'object']

A strategy to solve the problem
===============================

We should recall at this point a design decision that we took for the
hierarchy of classes derived from categories: *the semantic shall only
depend on the inheritance order*, not on the specific MRO, and in
particular not on the order of the bases (see
:ref:`On the order of super categories <category-primer-category-order>`).

If a choice needs to be made (for example for efficiency reasons),
then this should be done explicitly, on a method-by-method basis. In
practice this design goal is not yet met.

.. NOTE::

    When managing large hierarchies of classes in other contexts this
    may be too strong a design decision.

The strategy we use for hierarchies of classes derived from categories
is then:

1. To choose a global total order on the whole hierarchy of classes.
2. To control ``C3`` to get it to return MROs that follow this total order.

A basic approach for point 1., that will work for any hierarchy of
classes, is to enumerate the classes while they are constructed
(making sure that the bases of each class are enumerated before that
class), and to order the classes according to that enumeration. A more
conceptual ordering may be desirable, in particular to get
deterministic and reproducible results. In the context of Sage, this
is mostly relevant for those doctests displaying all the categories or
classes that an object inherits from.

Getting fine control on C3
==========================

This module is about point 2.

The natural approach would be to change the algorithm used by Python
to compute the MRO. However, changing Python's default algorithm just
for our needs is obviously not an option, and there is currently no
hook to customize specific classes to use a different algorithm.
Pushing the addition of such a hook into stock Python would take too
much time and effort.

Another approach would be to use the "adding bases" trick
straightforwardly, putting the list of *all* the super classes of a
class as its bases. However, this would have several drawbacks:

- It is not so elegant, in particular because it duplicates
  information: we already know through ``A5`` that ``A7`` is a
  subclass of ``A1``. This duplication could be acceptable in our
  context because the hierarchy of classes is generated automatically
  from a conceptual hierarchy (the categories) which serves as single
  point of truth for calculating the bases of each class.

- It increases the complexity of the calculation of the MRO with
  ``C3``. For example, for a linear hierarchy of classes, the
  complexity goes from `O(n^2)` to `O(n^3)` which is not acceptable.

- It increases the complexity of inspecting the classes. For example,
  the current implementation of the ``dir`` command in Python has no
  cache, and its complexity is linear in the number of maximal paths
  in the class hierarchy graph as defined by the bases. For a linear
  hierarchy, this is of complexity `O(p_n)` where `p_n` is the number
  of integer partitions of `n`, which is exponential. And indeed,
  running ``dir`` for a typical class like
  ``GradedHopfAlgebrasWithBasis(QQ).parent_class`` with ``37`` super
  classes took `18` seconds with this approach.

  Granted: this mostly affects the ``dir`` command and could be blamed
  on its current implementation. With appropriate caching, it could be
  reimplemented to have a complexity roughly linear in the number of
  classes in the hierarchy. But this won't happen any time soon in a
  stock Python.

This module refines this approach to make it acceptable, if not
seamless. Given a hierarchy and a total order on this hierarchy, it
calculates for each element of the hierarchy the smallest list of
additional bases that forces ``C3`` to return the desired MRO. This is
achieved by implementing an instrumented variant of the ``C3``
algorithm (which we call *instrumented* ``C3``) that detects when
``C3`` is about to take a wrong decision and adds one base to force
the right decision. Then, running the standard ``C3`` algorithm with
the updated list of bases (which we call *controlled* ``C3``) yields
the desired MRO.

EXAMPLES:

As an experimentation and testing tool, we use a class
:class:`HierarchyElement` whose instances can be constructed from a
hierarchy described by a poset, a digraph, or more generally a
successor relation. By default, the desired MRO is sorted
decreasingly. Another total order can be specified using a sorting
key.

We consider the smallest poset describing a class hierarchy admitting
no MRO whatsoever::

    sage: P = Poset({10: [9,8,7], 9: [6,1], 8: [5,2], 7: [4,3],                         # needs sage.graphs
    ....:            6: [3,2], 5: [3,1], 4: [2,1]},
    ....:           linear_extension=True, facade=True)

And build a :class:`HierarchyElement` from it::

    sage: from sage.misc.c3_controlled import HierarchyElement
    sage: x = HierarchyElement(10, P)                                                   # needs sage.graphs

Here are its bases::

    sage: HierarchyElement(10, P)._bases                                                # needs sage.graphs
    [9, 8, 7]

Using the standard ``C3`` algorithm fails::

    sage: x.mro_standard                                                                # needs sage.graphs
    Traceback (most recent call last):
    ...
    ValueError: Cannot merge the items 3, 3, 2.

We also get a failure when we relabel `P` according to another linear
extension. For easy relabelling, we first need to set an appropriate
default linear extension for `P`::

    sage: linear_extension = list(reversed(IntegerRange(1, 11)))
    sage: P = P.with_linear_extension(linear_extension)                                 # needs sage.graphs
    sage: list(P)                                                                       # needs sage.graphs
    [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]

Now we play with a specific linear extension of `P`::

    sage: # needs sage.graphs
    sage: Q = P.linear_extension([10, 9, 8, 7, 6, 5, 4, 1, 2, 3]).to_poset()
    sage: Q.cover_relations()
    [[10, 9], [10, 8], [10, 7], [9, 6], [9, 3], [8, 5], [8, 2], [7, 4],
     [7, 1], [6, 2], [6, 1], [5, 3], [5, 1], [4, 3], [4, 2]]
    sage: x = HierarchyElement(10, Q)
    sage: x.mro_standard
    Traceback (most recent call last):
    ...
    ValueError: Cannot merge the items 2, 3, 3.

On the other hand, both the instrumented ``C3`` algorithm, and the
controlled ``C3`` algorithm give the desired MRO::

    sage: x.mro                                                                         # needs sage.graphs
    [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
    sage: x.mro_controlled                                                              # needs sage.graphs
    [10, 9, 8, 7, 6, 5, 4, 3, 2, 1]

The above checks, and more, can be run with::

    sage: x._test_mro()                                                                 # needs sage.graphs

In practice, the control was achieved by adding the following bases::

    sage: x._bases                                                                      # needs sage.graphs
    [9, 8, 7]
    sage: x._bases_controlled                                                           # needs sage.graphs
    [9, 8, 7, 6, 5]

Altogether, four bases were added for control::

    sage: sum(len(HierarchyElement(q, Q)._bases) for q in Q)                            # needs sage.graphs
    15
    sage: sum(len(HierarchyElement(q, Q)._bases_controlled) for q in Q)                 # needs sage.graphs
    19

This information can also be recovered with::

    sage: x.all_bases_len()                                                             # needs sage.graphs
    15
    sage: x.all_bases_controlled_len()                                                  # needs sage.graphs
    19

We now check that the ``C3`` algorithm fails for all linear extensions
`l` of this poset, whereas both the instrumented and controlled ``C3``
algorithms succeed; along the way, we collect some statistics::

    sage: L = P.linear_extensions()                                                     # needs sage.graphs
    sage: stats = []
    sage: for l in L:                                                                   # needs sage.graphs sage.modules
    ....:     x = HierarchyElement(10, l.to_poset())
    ....:     try: # Check that x.mro_standard always fails with a ValueError
    ....:         x.mro_standard
    ....:     except ValueError:
    ....:         pass
    ....:     else:
    ....:         assert False
    ....:     assert x.mro            == list(P)
    ....:     assert x.mro_controlled == list(P)
    ....:     assert x.all_bases_len() == 15
    ....:     stats.append(x.all_bases_controlled_len()-x.all_bases_len())

Depending on the linear extension `l` it was necessary to add between
one and five bases for control; for example, `216` linear extensions
required the addition of four bases::

    sage: sorted(Word(stats).evaluation_sparse())                                       # needs sage.combinat sage.graphs sage.modules
    [(1, 36), (2, 108), (3, 180), (4, 216), (5, 180)]

We now consider a hierarchy of categories::

    sage: from operator import attrgetter
    sage: x = HierarchyElement(Groups(), attrcall("super_categories"), attrgetter("_cmp_key"))
    sage: x.mro
    [Category of groups, Category of monoids, Category of semigroups,
     Category of inverse unital magmas, Category of unital magmas, Category of magmas,
     Category of sets, Category of sets with partial maps, Category of objects]
    sage: x.mro_standard
    [Category of groups, Category of monoids, Category of semigroups,
     Category of inverse unital magmas, Category of unital magmas, Category of magmas,
     Category of sets, Category of sets with partial maps, Category of objects]

For a typical category, few bases, if any, need to be added to force
``C3`` to give the desired order::

    sage: C = FiniteFields()
    sage: x = HierarchyElement(C, attrcall("super_categories"), attrgetter("_cmp_key"))
    sage: x.mro == x.mro_standard
    False
    sage: x.all_bases_len()
    72
    sage: x.all_bases_controlled_len()
    76

    sage: C = GradedHopfAlgebrasWithBasis(QQ)
    sage: x = HierarchyElement(C, attrcall("super_categories"), attrgetter("_cmp_key"))
    sage: x._test_mro()
    sage: x.mro == x.mro_standard
    False
    sage: x.all_bases_len()
    114
    sage: x.all_bases_controlled_len()
    117

The following can be used to search through the Sage named categories
for any that requires the addition of some bases. The output may
change a bit when the category hierarchy is changed. As long as the
list below does not change radically, it's fine to just update this
doctest::

    sage: from sage.categories.category import category_sample
    sage: sorted([C for C in category_sample()                                          # needs sage.combinat sage.graphs sage.modules sage.rings.number_field
    ....:         if len(C._super_categories_for_classes) != len(C.super_categories())],
    ....:        key=str)
    [Category of affine Weyl groups,
     Category of fields,
     Category of finite Weyl groups,
     Category of finite dimensional Hopf algebras with basis over Rational Field,
     Category of finite dimensional algebras with basis over Rational Field,
     Category of finite enumerated permutation groups,
     Category of number fields]

AUTHOR:

- Nicolas M. Thiery (2012-09): initial version.
"""
# ****************************************************************************
#  Copyright (C) 2012-2013  Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall
from sage.misc.cachefunc import cached_function, cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.dynamic_class import dynamic_class

##############################################################################
# Implementation of the total order between categories
##############################################################################

cdef tuple atoms = ("FacadeSets",
                    "FiniteSets", "Sets.Infinite", "EnumeratedSets", "SetsWithGrading",
                    "Posets", "LatticePosets", "Crystals", "AdditiveMagmas",
                    "FiniteDimensionalModules", "GradedModules", "ModulesWithBasis",
                    "Magmas", "Semigroups", "Monoids", "PermutationGroups",
                    "MagmasAndAdditiveMagmas", "Rngs", "Domains", "HopfAlgebras")


cdef dict flags = { atom: 1 << i for i,atom in enumerate(atoms) }

cdef class CmpKey:
    r"""
    This class implements the lazy attribute ``Category._cmp_key``.

    The comparison key ``A._cmp_key`` of a category is used to define
    an (almost) total order on non-join categories by setting, for two
    categories `A` and `B`, `A<B` if ``A._cmp_key > B._cmp_key``. This
    order in turn is used to give a normal form to join's, and help
    toward having a consistent method resolution order for
    parent/element classes.

    The comparison key should satisfy the following properties:

    - If `A` is a subcategory of `B`, then `A < B` (so that
      ``A._cmp_key > B._cmp_key``). In particular,
      :class:`Objects() <Objects>` is the largest category.

    - If `A != B` and taking the join of `A` and `B` makes sense
      (e.g. taking the join of ``Algebras(GF(5))`` and
      ``Algebras(QQ)`` does not make sense), then `A<B` or `B<A`.

    The rationale for the inversion above between `A<B` and
    ``A._cmp_key > B._cmp_key`` is that we want the order to
    be compatible with inclusion of categories, yet it's easier in
    practice to create keys that get bigger and bigger while we go
    down the category hierarchy.

    This implementation applies to join-irreducible categories
    (i.e. categories that are not join categories). It returns a
    pair of integers ``(flags, i)``, where ``flags`` is to be
    interpreted as a bit vector. The first bit is set if ``self``
    is a facade set. The second bit is set if ``self`` is finite.
    And so on. The choice of the flags is adhoc and was primarily
    crafted so that the order between categories would not change
    too much upon integration of :issue:`13589` and would be
    reasonably session independent. The number ``i`` is there
    to resolve ambiguities; it is session dependent, and is
    assigned increasingly when new categories are created.

    .. NOTE::

        This is currently not implemented using a
        :class:`lazy_attribute` for speed reasons only (the code is in
        Cython and takes advantage of the fact that Category objects
        always have a ``__dict__`` dictionary)

    .. TODO::

        - Handle nicely (covariant) functorial constructions and axioms

    EXAMPLES::

        sage: Objects()._cmp_key
        (0, 0)
        sage: SetsWithPartialMaps()._cmp_key
        (0, 1)
        sage: Sets()._cmp_key
        (0, 2)
        sage: Sets().Facade()._cmp_key
        (1, ...)
        sage: Sets().Finite()._cmp_key
        (2, ...)
        sage: Sets().Infinite()._cmp_key
        (4, ...)
        sage: EnumeratedSets()._cmp_key
        (8, ...)
        sage: FiniteEnumeratedSets()._cmp_key
        (10, ...)
        sage: SetsWithGrading()._cmp_key
        (16, ...)
        sage: Posets()._cmp_key
        (32, ...)
        sage: LatticePosets()._cmp_key
        (96, ...)
        sage: Crystals()._cmp_key
        (136, ...)
        sage: AdditiveMagmas()._cmp_key
        (256, ...)
        sage: Magmas()._cmp_key
        (4096, ...)
        sage: CommutativeAdditiveSemigroups()._cmp_key
        (256, ...)
        sage: Rings()._cmp_key
        (225536, ...)
        sage: Algebras(QQ)._cmp_key
        (225536, ...)
        sage: AlgebrasWithBasis(QQ)._cmp_key
        (227584, ...)
        sage: GradedAlgebras(QQ)._cmp_key
        (226560, ...)
        sage: GradedAlgebrasWithBasis(QQ)._cmp_key
        (228608, ...)

    For backward compatibility we currently want the following comparisons::

        sage: EnumeratedSets()._cmp_key > Sets().Facade()._cmp_key
        True
        sage: AdditiveMagmas()._cmp_key > EnumeratedSets()._cmp_key
        True

        sage: Category.join([EnumeratedSets(), Sets().Facade()]).parent_class._an_element_.__module__
        'sage.categories.enumerated_sets'

        sage: CommutativeAdditiveSemigroups()._cmp_key < Magmas()._cmp_key
        True
        sage: VectorSpaces(QQ)._cmp_key < Rings()._cmp_key
        True
        sage: VectorSpaces(QQ)._cmp_key < Magmas()._cmp_key
        True
    """
    cdef int count

    def __init__(self):
        """
        Set the internal category counter to zero.

        EXAMPLES::

            sage: Objects()._cmp_key    # indirect doctest
            (0, 0)
        """
        self.count = -1

    def __get__(self, object inst, object cls):
        """
        Bind the comparison key to the given instance.

        EXAMPLES::

            sage: C = Algebras(FractionField(QQ['x']))
            sage: C._cmp_key
            (225536, ...)
            sage: '_cmp_key' in C.__dict__    # indirect doctest
            True
        """
        # assert not isinstance(inst, JoinCategory)
        # Note that cls is a DynamicClassMetaclass, hence not a type
        cdef str classname = cls.__base__.__name__
        cdef int flag = flags.get(classname, 0)
        cdef object cat
        for cat in inst._super_categories:
            flag = flag | <int>(<tuple>(cat._cmp_key)[0])
        self.count += 1
        inst._cmp_key = (flag, self.count)
        return flag, self.count

_cmp_key = CmpKey()


cdef class CmpKeyNamed:
    """
    This class implements the lazy attribute ``CategoryWithParameters._cmp_key``.

    .. SEEALSO::

        - :class:`CmpKey`
        - :class:`lazy_attribute`
        - :class:`sage.categories.category.CategoryWithParameters`.

    .. NOTE::

        - The value of the attribute depends only on the parameters of
          this category.

        - This is currently not implemented using a
          :class:`lazy_attribute` for speed reasons only.

    EXAMPLES::

        sage: Algebras(GF(3))._cmp_key == Algebras(GF(5))._cmp_key  # indirect doctest
        True
        sage: Algebras(ZZ)._cmp_key != Algebras(GF(5))._cmp_key
        True
    """
    def __get__(self, object inst, object cls):
        """
        EXAMPLES::

            sage: Algebras(GF(3))._cmp_key == Algebras(GF(5))._cmp_key  # indirect doctest
            True
            sage: Algebras(ZZ)._cmp_key != Algebras(GF(5))._cmp_key
            True
        """
        cdef dict D = cls._make_named_class_cache
        cdef str name = "_cmp_key"
        cdef tuple key = (cls.__base__, name, inst._make_named_class_key(name))
        try:
            result = D[key]
            inst._cmp_key = result
            return result
        except KeyError:
            pass
        result = _cmp_key.__get__(inst,cls)
        D[key] = result
        return result

_cmp_key_named = CmpKeyNamed()


##############################################################################

def C3_merge(list lists):
    r"""
    Return the input lists merged using the ``C3`` algorithm.

    EXAMPLES::

        sage: from sage.misc.c3_controlled import C3_merge
        sage: C3_merge([[3,2],[4,3,1]])
        [4, 3, 2, 1]
        sage: C3_merge([[3,2],[4,1]])
        [3, 2, 4, 1]

    This function is only used for testing and experimenting purposes,
    but exercised quite some by the other doctests in this file.

    It is an extract of :func:`sage.misc.c3.C3_algorithm`; the latter
    could be possibly rewritten to use this one to avoid duplication.
    """
    cdef list out = []
    # Data structure / invariants:
    # We will be working with the MROs of the super objects
    # together with the list of bases of ``self``.
    # Each list is split between its head (in ``heads``) and tail (in
    # ``tails'') . Each tail is stored reversed, so that we can use a
    # cheap pop() in lieue of pop(0). A duplicate of the tail is
    # stored as a set in ``tailsets`` for cheap membership testing.
    # Since we actually want comparison by identity, not equality,
    # what we store is the set of memory locations of objects
    cdef object O, X
    cdef list tail, l
    cdef set tailset

    cdef list tails    = [l[::-1]              for l in lists if l]
    cdef list heads    = [tail.pop()             for tail in tails]
    cdef list tailsets = [set(O for O in tail) for tail in tails] # <size_t><void *>

    cdef int i, j, nbheads
    nbheads = len(heads)
    cdef bint next_item_found

    while nbheads:
        for i in range(nbheads):
            O = heads[i]
            # Does O appear in none of the tails?  ``all(O not in tail for tail in tailsets)``
            next_item_found = True
            for j in range(nbheads):
                if j == i:
                    continue
                tailset = tailsets[j]
                if O in tailset: # <size_t><void *>O
                    next_item_found = False
                    break
            if next_item_found:
                out.append(O)
                # Clear O from other heads, removing the line altogether
                # if the tail is already empty.
                # j goes down so that ``del heads[j]`` does not screw up the numbering
                for j in range(nbheads-1, -1, -1):
                    if heads[j] == O: # is O
                        tail = tails[j]
                        if tail:
                            X = tail.pop()
                            heads[j] = X
                            tailset = tailsets[j]
                            tailset.remove(X) # <size_t><void *>X)
                        else:
                            del heads[j]
                            del tails[j]
                            del tailsets[j]
                            nbheads -= 1
                break
        if not next_item_found:
            # No head is available
            raise ValueError("Cannot merge the items %s."%', '.join(repr(head) for head in heads))
    return out


cpdef identity(x):
    r"""
    EXAMPLES::

        sage: from sage.misc.c3_controlled import identity
        sage: identity(10)
        10
    """
    return x

cpdef tuple C3_sorted_merge(list lists, key=identity):
    r"""
    Return the sorted input lists merged using the ``C3`` algorithm, with a twist.

    INPUT:

    - ``lists`` -- a non empty list (or iterable) of lists (or
      iterables), each sorted strictly decreasingly according
      to ``key``
    - ``key`` -- a function

    OUTPUT: a pair ``(result, suggestion)``

    ``result`` is the sorted list obtained by merging the lists in
    ``lists`` while removing duplicates, and ``suggestion`` is a list
    such that applying ``C3`` on ``lists`` with its last list replaced
    by ``suggestion`` would return ``result``.

    EXAMPLES:

    With the following input, :func:`C3_merge` returns right away a
    sorted list::

        sage: from sage.misc.c3_controlled import C3_merge
        sage: C3_merge([[2],[1]])
        [2, 1]

    In that case, :func:`C3_sorted_merge` returns the same result,
    with the last line unchanged::

        sage: from sage.misc.c3_controlled import C3_sorted_merge
        sage: C3_sorted_merge([[2],[1]])
        ([2, 1], [1])

    On the other hand, with the following input, :func:`C3_merge`
    returns a non sorted list::

        sage: C3_merge([[1],[2]])
        [1, 2]

    Then, :func:`C3_sorted_merge` returns a sorted list, and suggests
    to replace the last line by ``[2,1]``::

        sage: C3_sorted_merge([[1],[2]])
        ([2, 1], [2, 1])

    And indeed :func:`C3_merge` now returns the desired result::

        sage: C3_merge([[1],[2,1]])
        [2, 1]

    From now on, we use this little wrapper that checks that
    :func:`C3_merge`, with the suggestion of :func:`C3_sorted_merge`,
    returns a sorted list::

        sage: def C3_sorted_merge_check(lists):
        ....:     result, suggestion = C3_sorted_merge(lists)
        ....:     assert result == C3_merge(lists[:-1] + [suggestion])
        ....:     return result, suggestion

    Base cases::

        sage: C3_sorted_merge_check([])
        Traceback (most recent call last):
        ...
        ValueError: The input should be a non empty list of lists (or iterables)
        sage: C3_sorted_merge_check([[]])
        ([], [])
        sage: C3_sorted_merge_check([[1]])
        ([1], [1])
        sage: C3_sorted_merge_check([[3,2,1]])
        ([3, 2, 1], [3, 2, 1])
        sage: C3_sorted_merge_check([[1],[1]])
        ([1], [1])
        sage: C3_sorted_merge_check([[3,2,1],[3,2,1]])
        ([3, 2, 1], [3, 2, 1])

    Exercise different states for the last line::

        sage: C3_sorted_merge_check([[1],[2],[]])
        ([2, 1], [2, 1])
        sage: C3_sorted_merge_check([[1],[2], [1]])
        ([2, 1], [2, 1])

    Explore (all?) the different execution branches::

        sage: C3_sorted_merge_check([[3,1],[4,2]])
        ([4, 3, 2, 1], [4, 3, 2, 1])
        sage: C3_sorted_merge_check([[4,1],[3,2]])
        ([4, 3, 2, 1], [3, 2, 1])
        sage: C3_sorted_merge_check([[3,2],[4,1]])
        ([4, 3, 2, 1], [4, 3, 1])
        sage: C3_sorted_merge_check([[1],[4,3,2]])
        ([4, 3, 2, 1], [4, 3, 2, 1])
        sage: C3_sorted_merge_check([[1],[3,2], []])
        ([3, 2, 1], [2, 1])
        sage: C3_sorted_merge_check([[1],[4,3,2], []])
        ([4, 3, 2, 1], [2, 1])
        sage: C3_sorted_merge_check([[1],[4,3,2], [2]])
        ([4, 3, 2, 1], [2, 1])
        sage: C3_sorted_merge_check([[2],[1],[4],[3]])
        ([4, 3, 2, 1], [3, 2, 1])
        sage: C3_sorted_merge_check([[2],[1],[4],[]])
        ([4, 2, 1], [4, 2, 1])
        sage: C3_sorted_merge_check([[2],[1],[3],[4]])
        ([4, 3, 2, 1], [4, 3, 2, 1])
        sage: C3_sorted_merge_check([[2],[1],[3,2,1],[3]])
        ([3, 2, 1], [3])
        sage: C3_sorted_merge_check([[2],[1],[2,1],[3]])
        ([3, 2, 1], [3, 2])

    Exercises adding one item when the last list has a single element;
    the second example comes from an actual poset::

        sage: C3_sorted_merge_check([[5,4,2],[4,3],[5,4,1]])
        ([5, 4, 3, 2, 1], [5, 4, 3, 2, 1])
        sage: C3_sorted_merge_check([[6,4,2],[5,3],[6,5,1]])
        ([6, 5, 4, 3, 2, 1], [6, 5, 4, 3, 2, 1])
    """
    lists = list(lists)
    if not lists:
        raise ValueError("The input should be a non empty list of lists (or iterables)")
    #for l in lists:
    #    assert sorted(l, key = key, reverse=True) == l,\
    #        "Each input list should be sorted %s"%l

    cdef set suggestion = set(lists[-1])
    cdef bint last_list_non_empty = bool(lists[-1])
    cdef list out = []
    # Data structure / invariants:
    # - Each list remains sorted and duplicate free.
    # - Each list only evolves by popping its largest element
    #   Exception: elements can be inserted back into the last list.
    # - Whenever a list becomes empty, it's removed from the data structure.
    #   The order between the (remaining non empty) lists remains unchanged.
    # - nbheads contains the number of lists appearing in the data structure.
    # - The flag ``last_list_non_empty`` states whether the last
    #   list is currently non empty; if yes, by the above, this list is stored last.
    # - Each list is split between its head (in ``heads``) and tail (in ``tails'').
    # - Each tail is stored reversed, so that we can use a cheap ``pop()``
    #   in lieue of ``pop(0)``.
    # - A duplicate of this tail is stored as a set (of keys) in
    #   ``tailsets``, for cheap membership testing.

    cdef int i, j, max_i
    cdef bint cont
    cdef list tail, l
    cdef set tailset

    cdef list tails    = [l[::-1]                 for l in lists if l]
    cdef list heads    = [tail.pop()                for tail in tails]
    cdef set tmp_set
    cdef list tailsets = [] # remove closure [set(key(O) for O in tail) for tail in tails]
    for tail in tails:
        tmp_set = set()
        for O in tail:
            tmp_set.add(key(O))
        tailsets.append(tmp_set)
    # for i in range(len(tails)):
    #     assert len(tails[i]) == len(tailsets[i]), \
    #         "All objects should be distinct and have distinct sorting key!"+'\n'.join(" - %s: %s"%(key(O), O) for O in sorted(tails[i], key=key))

    cdef int nbheads = len(heads)
    cdef dict holder = {}

    # def print_state():
    #     print("-- %s -- %s ------"%(out,suggestion))
    #     for i in range(nbheads):
    #         print([heads[i]] + list(reversed(tails[i])))

    # def check_state():
    #     for i in range(nbheads):
    #         l = tails[i]
    #         if heads[i] is not None:
    #             l = l+[heads[i]]
    #         assert sorted(l, key=key) == l
    #         assert len(set(l)) == len(l)
    #         assert len(tails[i]) == len(set(tails[i])), \
    #             "C3's input list should have no repeats %s"%tails[i]
    #         assert set(key(O) for O in tails[i]) == tailsets[i], \
    #             "inconsistent tails[i] and tailsets[i]: %s %s"%(tails[i], tailsets[i])
    #         assert len(tails[i]) == len(tailsets[i]), \
    #             "keys should be distinct"%(tails[i])

    while nbheads:
        #print_state()
        #check_state()
        # Find the position of the largest head which will become the next item
        max_i   = 0
        max_key = key(heads[0])
        for i in range(1, nbheads):
            O = heads[i]
            O_key = key(O)
            if O_key > max_key:
                max_i = i
                max_key = O_key
        max_value = heads[max_i]

        # Find all the bad choices
        max_bad = None
        for i in range(max_i):
            O = heads[i]
            # Does O appear in none of the tails?
            O_key = key(O)
            # replace the closure
            # if any(O_key in tailsets[j] for j in range(nbheads) if j != i): continue
            cont = False
            for j in range(i):
                if O_key in tailsets[j]:
                    cont = True
                    break
            if cont:
                continue
            for j from i<j<nbheads:
                if O_key in tailsets[j]:
                    cont = True
                    break
            if cont:
                continue

            # The plain C3 algorithm would have chosen O as next item!
            if max_bad is None or O_key > key(max_bad):
                max_bad = O

            # We prevent this choice by inserting O in the tail of the
            # suggestions. At this stage, we only insert it into the
            # last list. Later, we will make sure that it is actually
            # in the tail of the last list.
            if not last_list_non_empty:
                # Reinstate the last list for the suggestion
                # if it had disappeared before
                heads.append(O)
                tails.append([])
                tailsets.append(set())
                nbheads += 1
                last_list_non_empty = True
            elif O_key > key(heads[-1]):
                tails[-1].append(heads[-1])
                tailsets[-1].add(key(heads[-1]))
                heads[-1] = O
            elif O != heads[-1]:
                assert O_key not in tailsets[-1], "C3 should not have chosen this O"
                # Use a heap or something for fast sorted insertion?
                # Since Python uses TimSort, that's probably not so bad.
                tails[-1].append(O)
                tails[-1].sort(key = key)
                tailsets[-1].add(O_key)
            suggestion.add(O)
            #check_state()

        # Insert max_value in the last list, if needed to hold off the bad items
        if max_bad is not None:
            last_head = heads[-1]
            if last_head is None or key(max_bad) >= key(last_head):
                if last_head is not None and last_head != max_bad:
                    tails[-1].append(last_head)
                    tailsets[-1].add(key(last_head))
                    #check_state()
                heads[-1] = max_value
                holder[max_bad] = max_value
                #check_state()

        out.append(max_value)
        # Clear O from other heads, removing the line altogether
        # if the tail is already empty.
        # j goes down so that ``del heads[j]`` does not screw up the numbering
        for j in range(nbheads-1, -1, -1):
            if heads[j] == max_value:
                tail = tails[j]
                if tail:
                    X = tail.pop()
                    heads[j] = X
                    tailset = tailsets[j]
                    tailset.remove(key(X))
                else:
                    del heads[j]
                    del tails[j]
                    del tailsets[j]
                    nbheads -= 1
                    if last_list_non_empty and j == nbheads:
                        last_list_non_empty = False
                #check_state()
    suggestion.update(holder.values())
    cdef list suggestion_list = sorted(suggestion, key = key, reverse=True)
    #assert C3_merge(lists[:-1]+[suggestion_list]) == out
    return (out, suggestion_list)


class HierarchyElement(object, metaclass=ClasscallMetaclass):
    """
    A class for elements in a hierarchy.

    This class is for testing and experimenting with various variants
    of the ``C3`` algorithm to compute a linear extension of the
    elements above an element in a hierarchy. Given the topic at hand,
    we use the following naming conventions. For `x` an element of the
    hierarchy, we call the elements just above `x` its *bases*, and
    the linear extension of all elements above `x` its *MRO*.

    By convention, the bases are given as lists of instances of
    :class:`HierarchyElement`, and MROs are given a list of the
    corresponding values.

    INPUT:

    - ``value`` -- an object
    - ``succ`` -- a successor function, poset or digraph from which
      one can recover the successors of ``value``
    - ``key`` -- a function taking values as input (default: the
      identity) this function is used to compute comparison keys for
      sorting elements of the hierarchy.

    .. NOTE::

        Constructing a :class:`HierarchyElement` immediately constructs the
        whole hierarchy above it.

    EXAMPLES:

    See the introduction of this module :mod:`sage.misc.c3_controlled`
    for many examples. Here we consider a large example, originally
    taken from the hierarchy of categories above
    :class:`HopfAlgebrasWithBasis`::

        sage: from sage.misc.c3_controlled import HierarchyElement
        sage: G = DiGraph({                                                             # needs sage.graphs
        ....:     44 :  [43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
        ....:     43 :  [42, 41, 40, 36, 35, 39, 38, 37, 33, 32, 31, 30, 29, 28, 27, 26, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
        ....:     42 :  [36, 35, 37, 30, 29, 28, 27, 26, 15, 14, 12, 11, 9, 8, 5, 3, 2, 1, 0],
        ....:     41 :  [40, 36, 35, 33, 32, 31, 30, 29, 28, 27, 26, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
        ....:     40 :  [36, 35, 32, 31, 30, 29, 28, 27, 26, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
        ....:     39 :  [38, 37, 33, 32, 31, 30, 29, 28, 27, 26, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
        ....:     38 :  [37, 33, 32, 31, 30, 29, 28, 27, 26, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
        ....:     37 :  [30, 29, 28, 27, 26, 15, 14, 12, 11, 9, 8, 5, 3, 2, 1, 0],
        ....:     36 :  [35, 30, 29, 28, 27, 26, 15, 14, 12, 11, 9, 8, 5, 3, 2, 1, 0],
        ....:     35 :  [29, 28, 27, 26, 15, 14, 12, 11, 9, 8, 5, 3, 2, 1, 0],
        ....:     34 :  [33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
        ....:     33 :  [32, 31, 30, 29, 28, 27, 26, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
        ....:     32 :  [31, 30, 29, 28, 27, 26, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
        ....:     31 :  [30, 29, 28, 27, 26, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
        ....:     30 :  [29, 28, 27, 26, 15, 14, 12, 11, 9, 8, 5, 3, 2, 1, 0],
        ....:     29 :  [28, 27, 26, 15, 14, 12, 11, 9, 8, 5, 3, 2, 1, 0],
        ....:     28 :  [27, 26, 15, 14, 12, 11, 9, 8, 5, 3, 2, 1, 0],
        ....:     27 :  [15, 14, 12, 11, 9, 8, 5, 3, 2, 1, 0],
        ....:     26 :  [15, 14, 12, 11, 9, 8, 5, 3, 2, 1, 0],
        ....:     25 :  [24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
        ....:     24 :  [4, 2, 1, 0],
        ....:     23 :  [22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
        ....:     22 :  [21, 20, 18, 17, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
        ....:     21 :  [20, 17, 4, 2, 1, 0],
        ....:     20 :  [4, 2, 1, 0],
        ....:     19 :  [18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
        ....:     18 :  [17, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
        ....:     17 :  [4, 2, 1, 0],
        ....:     16 :  [15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
        ....:     15 :  [14, 12, 11, 9, 8, 5, 3, 2, 1, 0],
        ....:     14 :  [11, 3, 2, 1, 0],
        ....:     13 :  [12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
        ....:     12 :  [11, 9, 8, 5, 3, 2, 1, 0],
        ....:     11 :  [3, 2, 1, 0],
        ....:     10 :  [9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
        ....:     9 :  [8, 5, 3, 2, 1, 0],
        ....:     8 :  [3, 2, 1, 0],
        ....:     7 :  [6, 5, 4, 3, 2, 1, 0],
        ....:     6 :  [4, 3, 2, 1, 0],
        ....:     5 :  [3, 2, 1, 0],
        ....:     4 :  [2, 1, 0],
        ....:     3 :  [2, 1, 0],
        ....:     2 :  [1, 0],
        ....:     1 :  [0],
        ....:     0 :  [],
        ....:     })

        sage: # needs sage.combinat sage.graphs
        sage: x = HierarchyElement(44, G)
        sage: x.mro
        [44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0]
        sage: x.cls
        <class '44.cls'>
        sage: x.cls.mro()
        [<class '44.cls'>, <class '43.cls'>, <class '42.cls'>, <class '41.cls'>, <class '40.cls'>, <class '39.cls'>, <class '38.cls'>, <class '37.cls'>, <class '36.cls'>, <class '35.cls'>, <class '34.cls'>, <class '33.cls'>, <class '32.cls'>, <class '31.cls'>, <class '30.cls'>, <class '29.cls'>, <class '28.cls'>, <class '27.cls'>, <class '26.cls'>, <class '25.cls'>, <class '24.cls'>, <class '23.cls'>, <class '22.cls'>, <class '21.cls'>, <class '20.cls'>, <class '19.cls'>, <class '18.cls'>, <class '17.cls'>, <class '16.cls'>, <class '15.cls'>, <class '14.cls'>, <class '13.cls'>, <class '12.cls'>, <class '11.cls'>, <class '10.cls'>, <class '9.cls'>, <class '8.cls'>, <class '7.cls'>, <class '6.cls'>, <class '5.cls'>, <class '4.cls'>, <class '3.cls'>, <class '2.cls'>, <class '1.cls'>, <class '0.cls'>, <... 'object'>]
    """
    @staticmethod
    def __classcall__(cls, value, succ, key=None):
        """
        EXAMPLES::

            sage: # needs sage.combinat sage.graphs
            sage: from sage.misc.c3_controlled import HierarchyElement
            sage: P = Poset((divisors(30), lambda x, y: y.divides(x)), facade=True)
            sage: x = HierarchyElement(10, P)
            sage: x
            10
            sage: x.bases
            [5, 2]
            sage: x.mro
            [10, 5, 2, 1]
        """
        from sage.categories.sets_cat import Sets

        try:
            from sage.combinat.posets.poset_examples import Posets
        except ImportError:
            pass
        else:
            if succ in Posets():
                assert succ in Sets().Facade()
                succ = succ.upper_covers

        try:
            from sage.graphs.digraph import DiGraph
        except ImportError:
            pass
        else:
            if isinstance(succ, DiGraph):
                succ = succ.copy()
                succ._immutable = True
                succ = succ.neighbors_out

        if key is None:
            key = identity

        @cached_function
        def f(x):
            return typecall(cls, x, [f(y) for y in succ(x)], key, f)

        return f(value)

    def __init__(self, value, bases, key, from_value):
        """
        EXAMPLES::

            sage: # needs sage.graphs
            sage: from sage.misc.c3_controlled import HierarchyElement
            sage: P = Poset((divisors(30), lambda x, y: y.divides(x)), facade=True)
            sage: x = HierarchyElement(10, P)
            sage: x
            10
            sage: x.value
            10
            sage: x._bases
            [5, 2]
            sage: x._key
            <built-in function identity>
            sage: x._key(10)
            10

        The ``_from_value`` attribute is a function that can be used
        to reconstruct an element of the hierarchy from its value::

            sage: x._from_value                                                         # needs sage.graphs
            Cached version of <cyfunction HierarchyElement.__classcall__.<locals>.f at ...>
            sage: x._from_value(x.value) is x                                           # needs sage.graphs
            True
        """
        self.value = value
        self._bases = sorted(bases, key=lambda x: key(x.value), reverse=True)
        self._key = key
        self._from_value = from_value

    def __repr__(self):
        """
        Return the representation of ``self`` which is that of its value.

        EXAMPLES::

            sage: from sage.misc.c3_controlled import HierarchyElement
            sage: P = Poset((divisors(30), lambda x, y: y.divides(x)), facade=True)     # needs sage.graphs
            sage: x = HierarchyElement(10, P)                                           # needs sage.graphs
            sage: x                                                                     # needs sage.graphs
            10
        """
        return repr(self.value)

    @lazy_attribute
    def bases(self):
        """
        The bases of ``self``.

        The bases are given as a list of instances of
        :class:`HierarchyElement`, sorted decreasingly according to
        the ``key`` function.

        EXAMPLES::

            sage: # needs sage.combinat sage.graphs
            sage: from sage.misc.c3_controlled import HierarchyElement
            sage: P = Poset((divisors(30), lambda x, y: y.divides(x)), facade=True)
            sage: x = HierarchyElement(10, P)
            sage: x.bases
            [5, 2]
            sage: type(x.bases[0])
            <class 'sage.misc.c3_controlled.HierarchyElement'>
            sage: x.mro
            [10, 5, 2, 1]
            sage: x._bases_controlled
            [5, 2]
        """
        return self._bases

    @lazy_attribute
    def mro(self):
        """
        The MRO for this object, calculated with :meth:`C3_sorted_merge`.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: from sage.misc.c3_controlled import HierarchyElement, C3_sorted_merge, identity
            sage: P = Poset({7: [5, 6], 5: [1, 2], 6: [3, 4]}, facade=True)
            sage: x = HierarchyElement(5, P)
            sage: x.mro
            [5, 2, 1]
            sage: x = HierarchyElement(6, P)
            sage: x.mro
            [6, 4, 3]
            sage: x = HierarchyElement(7, P)
            sage: x.mro
            [7, 6, 5, 4, 3, 2, 1]

            sage: C3_sorted_merge([[6, 4, 3], [5, 2, 1], [6, 5]], identity)             # needs sage.graphs
            ([6, 5, 4, 3, 2, 1], [6, 5, 4])

        TESTS::

            sage: assert all(isinstance(v, Integer) for v in x.mro)                     # needs sage.graphs
        """
        bases = self._bases
        result, suggestion = C3_sorted_merge([base.mro for base in bases]+[[base.value for base in bases]], key=self._key)
        result = [self.value] + result
        self._bases_controlled = suggestion
        return result

    @lazy_attribute
    def _bases_controlled(self):
        """
        A list of bases controlled by :meth:`C3_sorted_merge`.

        This triggers the calculation of the MRO using
        :meth:`C3_sorted_merge`, which sets this attribute as a side
        effect.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: from sage.misc.c3_controlled import HierarchyElement
            sage: P = Poset({7: [5, 6], 5: [1, 2], 6: [3, 4]}, facade=True)
            sage: x = HierarchyElement(7, P)
            sage: x._bases
            [6, 5]
            sage: x._bases_controlled
            [6, 5, 4]
        """
        self.mro
        return self._bases_controlled

    @lazy_attribute
    def mro_standard(self):
        """
        The MRO for this object, calculated with :meth:`C3_merge`.

        EXAMPLES::

            sage: from sage.misc.c3_controlled import HierarchyElement, C3_merge

            sage: # needs sage.graphs
            sage: P = Poset({7: [5, 6], 5: [1, 2], 6: [3, 4]}, facade=True)
            sage: x = HierarchyElement(5, P)
            sage: x.mro_standard
            [5, 2, 1]
            sage: x = HierarchyElement(6, P)
            sage: x.mro_standard
            [6, 4, 3]
            sage: x = HierarchyElement(7, P)
            sage: x.mro_standard
            [7, 6, 4, 3, 5, 2, 1]

            sage: C3_merge([[6, 4, 3], [5, 2, 1], [6, 5]])
            [6, 4, 3, 5, 2, 1]

        TESTS::

            sage: assert all(isinstance(v, Integer) for v in x.mro_standard)            # needs sage.graphs
        """
        bases = self._bases
        return [self.value] + C3_merge([base.mro_standard for base in bases]+[[base.value for base in bases]])

    @lazy_attribute
    def mro_controlled(self):
        """
        The MRO for this object, calculated with :meth:`C3_merge`, under
        control of :func:`C3_sorted_merge`

        EXAMPLES::

            sage: from sage.misc.c3_controlled import HierarchyElement, C3_merge

            sage: # needs sage.graphs
            sage: P = Poset({7: [5, 6], 5: [1, 2], 6: [3, 4]}, facade=True)
            sage: x = HierarchyElement(5, P)
            sage: x.mro_controlled
            [5, 2, 1]
            sage: x = HierarchyElement(6, P)
            sage: x.mro_controlled
            [6, 4, 3]
            sage: x = HierarchyElement(7, P)
            sage: x.mro_controlled
            [7, 6, 5, 4, 3, 2, 1]
            sage: x._bases
            [6, 5]
            sage: x._bases_controlled
            [6, 5, 4]

            sage: C3_merge([[6, 4, 3], [5, 2, 1], [6, 5]])
            [6, 4, 3, 5, 2, 1]
            sage: C3_merge([[6, 4, 3], [5, 2, 1], [6, 5, 4]])
            [6, 5, 4, 3, 2, 1]

        TESTS::

            sage: assert all(isinstance(v, Integer) for v in x.mro_controlled)          # needs sage.graphs
        """
        return [self.value] + C3_merge([base.mro_controlled for base in self._bases]+[self._bases_controlled])

    @cached_method
    def _test_mro(self):
        r"""
        Run consistency tests.

        This checks in particular that the instrumented ``C3`` and
        controlled ``C3`` algorithms give, as desired, the
        decreasingly sorted list of the objects above in the
        hierarchy. For the controlled ``C3`` algorithm, this includes
        both Sage's implementation, and Python's implementation (by
        constructing an appropriate hierarchy of classes).

        It is cached because it is run recursively on the elements
        above ``self``.

        EXAMPLES::

            sage: from sage.misc.c3_controlled import HierarchyElement
            sage: P = Poset({7: [5, 6], 5: [1, 2], 6: [3, 4]}, facade=True)             # needs sage.graphs
            sage: x = HierarchyElement(7, P)                                            # needs sage.graphs
            sage: x._test_mro()                                                         # needs sage.graphs
        """
        for b in self._bases:
            b._test_mro()
        try:
            assert self.mro_standard[0] == self.value
        except ValueError:
            # standard C3 failed to compute a mro; that's ok
            pass
        assert self.mro[0] == self.value
        assert self.mro_controlled[0] == self.value
        assert sorted([x.value for x in self.all_bases()], key=self._key, reverse = True) == self.mro
        assert self.mro == self.mro_controlled
        assert self.cls.mro() == [self._from_value(b).cls for b in self.mro]+[object]

    @lazy_attribute
    def cls(self):
        """
        Return a Python class with inheritance graph parallel to the hierarchy above ``self``.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: from sage.misc.c3_controlled import HierarchyElement
            sage: P = Poset((divisors(30), lambda x, y: y.divides(x)), facade=True)
            sage: x = HierarchyElement(1, P)
            sage: x.cls
            <class '1.cls'>
            sage: x.cls.mro()
            [<class '1.cls'>, <... 'object'>]
            sage: x = HierarchyElement(30, P)
            sage: x.cls
            <class '30.cls'>
            sage: x.cls.mro()
            [<class '30.cls'>, <class '15.cls'>, <class '10.cls'>, <class '6.cls'>, <class '5.cls'>, <class '3.cls'>, <class '2.cls'>, <class '1.cls'>, <... 'object'>]
        """
        super_classes = tuple(self._from_value(base).cls for base in self._bases_controlled)
        if not super_classes:
            super_classes = (object,)
        return dynamic_class("%s.cls"%self, super_classes)

    @cached_method
    def all_bases(self):
        """
        Return the set of all instances of :class:`HierarchyElement` above
        ``self``, ``self`` included.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: from sage.misc.c3_controlled import HierarchyElement
            sage: P = Poset((divisors(30), lambda x, y: y.divides(x)), facade=True)
            sage: HierarchyElement(1, P).all_bases()
            {1}
            sage: HierarchyElement(10, P).all_bases()  # random output
            {10, 5, 2, 1}
            sage: sorted([x.value for x in HierarchyElement(10, P).all_bases()])
            [1, 2, 5, 10]
        """
        return {self} | { x for base in self._bases for x in base.all_bases()  }

    def all_bases_len(self):
        """
        Return the cumulated size of the bases of the elements above ``self`` in the hierarchy.

        EXAMPLES::

            sage: from sage.misc.c3_controlled import HierarchyElement
            sage: P = Poset((divisors(30), lambda x, y: y.divides(x)), facade=True)     # needs sage.graphs
            sage: HierarchyElement(30, P).all_bases_len()                               # needs sage.graphs
            12
        """
        return sum( len(x._bases) for x in self.all_bases())

    def all_bases_controlled_len(self):
        """
        Return the cumulated size of the controlled bases of the elements above ``self`` in the hierarchy.

        EXAMPLES::

            sage: from sage.misc.c3_controlled import HierarchyElement
            sage: P = Poset((divisors(30), lambda x, y: y.divides(x)), facade=True)     # needs sage.graphs
            sage: HierarchyElement(30, P).all_bases_controlled_len()                    # needs sage.graphs
            13
        """
        return sum( len(x._bases_controlled) for x in self.all_bases())
