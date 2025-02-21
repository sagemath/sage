# sage_setup: distribution = sagemath-objects
"""
Fast functions for the category framework


AUTHOR:

- Simon King (initial version)
"""

# ****************************************************************************
#       Copyright (C) 2014 Simon King <simon.king@uni-jena.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

#######################################
#  Sorting

cpdef inline tuple category_sort_key(object category):
    """
    Return ``category._cmp_key``.

    This helper function is used for sorting lists of categories.

    It is semantically equivalent to
    :func:`operator.attrgetter` ``("_cmp_key")``, but currently faster.

    EXAMPLES::

        sage: from sage.categories.category_cy_helper import category_sort_key
        sage: category_sort_key(Rings()) is Rings()._cmp_key
        True
    """
    return category._cmp_key

cpdef tuple _sort_uniq(categories):
    """
    Return the categories after sorting them and removing redundant categories.

    Redundant categories include duplicates and categories which
    are super categories of other categories in the input.

    INPUT:

    - ``categories`` -- list (or iterable) of categories

    OUTPUT: a sorted tuple of mutually incomparable categories

    EXAMPLES::

        sage: Category._sort_uniq([Rings(), Monoids(), Coalgebras(QQ)])
        (Category of rings, Category of coalgebras over Rational Field)

    Note that, in the above example, ``Monoids()`` does not appear
    in the result because it is a super category of ``Rings()``.
    """
    cdef tuple cats = tuple(sorted(categories, key=category_sort_key, reverse=True))
    cdef list result = []
    cdef bint append
    for category in cats:
        append = True
        for cat in result:
            if cat.is_subcategory(category):
                append = False
                break
        if append:
            result.append(category)
    return tuple(result)

cpdef tuple _flatten_categories(categories, ClasscallMetaclass JoinCategory):
    """
    Return the tuple of categories in ``categories``, while
    flattening join categories.

    INPUT:

    - ``categories`` -- list (or iterable) of categories

    - ``JoinCategory`` -- a type such that instances of that type will be
      replaced by its super categories. Usually, this type is
      :class:`JoinCategory`.

    .. NOTE::

        It is needed to provide :class:`~sage.categories.category.JoinCategory` as
        an argument, since we need to prevent a circular import.

    EXAMPLES::

        sage: Category._flatten_categories([Algebras(QQ), Category.join([Monoids(), Coalgebras(QQ)]), Sets()], sage.categories.category.JoinCategory)
        (Category of algebras over Rational Field, Category of monoids,
         Category of coalgebras over Rational Field, Category of sets)
    """
    # Invariant: the super categories of a JoinCategory are not JoinCategories themselves
    cdef list out = []
    for category in categories:
        if isinstance(category, JoinCategory):
            out.extend(category.super_categories())
        else:
            out.append(category)
    return tuple(out)

#############################################
#  Join

cdef bint is_supercategory_of_done(new_cat, dict done) noexcept:
    # This is a helper function. It replaces the closure
    # any(cat.is_subcategory(new_cat) for cat in done)
    for cat in done:
        if cat.is_subcategory(new_cat):
            return True
    return False

cpdef tuple join_as_tuple(tuple categories, tuple axioms, tuple ignore_axioms):
    """
    Helper for :meth:`~sage.categories.category.Category.join`.

    INPUT:

    - ``categories`` -- tuple of categories to be joined
    - ``axioms`` -- tuple of strings; the names of some
      supplementary axioms
    - ``ignore_axioms`` -- tuple of pairs ``(cat, axiom)``, such
      that ``axiom`` will not be applied to ``cat``, should ``cat``
      occur in the algorithm

    EXAMPLES::

        sage: from sage.categories.category_cy_helper import join_as_tuple
        sage: T = (Coalgebras(QQ), Sets().Finite(), Algebras(ZZ), SimplicialComplexes())
        sage: join_as_tuple(T,(),())
        (Category of algebras over Integer Ring,
         Category of finite monoids,
         Category of finite additive groups,
         Category of coalgebras over Rational Field,
         Category of finite simplicial complexes)
        sage: join_as_tuple(T,('WithBasis',),())
        (Category of algebras with basis over Integer Ring,
         Category of finite monoids,
         Category of coalgebras with basis over Rational Field,
         Category of finite additive groups,
         Category of finite simplicial complexes)
        sage: join_as_tuple(T,(),((Monoids(),'Finite'),))
        (Category of algebras over Integer Ring,
         Category of finite additive groups,
         Category of coalgebras over Rational Field,
         Category of finite simplicial complexes)
    """
    cdef set axiomsS = set(axioms)
    for category in categories:
        axiomsS.update(category.axioms())
    cdef dict done = dict()
    cdef set todo = set()
    cdef frozenset axs
    for category in categories:
        axs = category.axioms()
        for (cat, axiom) in ignore_axioms:
            if category.is_subcategory(cat):
                axs = axs | {axiom}
        done[category] = axs
        for axiom in axiomsS.difference(axs):
            todo.add((category, axiom))

    # Invariants:
    # - the current list of categories is stored in the keys of ``done``
    # - todo contains the ``complement`` of done; i.e.
    #   for category in the keys of done,
    #   (category, axiom) is in todo iff axiom is not in done[category]
    cdef list new_cats
    cdef set new_axioms
    while todo:
        (category, axiom) = todo.pop()
        # It's easier to remove categories from done than from todo
        # So we check that ``category`` had not been removed
        if category not in done:
            continue

        # Removes redundant categories
        new_cats = [new_cat for new_cat in <tuple>(category._with_axiom_as_tuple(axiom))
                    if not is_supercategory_of_done(new_cat, done)]
        for cat in list(done.keys()):
            for new_cat in new_cats:
                if new_cat.is_subcategory(cat):
                    del done[cat]
                    break

        new_axioms = set()
        for new_cat in new_cats:
            for axiom in new_cat.axioms():
                if axiom not in axiomsS:
                    new_axioms.add(axiom)

        # Mark old categories with new axioms as todo
        for category in done:
            for axiom in new_axioms:
                todo.add((category, axiom))
        for new_cat in new_cats:
            axs = new_cat.axioms()
            for (cat, axiom) in ignore_axioms:
                if new_cat.is_subcategory(cat):
                    axs = axs | {axiom}
            done[new_cat] = axs
            for axiom in axiomsS.difference(axs):
                todo.add((new_cat, axiom))

    return _sort_uniq(done)


#############################################
#  Axiom related functions

cdef class AxiomContainer(dict):
    """
    A fast container for axioms.

    This is derived from :class:`dict`. A key is the name of an axiom. The
    corresponding value is the "rank" of this axiom, that is used to order the
    axioms in :func:`canonicalize_axioms`.

    EXAMPLES::

        sage: all_axioms = sage.categories.category_with_axiom.all_axioms
        sage: isinstance(all_axioms, sage.categories.category_with_axiom.AxiomContainer)
        True
    """
    def add(self, axiom):
        """
        Add a new axiom name, of the next rank.

        EXAMPLES::

            sage: all_axioms = sage.categories.category_with_axiom.all_axioms
            sage: m = max(all_axioms.values())
            sage: all_axioms.add('Awesome')
            sage: all_axioms['Awesome'] == m + 1
            True

        To avoid side effects, we remove the added axiom::

            sage: del all_axioms['Awesome']
        """
        self[axiom] = len(self)

    def __iadd__(self, L):
        """
        Inline addition, which means to add a list of axioms to the container.

        EXAMPLES::

            sage: all_axioms = sage.categories.category_with_axiom.all_axioms
            sage: m = max(all_axioms.values())
            sage: all_axioms += ('Fancy', 'Awesome')
            sage: all_axioms['Awesome'] == m + 2
            True

        To avoid side effects, we delete the axioms that we just added::

            sage: del all_axioms['Awesome'], all_axioms['Fancy']
        """
        for axiom in L:
            self.add(axiom)
        return self


cpdef inline get_axiom_index(AxiomContainer all_axioms, str axiom):
    """
    Helper function: Return the rank of an axiom.

    INPUT:

    - ``all_axioms`` -- the axiom collection
    - ``axiom`` -- string, name of an axiom

    EXAMPLES::

        sage: all_axioms = sage.categories.category_with_axiom.all_axioms
        sage: from sage.categories.category_cy_helper import get_axiom_index
        sage: get_axiom_index(all_axioms, 'AdditiveCommutative') == all_axioms['AdditiveCommutative']
        True
    """
    return (<dict>all_axioms)[axiom]


cpdef tuple canonicalize_axioms(AxiomContainer all_axioms, axioms):
    r"""
    Canonicalize a set of axioms.

    INPUT:

    - ``all_axioms`` -- all available axioms

    - ``axioms`` -- set (or iterable) of axioms

    .. NOTE::

        :class:`AxiomContainer` provides a fast container for axioms, and the
        collection of axioms is stored in
        :mod:`sage.categories.category_with_axiom`. In order to avoid circular
        imports, we expect that the collection of all axioms is provided as an
        argument to this auxiliary function.

    OUTPUT:

    A set of axioms as a tuple sorted according to the order of the
    tuple ``all_axioms`` in :mod:`sage.categories.category_with_axiom`.

    EXAMPLES::

        sage: from sage.categories.category_with_axiom import canonicalize_axioms, all_axioms
        sage: canonicalize_axioms(all_axioms, ["Commutative", "Connected", "WithBasis", "Finite"])
        ('Finite', 'Connected', 'WithBasis', 'Commutative')
        sage: canonicalize_axioms(all_axioms, ["Commutative", "Connected", "Commutative", "WithBasis", "Finite"])
        ('Finite', 'Connected', 'WithBasis', 'Commutative')
    """
    cdef list L = list(set(axioms))
    L.sort(key=(all_axioms).__getitem__)
    return tuple(L)
