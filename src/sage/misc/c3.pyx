"""
The C3 algorithm

The C3 algorithm is used as method resolution order for new style
classes in Python. The implementation here is used to order the list
of super categories of a category.

AUTHOR:

- Simon King (2011-11): initial version.
"""

# ****************************************************************************
#       Copyright (C) 2011 Simon King <simon.king@uni-jena.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************


cpdef list C3_algorithm(object start, str bases, str attribute, bint proper):
    """
    An implementation of the C3 algorithm.

    C3 is the algorithm used by Python to construct the method
    resolution order for new style classes involving multiple
    inheritance.

    After :issue:`11943` this implementation was used to compute the
    list of super categories of a category; see
    :meth:`~sage.categories.category.Category.all_super_categories`.
    The purpose is to ensure that list of super categories matches
    with the method resolution order of the parent or element classes
    of a category.

    Since :issue:`13589`, this implementation is superseded by that in
    :mod:`sage.misc.c3_controlled`, that puts the ``C3`` algorithm
    under control of some total order on categories.  This guarantees
    that ``C3`` always finds a consistent Method Resolution Order. For
    background, see :mod:`sage.misc.c3_controlled`.

    INPUT:

    - ``start`` -- an object; the returned list is built upon data
      provided by certain attributes of ``start``
    - ``bases`` -- string; the name of an attribute of ``start``
      providing a list of objects
    - ``attribute`` -- string; the name of an attribute of the
      objects provided in ``getattr(start,bases)``. That attribute is
      supposed to provide a list.

    ASSUMPTIONS:

    Our implementation of the algorithm only works on lists of
    objects that compare equal if and only if they are identical.

    OUTPUT:

    A list, the result of the C3 algorithm applied to the list
    ``[getattr(X,attribute) for X in getattr(start,bases)]``.

    EXAMPLES:

    We create a class for elements in a hierarchy that uses the ``C3``
    algorithm to compute, for each element, a linear extension of the
    elements above it::

    .. TODO:: Move back the __init__ at the beginning

        sage: from sage.misc.c3 import C3_algorithm
        sage: class HierarchyElement(UniqueRepresentation):
        ....:     @lazy_attribute
        ....:     def _all_bases(self):
        ....:         return C3_algorithm(self, '_bases', '_all_bases', False)
        ....:     def __repr__(self):
        ....:         return self._name
        ....:     def __init__(self, name, bases):
        ....:         self._name = name
        ....:         self._bases = list(bases)

    We construct a little hierarchy::

        sage: T = HierarchyElement("T", ())
        sage: X = HierarchyElement("X", (T,))
        sage: Y = HierarchyElement("Y", (T,))
        sage: A = HierarchyElement("A", (X, Y))
        sage: B = HierarchyElement("B", (Y, X))
        sage: Foo = HierarchyElement("Foo", (A, B))

    And inspect the linear extensions associated to each element::

        sage: T._all_bases
        [T]
        sage: X._all_bases
        [X, T]
        sage: Y._all_bases
        [Y, T]
        sage: A._all_bases
        [A, X, Y, T]
        sage: B._all_bases
        [B, Y, X, T]

    So far so good. However::

        sage: Foo._all_bases
        Traceback (most recent call last):
        ...
        ValueError: Cannot merge the items X, Y.

    The ``C3`` algorithm is not able to create a consistent linear
    extension. Indeed, its specifications impose that, if ``X`` and
    ``Y`` appear in a certain order in the linear extension for an
    element of the hierarchy, then they should appear in the same
    order for any lower element. This is clearly not possibly for
    ``Foo``, since ``A`` and ``B`` impose incompatible orders. If the
    above was a hierarchy of classes, Python would complain that it
    cannot calculate a consistent Method Resolution Order.

    TESTS:

    Regression test for bug #1 of :issue:`13501`::

        sage: class C(): pass
        sage: class F(): pass
        sage: class G(): pass
        sage: class B(C,F):    pass
        sage: class D(F,G):    pass
        sage: class E(F):      pass
        sage: class A(B,D,E):  pass
        sage: [cls.__name__ for cls in A.mro()]
        ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'object']

        sage: C = HierarchyElement("C", ())
        sage: F = HierarchyElement("F", ())
        sage: G = HierarchyElement("G", ())
        sage: B = HierarchyElement("B", (C, F))
        sage: D = HierarchyElement("D", (F, G))
        sage: E = HierarchyElement("E", (F,))
        sage: A = HierarchyElement("A", (B, D, E))
        sage: A._all_bases
        [A, B, C, D, E, F, G]

    Regression test for bug #2 of :issue:`13501`. The following should
    fail since ``A`` asks for ``B`` to come before ``C``, where as
    ``B`` is a super class of ``C``::

        sage: class B(): pass
        sage: class C(B): pass
        sage: class A(B, C): pass
        Traceback (most recent call last):
        ...
        TypeError: ...Cannot create a consistent method resolution order (MRO)
        for bases ...

        sage: B = HierarchyElement("B", ())
        sage: C = HierarchyElement("C", (B,))
        sage: A = HierarchyElement("A", (B,C))
        sage: A._all_bases
        Traceback (most recent call last):
        ...
        ValueError: Cannot merge the items B, C, B.

    Since :issue:`11943`, the following consistency tests are part
    of the test suites of categories (except for hom categories)::

        sage: C = Category.join([HopfAlgebrasWithBasis(QQ), FiniteEnumeratedSets()])
        sage: C.parent_class.mro() == [x.parent_class for x in C.all_super_categories()]+[object]
        True
        sage: C.element_class.mro() == [x.element_class for x in C.all_super_categories()]+[object]
        True
    """
    cdef list out
    if proper:
        out = []
    else:
        out = [start]
    cdef list args = getattr(start, bases)
    if not args:
        return out
    # Data structure / invariants:
    # We will be working with the MRO's of the super objects
    # together with the list of bases of ``self``.
    # Each list is split between its head (in ``heads``) and tail (in
    # ``tails'') . Each tail is stored reversed, so that we can use a
    # cheap pop() in lieue of pop(0). A duplicate of the tail is
    # stored as a set in ``tailsets`` for cheap membership testing.
    # Since we actually want comparison by identity, not equality,
    # what we store is the set of memory locations of objects
    cdef object O, X
    cdef list tails = [getattr(obj, attribute) for obj in args]
    tails.append(args)
    tails              = [list(reversed(tail))                   for tail in tails]
    cdef list heads    = [tail.pop()                             for tail in tails]
    cdef list tailsets = [set([<size_t><void *>O for O in tail]) for tail in tails]

    cdef int i, j, nbheads
    nbheads = len(heads)
    cdef bint next_item_found
    cdef list tail_list

    while nbheads:
        for i from 0 <= i < nbheads:
            O = heads[i]
            # Does O appear in none of the tails?  ``all(O not in tail for tail in tailsets)``
            next_item_found = True
            for j from 0 <= j < nbheads:
                if j == i:
                    continue
                if <size_t><void *>O in <set>tailsets[j]:
                    next_item_found = False
                    break
            if next_item_found:
                out.append(O)
                # Clear O from other heads, removing the line altogether
                # if the tail is already empty.
                # j goes down so that ``del heads[j]`` does not screw up the numbering
                for j from nbheads > j >= 0:
                    if heads[j] is O:
                        tail_list = tails[j]
                        if tail_list:
                            X = tail_list.pop()
                            heads[j] = X
                            <set>tailsets[j].remove(<size_t><void *>X)
                        else:
                            del heads[j]
                            del tails[j]
                            del tailsets[j]
                            nbheads -= 1
                break
        if not next_item_found:
            # No head is available
            raise ValueError("Cannot merge the items {}.".format(', '.join([repr(head) for head in heads])))
    return out
