"""
Constructor for permutations

This module contains the generic constructor to build element of the symmetric
groups (or more general permutation groups) called
:class:`~sage.groups.perm_gps.permgroup_element.PermutationGroupElement`. These
objects have a more group theoretic flavor than the more combinatorial
:class:`~sage.combinat.permutation.Permutation`.
"""
# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2006 David Joyner
#       Copyright (C) 2019 Vincent Delecroix <20100.delecroix@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from . import permgroup_element
from sage.misc.sage_eval import sage_eval
from sage.misc.lazy_import import lazy_import
from sage.interfaces.gap import GapElement
from sage.libs.gap.element import GapElement_Permutation

lazy_import('sage.combinat.permutation', ['Permutation', 'from_cycles'])

try:
    from sage.libs.pari.all import pari_gen
except ImportError:
    pari_gen = ()


def PermutationGroupElement(g, parent=None, check=True):
    r"""
    Build a permutation from ``g``.

    INPUT:

    - ``g`` -- either

      - a list of images

      - a tuple describing a single cycle

      - a list of tuples describing the cycle decomposition

      - a string describing the cycle decomposition

    - ``parent`` -- (optional) an ambient permutation group for the result;
      it is mandatory if you want a permutation on a domain different
      from `\{1, \ldots, n\}`

    - ``check`` -- boolean (default: ``True``); whether additional check are
      performed. Setting it to ``False`` is likely to result in faster code.

    EXAMPLES:

    Initialization as a list of images::

        sage: p = PermutationGroupElement([1,4,2,3])
        sage: p
        (2,4,3)
        sage: p.parent()
        Symmetric group of order 4! as a permutation group

    Initialization as a list of cycles::

        sage: p = PermutationGroupElement([(3,5),(4,6,9)])
        sage: p
        (3,5)(4,6,9)
        sage: p.parent()
        Symmetric group of order 9! as a permutation group

    Initialization as a string representing a cycle decomposition::

        sage: p = PermutationGroupElement('(2,4)(3,5)')
        sage: p
        (2,4)(3,5)
        sage: p.parent()
        Symmetric group of order 5! as a permutation group

    By default the constructor assumes that the domain is `\{1, \dots, n\}`
    but it can be set to anything via its second ``parent`` argument::

        sage: S = SymmetricGroup(['a', 'b', 'c', 'd', 'e'])
        sage: PermutationGroupElement(['e', 'c', 'b', 'a', 'd'], S)
        ('a','e','d')('b','c')
        sage: PermutationGroupElement(('a', 'b', 'c'), S)
        ('a','b','c')
        sage: PermutationGroupElement([('a', 'c'), ('b', 'e')], S)
        ('a','c')('b','e')
        sage: PermutationGroupElement("('a','b','e')('c','d')", S)
        ('a','b','e')('c','d')

    But in this situation, you might want to use the more direct::

        sage: S(['e', 'c', 'b', 'a', 'd'])
        ('a','e','d')('b','c')
        sage: S(('a', 'b', 'c'))
        ('a','b','c')
        sage: S([('a', 'c'), ('b', 'e')])
        ('a','c')('b','e')
        sage: S("('a','b','e')('c','d')")
        ('a','b','e')('c','d')
    """
    if isinstance(g, permgroup_element.PermutationGroupElement):
        if parent is None or g.parent() is parent:
            return g

    if parent is None:
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup

        try:
            v = standardize_generator(g, None)
        except KeyError:
            raise ValueError("invalid permutation vector: %s" % g)
        parent = SymmetricGroup(len(v))
        # We have constructed the parent from the element and already checked
        #   that it is a valid permutation
        check = False

    return parent.element_class(g, parent, check)


def string_to_tuples(g):
    """
    EXAMPLES::

        sage: from sage.groups.perm_gps.constructor import string_to_tuples
        sage: string_to_tuples('(1,2,3)')
        [(1, 2, 3)]
        sage: string_to_tuples('(1,2,3)(4,5)')
        [(1, 2, 3), (4, 5)]
        sage: string_to_tuples(' (1,2, 3) (4,5)')
        [(1, 2, 3), (4, 5)]
        sage: string_to_tuples('(1,2)(3)')
        [(1, 2), (3,)]
    """
    if not isinstance(g, str):
        raise ValueError("g (= %s) must be a string" % g)
    elif g == '()':
        return []
    g = g.replace('\n', '').replace(' ', '').replace(')(', '),(').replace(')', ',)')
    g = '[' + g + ']'
    return sage_eval(g, preparse=False)


def standardize_generator(g, convert_dict=None, as_cycles=False):
    r"""
    Standardize the input for permutation group elements to a list
    or a list of tuples.

    This was factored out of the
    ``PermutationGroupElement.__init__`` since
    ``PermutationGroup_generic.__init__`` needs to do the same computation
    in order to compute the domain of a group when it's not explicitly
    specified.

    INPUT:

    - ``g`` -- a :class:`list`, :class:`tuple`, :class:`string`, :class:`GapElement`,
      :class:`PermutationGroupElement`, or :class:`Permutation`

    - ``convert_dict`` -- (optional) a dictionary used to convert the
      points to a number compatible with GAP

    - ``as_cycles`` -- boolean (default: ``False``); whether the output should be
      as cycles or in one-line notation

    OUTPUT: the permutation in as a list in one-line notation or a list of
    cycles as tuples

    EXAMPLES::

        sage: from sage.groups.perm_gps.constructor import standardize_generator
        sage: standardize_generator('(1,2)')
        [2, 1]

        sage: p = PermutationGroupElement([(1,2)])
        sage: standardize_generator(p)
        [2, 1]
        sage: standardize_generator(p._gap_())
        [2, 1]
        sage: standardize_generator((1,2))
        [2, 1]
        sage: standardize_generator([(1,2)])
        [2, 1]

        sage: standardize_generator(p, as_cycles=True)
        [(1, 2)]
        sage: standardize_generator(p._gap_(), as_cycles=True)
        [(1, 2)]
        sage: standardize_generator((1,2), as_cycles=True)
        [(1, 2)]
        sage: standardize_generator([(1,2)], as_cycles=True)
        [(1, 2)]

        sage: standardize_generator(Permutation([2,1,3]))
        [2, 1, 3]
        sage: standardize_generator(Permutation([2,1,3]), as_cycles=True)
        [(1, 2), (3,)]

    ::

        sage: d = {'a': 1, 'b': 2}
        sage: p = SymmetricGroup(['a', 'b']).gen(0); p
        ('a','b')
        sage: standardize_generator(p, convert_dict=d)
        [2, 1]
        sage: standardize_generator(p._gap_(), convert_dict=d)
        [2, 1]
        sage: standardize_generator(('a','b'), convert_dict=d)
        [2, 1]
        sage: standardize_generator([('a','b')], convert_dict=d)
        [2, 1]

        sage: standardize_generator(p, convert_dict=d, as_cycles=True)
        [(1, 2)]
        sage: standardize_generator(p._gap_(), convert_dict=d, as_cycles=True)
        [(1, 2)]
        sage: standardize_generator(('a','b'), convert_dict=d, as_cycles=True)
        [(1, 2)]
        sage: standardize_generator([('a','b')], convert_dict=d, as_cycles=True)
        [(1, 2)]
    """
    if isinstance(g, pari_gen):
        g = list(g)

    needs_conversion = True

    if isinstance(g, GapElement_Permutation):
        g = g.sage()
        needs_conversion = False
    if isinstance(g, GapElement):
        g = str(g)
        needs_conversion = False
    if isinstance(g, Permutation):
        if as_cycles:
            return g.cycle_tuples()
        return g._list
    elif isinstance(g, permgroup_element.PermutationGroupElement):
        if not as_cycles:
            l = list(range(1, g.parent().degree() + 1))
            return g._act_on_list_on_position(l)
        g = g.cycle_tuples()
    elif isinstance(g, str):
        g = string_to_tuples(g)
    elif isinstance(g, tuple) and (len(g) == 0 or not isinstance(g[0], tuple)):
        g = [g]

    # Get the permutation in list notation
    if isinstance(g, list) and not (g and isinstance(g[0], tuple)):
        if convert_dict is not None and needs_conversion:
            for i, x in enumerate(g):
                g[i] = convert_dict[x]
        if as_cycles:
            return Permutation(g).cycle_tuples()
        return g

    # Otherwise it is in cycle notation
    if convert_dict is not None and needs_conversion:
        g = [tuple([convert_dict[x] for x in cycle]) for cycle in g]
    if not as_cycles:
        degree = max([1] + [max(cycle + (1,)) for cycle in g])
        g = from_cycles(degree, g)
    return g
