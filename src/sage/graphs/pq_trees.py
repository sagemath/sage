# sage_setup: distribution = sagemath-graphs
r"""
PQ-Trees

This module implements PQ-Trees, a data structure use to represent all
permutations of the columns of a matrix which satisfy the *consecutive ones*
*property*:

  A binary matrix satisfies the *consecutive ones property* if the 1s are
  contiguous in each of its rows (or equivalently, if no row contains the regexp
  pattern `10^+1`).

  Alternatively, one can say that a sequence of sets `S_1,...,S_n` satisfies the
  *consecutive ones property* if for any `x` the indices of the sets containing
  `x` is an interval of `[1,n]`.

This module is used for the recognition of Interval Graphs (see
:meth:`~sage.graphs.generic_graph.GenericGraph.is_interval`).

**P-tree and Q-tree**


- A `P`-tree with children `c_1,...,c_k` (which can be `P`-trees, `Q`-trees, or
  actual sets of points) indicates that all `k!` permutations of the children
  are allowed.

  Example: `\{1,2\},\{3,4\},\{5,6\}` (disjoint sets can be permuted in any way)

- A `Q`-tree with children `c_1,...,c_k` (which can be `P`-trees, `Q`-trees, or
  actual sets of points) indicates that only two permutations of its children
  are allowed: `c_1,...,c_k` or `c_k,...,c_1`.

  Example: `\{1,2\},\{2,3\},\{3,4\},\{4,5\},\{5,6\}` (only two permutations of
  these sets have the *consecutive ones property*).

**Computation of all possible orderings**

#. In order to compute all permutations of a sequence of sets `S_1,...,S_k`
   satisfying the *consecutive ones property*, we initialize `T` as a `P`-tree
   whose children are all the `S_1,...,S_k`, thus representing the set of all
   `k!` permutations of them.

#. We select some element `x` and update the data structure `T` to restrict the
   permutations it describes to those that keep the occurrences of `x` on an
   interval of `[1,...,k]`. This will result in a new `P`-tree whose children
   are:

   * all `\bar c_x` sets `S_i` which do *not* contain `x`.
   * a new `P`-tree whose children are the `c_x` sets `S_i` containing `x`.

   This describes the set of all `c_x!\times \bar c'_x!` permutations of
   `S_1,...,S_k` that keep the sets containing `x` on an interval.

#. We take a second element `x'` and update the data structure `T` to restrict
   the permutations it describes to those that keep `x'` on an interval of
   `[1,...,k]`. The sets `S_1,...,S_k` belong to 4 categories:

   * The family `S_{00}` of sets which do not contain any of
     `x,x'`.

   * The family `S_{01}` of sets which contain `x'` but do not contain
     `x`.

   * The family `S_{10}` of sets which contain `x` but do not contain
     `x'`.

   * The family `S_{11}` of sets which contain `x'` and `x'`.

   With these notations, the permutations of `S_1,...,S_k` which keep the
   occurrences of `x` and `x'` on an interval are of two forms:

   * <some sets `S_{00}`>, <sets from `S_{10}`>, <sets from `S_{11}`>, <sets from `S_{01}`>, <other sets from `S_{00}`>
   * <some sets `S_{00}`>, <sets from `S_{01}`>, <sets from `S_{11}`>, <sets from `S_{10}`>, <other sets from `S_{00}`>

   These permutations can be modeled with the following `PQ`-tree:

   * A `P`-tree whose children are:

     * All sets from `S_{00}`
     * A `Q`-tree whose children are:

       * A `P`-tree with whose children are the sets from `S_{10}`
       * A `P`-tree with whose children are the sets from `S_{11}`
       * A `P`-tree with whose children are the sets from `S_{01}`

#. One at a time, we update the data structure with each element until they are
   all exhausted, or until we reach a proof that no permutation satisfying the
   *consecutive ones property* exists.

   Using these two types of tree, and exploring the different cases of
   intersection, it is possible to represent all the possible permutations of
   our sets satisfying our constraints, or to prove that no such ordering
   exists. This is the whole purpose of this module, and is explained with more
   details in many places, for example in the following document from Hajiaghayi
   [Haj2000]_.

Authors:

Nathann Cohen (initial implementation)


Methods and functions
---------------------
"""

# ****************************************************************************
#       Copyright (C) 2012 Nathann Cohen <nathann.cohen@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

# Constants, to make the code more readable

FULL = 2
PARTIAL = 1
EMPTY = 0
ALIGNED = True
UNALIGNED = False


##########################################################################
# Some Helper Functions                                                  #
#                                                                        #
# As the elements of a PQ-Tree can be either P-Trees, Q-Trees, or the    #
# sets themselves (the leaves), the following functions are              #
# meant to be applied both on PQ-Trees and Sets, and mimic for the       #
# latter the behaviour we expect from the corresponding methods          #
# defined in class PQ                                                    #
##########################################################################

def _set_contiguous(tree, x):
    """
    Helper function for updating ``tree``.

    The objective is to ensure that the sets containing ``x`` are contiguous for
    any admissible permutation of its subtrees.

    TESTS::

        sage: from sage.graphs.pq_trees import _set_contiguous, P
        sage: p1 = P([[0, 1], [1, 2], [2, 3], [3, 0]])
        sage: _set_contiguous(p1, 0)
        (1, True)
        sage: p1
        ('P', [{1, 2}, {2, 3}, ('P', [{0, 1}, {0, 3}])])
        sage: _set_contiguous(p1, 2)
        (1, True)
        sage: p1
        ('P', [('P', [{0, 1}, {0, 3}]), ('P', [{1, 2}, {2, 3}])])
        sage: _set_contiguous(p1, 1)
        (1, False)
        sage: p1
        ('P', [('Q', [{0, 3}, {0, 1}, {1, 2}, {2, 3}])])
        sage: p2 = P([[0, 1], [0, 2], [0, 3]])
        sage: _set_contiguous(p2, 0)
        (2, True)
        sage: p2
        ('P', [{0, 1}, {0, 2}, {0, 3}])
        sage: _set_contiguous(p2, Set([1, 2]))
        (0, True)
    """
    if isinstance(tree, PQ):
        return tree.set_contiguous(x)
    elif x in tree:
        return (FULL, ALIGNED)
    return (EMPTY, ALIGNED)


def _new_P(liste):
    """
    Helper function returning a new P-tree.

    TESTS::

        sage: from sage.graphs.pq_trees import _new_P
        sage: _new_P([[1,2], [2,3]])
        ('P', [{1, 2}, {2, 3}])
        sage: _new_P([[1,2]])
        [1, 2]
    """
    if len(liste) > 1:
        return P(liste)
    return liste[0]


def _new_Q(liste):
    """
    Helper function returning a new Q-tree.

    TESTS::

        sage: from sage.graphs.pq_trees import _new_Q
        sage: _new_Q([[1,2], [2,3]])
        ('Q', [{1, 2}, {2, 3}])
        sage: _new_Q([[1,2]])
        [1, 2]
    """
    if len(liste) > 1:
        return Q(liste)
    return liste[0]


def _flatten(x):
    """
    Helper function returning a flatten version of ``x``, if ``x`` is a PQ-tree.

    TESTS::

        sage: from sage.graphs.pq_trees import P, Q, _flatten
        sage: p = Q([P([[1,2], [2,3]])])
        sage: _flatten(p)
        ('P', [{1, 2}, {2, 3}])
        sage: _flatten([p, p])
        [('Q', [('P', [{1, 2}, {2, 3}])]), ('Q', [('P', [{1, 2}, {2, 3}])])]
    """
    if isinstance(x, PQ):
        return x.flatten()
    return x


impossible_msg = "Impossible"


def reorder_sets(sets):
    r"""
    Reorders a collection of sets such that each element appears on an
    interval.

    Given a collection of sets `C = S_1,...,S_k` on a ground set `X`,
    this function attempts to reorder them in such a way that `\forall
    x \in X` and `i<j` with `x\in S_i, S_j`, then `x\in S_l` for every
    `i<l<j` if it exists.

    INPUT:

    - ``sets`` -- list of instances of ``list, Set`` or ``set``

    ALGORITHM: PQ-Trees

    EXAMPLES:

    There is only one way (up to reversal) to represent contiguously
    the sequence of sets `\{i-1, i, i+1\}`::

        sage: from sage.graphs.pq_trees import reorder_sets
        sage: seq = [Set([i-1,i,i+1]) for i in range(1,15)]

    We apply a random permutation::

        sage: p = Permutations(len(seq)).random_element()
        sage: seq = [ seq[p(i+1)-1] for i in range(len(seq)) ]
        sage: ordered = reorder_sets(seq)
        sage: if not 0 in ordered[0]:
        ....:    ordered = ordered.reverse()
        sage: print(ordered)
        [{0, 1, 2}, {1, 2, 3}, {2, 3, 4}, {3, 4, 5}, {4, 5, 6}, {5, 6, 7},
         {8, 6, 7}, {8, 9, 7}, {8, 9, 10}, {9, 10, 11}, {10, 11, 12},
         {11, 12, 13}, {12, 13, 14}, {13, 14, 15}]
    """
    if len(sets) <= 2:
        return sets

    s = set().union(*sets)  # union of the sets

    tree = P(sets)

    for i in s:
        tree.set_contiguous(i)
        tree = _flatten(tree)

    return tree.ordering()


class PQ:
    r"""
    PQ-Trees

    This class should not be instantiated by itself: it is extended by
    :class:`P` and :class:`Q`. See the documentation of
    :mod:`sage.graphs.pq_trees` for more information.

    AUTHOR : Nathann Cohen
    """

    def __init__(self, seq):
        r"""
        Construction of a PQ-Tree.

        EXAMPLES::

            sage: from sage.graphs.pq_trees import P, Q
            sage: p = Q([[1,2], [2,3], P([[2,4], [2,8], [2,9]])])

        :issue:`17787`::

            sage: Graph('GvGNp?').is_interval()
            False
        """
        from sage.sets.set import Set

        self._children = []
        for e in seq:
            if isinstance(e, list):
                e = Set(e)

            if e not in self._children:
                self._children.append(e)

    def reverse(self):
        r"""
        Recursively reverse ``self`` and its children.

        EXAMPLES::

            sage: from sage.graphs.pq_trees import P, Q
            sage: p = Q([[1,2], [2,3], P([[2,4], [2,8], [2,9]])])
            sage: p.ordering()
            [{1, 2}, {2, 3}, {2, 4}, {8, 2}, {9, 2}]
            sage: p.reverse()
            sage: p.ordering()
            [{9, 2}, {8, 2}, {2, 4}, {2, 3}, {1, 2}]
        """
        for i in self._children:
            if isinstance(i, PQ):
                i.reverse()

        self._children.reverse()

    def __contains__(self, v):
        r"""
        Test whether there exists an element of ``self`` containing
        an element ``v``.

        INPUT:

        - ``v`` -- an element of the ground set

        EXAMPLES::

            sage: from sage.graphs.pq_trees import P, Q
            sage: p = Q([[1,2], [2,3], P([[2,4], [2,8], [2,9]])])
            sage: 5 in p
            False
            sage: 9 in p
            True
        """
        return any(v in i for i in self)

    def __iter__(self):
        r"""
        Iterate over the children of ``self``.

        EXAMPLES::

            sage: from sage.graphs.pq_trees import P, Q
            sage: p = Q([[1,2], [2,3], P([[2,4], [2,8], [2,9]])])
            sage: for i in p:
            ....:     print(i)
            {1, 2}
            {2, 3}
            ('P', [{2, 4}, {8, 2}, {9, 2}])
        """
        yield from self._children

    def number_of_children(self):
        r"""
        Return the number of children of ``self``.

        EXAMPLES::

            sage: from sage.graphs.pq_trees import P, Q
            sage: p = Q([[1,2], [2,3], P([[2,4], [2,8], [2,9]])])
            sage: p.number_of_children()
            3
        """
        return len(self._children)

    def ordering(self):
        r"""
        Return the current ordering given by listing the leaves from
        left to right.

        EXAMPLES::

            sage: from sage.graphs.pq_trees import P, Q
            sage: p = Q([[1,2], [2,3], P([[2,4], [2,8], [2,9]])])
            sage: p.ordering()
            [{1, 2}, {2, 3}, {2, 4}, {8, 2}, {9, 2}]
        """
        value = []
        for i in self:
            if isinstance(i, PQ):
                value.extend(i.ordering())
            else:
                value.append(i)

        return value

    def __repr__(self) -> str:
        r"""
        Succinctly represent ``self``.

        EXAMPLES::

            sage: from sage.graphs.pq_trees import P, Q
            sage: p = Q([[1,2], [2,3], P([[2,4], [2,8], [2,9]])])
            sage: print(p)
            ('Q', [{1, 2}, {2, 3}, ('P', [{2, 4}, {8, 2}, {9, 2}])])
        """
        return str((("P" if isinstance(self, P) else "Q"), self._children))

    def simplify(self, v, left=False, right=False):
        r"""
        Return a simplified copy of ``self`` according to the element ``v``.

        If ``self`` is a partial P-tree for ``v``, we would like to
        restrict the permutations of its children to permutations
        keeping the children containing ``v`` contiguous. This
        function also "locks" all the elements not containing ``v``
        inside a `P`-tree, which is useful when one want to keep the
        elements containing ``v`` on one side (which is the case when
        this method is called).

        INPUT:

        - ``left``, ``right`` -- booleans; whether ``v`` is aligned to the
          right or to the left

        - ``v`` -- an element of the ground set

        OUTPUT:

        If ``self`` is a `Q`-Tree, the sequence of its children is
        returned. If ``self`` is a `P`-tree, 2 `P`-tree are returned,
        namely the two `P`-tree defined above and restricting the
        permutations, in the order implied by ``left, right`` (if
        ``right =True``, the second `P`-tree will be the one gathering
        the elements containing ``v``, if ``left=True``, the
        opposite).

        .. NOTE::

           This method is assumes that ``self`` is partial for ``v``,
           and aligned to the side indicated by ``left, right``.

        EXAMPLES:

        A `P`-Tree ::

            sage: from sage.graphs.pq_trees import P, Q
            sage: p = P([[2,4], [1,2], [0,8], [0,5]])
            sage: p.simplify(0, right = True)
            [('P', [{2, 4}, {1, 2}]), ('P', [{0, 8}, {0, 5}])]

        A `Q`-Tree ::

            sage: q = Q([[2,4], [1,2], [0,8], [0,5]])
            sage: q.simplify(0, right = True)
            [{2, 4}, {1, 2}, {0, 8}, {0, 5}]
        """
        if sum([left, right]) != 1:
            raise ValueError("Exactly one of left or right must be specified")

        if isinstance(self, Q):
            L = []
            for c in self._children:
                if (isinstance(c, PQ) and               # Is c partial?
                        v in c and                      # (does c contain sets with
                        any(v not in cc for cc in c)):  # and without v ?)
                    L.extend(c.simplify(v, right=right, left=left))
                else:
                    L.append(c)
            return L
        else:
            empty = []
            full = []
            partial = []

            for c in self._children:
                if v in c:
                    if (isinstance(c, PQ) and               # Is c partial? (does c contain
                            any(v not in cc for cc in c)):  # sets with and without v ?)
                        partial = c.simplify(v, right=right, left=left)
                    else:
                        full.append(c)
                else:
                    empty.append(c)
            if empty:
                empty = [_new_P(empty)]
            if full:
                full = [_new_P(full)]

            if right:
                return empty + partial + full
            else:
                return full + partial + empty

    def flatten(self):
        r"""
        Return a flattened copy of ``self``.

        If ``self`` has only one child, we may as well consider its
        child's children, as ``self`` encodes no information. This
        method recursively "flattens" trees having only on PQ-tree
        child, and returns it.

        EXAMPLES::

            sage: from sage.graphs.pq_trees import P, Q
            sage: p = Q([P([[2,4], [2,8], [2,9]])])
            sage: p.flatten()
            ('P', [{2, 4}, {8, 2}, {9, 2}])
        """
        if self.number_of_children() == 1:
            return _flatten(self._children[0])
        else:
            self._children = [_flatten(x) for x in self._children]
            return self


class P(PQ):
    r"""
    A P-Tree is a PQ-Tree whose children can be permuted in any way.

    For more information, see the documentation of :mod:`sage.graphs.pq_trees`.
    """
    def set_contiguous(self, v):
        r"""
        Update ``self`` so that the sets containing ``v`` are
        contiguous for any admissible permutation of its subtrees.

        INPUT:

        - ``v`` -- an element of the ground set

        OUTPUT:

        According to the cases:

            * ``(EMPTY, ALIGNED)`` if no set of the tree contains
              an occurrence of ``v``

            * ``(FULL, ALIGNED)`` if all the sets of the tree contain
              ``v``

            * ``(PARTIAL, ALIGNED)`` if some (but not all) of the sets
              contain ``v``, all of which are aligned
              to the right of the ordering at the end when the function ends

            * ``(PARTIAL, UNALIGNED)`` if some (but not all) of the
              sets contain ``v``, though it is impossible to align them
              all to the right

        In any case, the sets containing ``v`` are contiguous when this
        function ends. If there is no possibility of doing so, the function
        raises a :exc:`ValueError` exception.

        EXAMPLES:

        Ensuring the sets containing ``0`` are continuous::

            sage: from sage.graphs.pq_trees import P, Q
            sage: p = P([[0,3], [1,2], [2,3], [2,4], [4,0],[2,8], [2,9]])
            sage: p.set_contiguous(0)
            (1, True)
            sage: print(p)
            ('P', [{1, 2}, {2, 3}, {2, 4}, {8, 2}, {9, 2}, ('P', [{0, 3}, {0, 4}])])

        Impossible situation::

            sage: p = P([[0,1], [1,2], [2,3], [3,0]])
            sage: p.set_contiguous(0)
            (1, True)
            sage: p.set_contiguous(1)
            (1, True)
            sage: p.set_contiguous(2)
            (1, True)
            sage: p.set_contiguous(3)
            Traceback (most recent call last):
            ...
            ValueError: Impossible
        """

        ###############################################################
        # Defining Variables :                                        #
        #                                                             #
        # Collecting the information of which children are FULL of v, #
        # which ones are EMPTY, PARTIAL_ALIGNED and PARTIAL_UNALIGNED #
        #                                                             #
        # Defining variables for their cardinals, just to make the    #
        # code slightly more readable :-)                             #
        ###############################################################

        for x in self:
            _set_contiguous(x, v)
        self.flatten()
        seq = [_set_contiguous(x, v) for x in self]

        f_seq = dict(zip(self, seq))

        set_FULL = []
        set_EMPTY = []
        set_PARTIAL_ALIGNED = []
        set_PARTIAL_UNALIGNED = []

        sorting = {
            (FULL, ALIGNED): set_FULL,
            (EMPTY, ALIGNED): set_EMPTY,
            (PARTIAL, ALIGNED): set_PARTIAL_ALIGNED,
            (PARTIAL, UNALIGNED): set_PARTIAL_UNALIGNED
            }

        for i in self:
            sorting[f_seq[i]].append(i)

        n_FULL = len(set_FULL)
        n_EMPTY = len(set_EMPTY)
        n_PARTIAL_ALIGNED = len(set_PARTIAL_ALIGNED)
        n_PARTIAL_UNALIGNED = len(set_PARTIAL_UNALIGNED)

        # Excludes the situation where there is no solution.
        # read next comment for more explanations

        if (n_PARTIAL_ALIGNED > 2 or
                (n_PARTIAL_UNALIGNED >= 1 and n_EMPTY != self.number_of_children() - 1)):
            raise ValueError(impossible_msg)

        # From now on, there are at most two pq-trees which are partially filled
        # If there is one which is not aligned to the right, all the others are empty

        #########################################################
        # 1/2                                                   #
        #                                                       #
        # Several easy cases where we can decide without paying #
        # attention                                             #
        #########################################################

        # All the children are FULL
        elif n_FULL == self.number_of_children():
            return FULL, True

        # All the children are empty
        elif n_EMPTY == self.number_of_children():
            return EMPTY, True

        # There is a PARTIAL UNALIGNED element (and all the others are
        # empty as we checked before

        elif n_PARTIAL_UNALIGNED == 1:
            return (PARTIAL, UNALIGNED)

        # If there is just one partial element and all the others are
        # empty, we just reorder the set to put it at the right end

        elif (n_PARTIAL_ALIGNED == 1 and
              n_EMPTY == self.number_of_children()-1):

            self._children = set_EMPTY + set_PARTIAL_ALIGNED
            return (PARTIAL, ALIGNED)

        ################################################################
        # 2/2                                                          #
        #                                                              #
        # From now on, there are at most two partial pq-trees and all  #
        # of them have v aligned to their right                        #
        #                                                              #
        # We now want to order them in such a way that all the         #
        # elements containing v are located on the right               #
        ################################################################

        else:

            self._children = []

            # We first move the empty elements to the left, if any

            if n_EMPTY > 0:
                self._children.extend(set_EMPTY)

            # If there is one partial element we but have to add it to
            # the sequence, then add all the full elements

            # We must also make sure these elements will not be
            # reordered in such a way that the elements containing v
            # are not contiguous

            # ==> We create a Q-tree

            if n_PARTIAL_ALIGNED < 2:

                new = []

                # add the partial element, if any
                if n_PARTIAL_ALIGNED == 1:

                    subtree = set_PARTIAL_ALIGNED[0]
                    new.extend(subtree.simplify(v, right=ALIGNED))

                # Then the full elements, if any, in a P-tree (we can
                # permute any two of them while keeping all the
                # elements containing v on an interval

                if n_FULL > 0:

                    new.append(_new_P(set_FULL))

                # We lock all of them in a Q-tree

                self._children.append(_new_Q(new))

                return PARTIAL, True

            # If there are 2 partial elements, we take care of both
            # ends. We also know it will not be possible to align the
            # interval of sets containing v to the right

            else:
                new = []

                # The second partial element is aligned to the right
                # while, as we want to put it at the end of the
                # interval, it should be aligned to the left
                set_PARTIAL_ALIGNED[1].reverse()

                # 1/3
                # Left partial subtree
                subtree = set_PARTIAL_ALIGNED[0]
                new.extend(subtree.simplify(v, right=ALIGNED))

                # 2/3
                # Center (Full elements, in a P-tree, as they can be
                # permuted)

                if n_FULL > 0:
                    new.append(_new_P(set_FULL))

                # 3/3
                # Right partial subtree
                subtree = set_PARTIAL_ALIGNED[1]
                new.extend(subtree.simplify(v, left=ALIGNED))

                # We add all of it, locked in a Q-Tree
                self._children.append(_new_Q(new))

                return PARTIAL, False

    def cardinality(self):
        r"""
        Return the number of orderings allowed by the structure.

        .. SEEALSO::

            :meth:`orderings` -- iterate over all admissible orderings

        EXAMPLES::

            sage: from sage.graphs.pq_trees import P, Q
            sage: p = P([[0,3], [1,2], [2,3], [2,4], [4,0],[2,8], [2,9]])
            sage: p.cardinality()
            5040
            sage: p.set_contiguous(3)
            (1, True)
            sage: p.cardinality()
            1440
        """
        from math import factorial
        n = factorial(self.number_of_children())
        for c in self._children:
            if isinstance(c, PQ):
                n = n*c.cardinality()
        return n

    def orderings(self):
        r"""
        Iterate over all orderings of the sets allowed by the structure.

        .. SEEALSO::

            :meth:`cardinality` -- return the number of orderings

        EXAMPLES::

            sage: from sage.graphs.pq_trees import P, Q
            sage: p = P([[2,4], [1,2], [0,8], [0,5]])
            sage: for o in p.orderings():
            ....:    print(o)
            ({2, 4}, {1, 2}, {0, 8}, {0, 5})
            ({2, 4}, {1, 2}, {0, 5}, {0, 8})
            ({2, 4}, {0, 8}, {1, 2}, {0, 5})
            ({2, 4}, {0, 8}, {0, 5}, {1, 2})
            ...
        """
        from itertools import permutations, product
        for p in permutations(self._children):
            yield from product(*[x.orderings() if isinstance(x, PQ) else [x]
                                 for x in p])


class Q(PQ):
    r"""
    A Q-Tree is a PQ-Tree whose children are ordered up to reversal.

    For more information, see the documentation of :mod:`sage.graphs.pq_trees`.
    """

    def set_contiguous(self, v):
        r"""
        Update ``self`` so that the sets containing ``v`` are
        contiguous for any admissible permutation of its subtrees.

        INPUT:

        - ``v`` -- an element of the ground set

        OUTPUT:

        According to the cases:

            * ``(EMPTY, ALIGNED)`` if no set of the tree contains
              an occurrence of ``v``

            * ``(FULL, ALIGNED)`` if all the sets of the tree contain
              ``v``

            * ``(PARTIAL, ALIGNED)`` if some (but not all) of the sets
              contain ``v``, all of which are aligned
              to the right of the ordering at the end when the function ends

            * ``(PARTIAL, UNALIGNED)`` if some (but not all) of the
              sets contain ``v``, though it is impossible to align them
              all to the right

        In any case, the sets containing ``v`` are contiguous when this
        function ends. If there is no possibility of doing so, the function
        raises a :exc:`ValueError` exception.

        EXAMPLES:

        Ensuring the sets containing ``0`` are continuous::

            sage: from sage.graphs.pq_trees import P, Q
            sage: q = Q([[2,3], Q([[3,0],[3,1]]), Q([[4,0],[4,5]])])
            sage: q.set_contiguous(0)
            (1, False)
            sage: print(q)
            ('Q', [{2, 3}, {1, 3}, {0, 3}, {0, 4}, {4, 5}])

        Impossible situation::

            sage: p = Q([[0,1], [1,2], [2,0]])
            sage: p.set_contiguous(0)
            Traceback (most recent call last):
            ...
            ValueError: Impossible
        """
        #################################################################
        # Guidelines :                                                  #
        #                                                               #
        # As the tree is a Q-Tree, we can but reverse the order in      #
        # which the elements appear. It means that we can but check     #
        # the elements containing v are already contiguous (even        #
        # though we have to take special care of partial elements --    #
        # the endpoints of the interval), and answer accordingly        #
        # (partial, full, empty, aligned..). We also want to align the  #
        # elements containing v to the right if possible.               #
        ################################################################

        ###############################################################
        # Defining Variables :                                        #
        #                                                             #
        # Collecting the information of which children are FULL of v, #
        # which ones are EMPTY, PARTIAL_ALIGNED and PARTIAL_UNALIGNED #
        #                                                             #
        # Defining variables for their cardinals, just to make the    #
        # code slightly more readable :-)                             #
        ###############################################################

        for x in self:
            _set_contiguous(x, v)
        self.flatten()
        seq = [_set_contiguous(x, v) for x in self]

        f_seq = dict(zip(self, seq))

        set_FULL = []
        set_EMPTY = []
        set_PARTIAL_ALIGNED = []
        set_PARTIAL_UNALIGNED = []

        sorting = {
            (FULL, ALIGNED): set_FULL,
            (EMPTY, ALIGNED): set_EMPTY,
            (PARTIAL, ALIGNED): set_PARTIAL_ALIGNED,
            (PARTIAL, UNALIGNED): set_PARTIAL_UNALIGNED
            }

        for i in self:
            sorting[f_seq[i]].append(i)

        n_FULL = len(set_FULL)
        n_EMPTY = len(set_EMPTY)
        n_PARTIAL_ALIGNED = len(set_PARTIAL_ALIGNED)
        n_PARTIAL_UNALIGNED = len(set_PARTIAL_UNALIGNED)

        ###################################################################
        #                                                                 #
        # Picking the good ordering for the children :                    #
        #                                                                 #
        #                                                                 #
        # There is a possibility of aligning to the right iif             #
        # the vector can assume the form (as a regular expression) :      #
        #                                                                 #
        # (EMPTY *) PARTIAL (FULL *) Of course, each of these three       #
        # members could be empty                                          #
        #                                                                 #
        # Hence, in the following case we reverse the vector :            #
        #                                                                 #
        # * if the last element is empty (as we checked the whole         #
        #   vector is not empty                                           #
        #                                                                 #
        # * if the last element is partial, aligned, and all the          #
        #   others are full                                               #
        ###################################################################

        if (f_seq[self._children[-1]] == (EMPTY, ALIGNED) or
            (f_seq[self._children[-1]] == (PARTIAL, ALIGNED) and
             n_FULL == self.number_of_children() - 1)):
            # We reverse the order of the elements in the SET only.
            # Which means that they are still aligned to the right !
            self._children.reverse()

        #########################################################
        # 1/2                                                   #
        #                                                       #
        # Several easy cases where we can decide without paying #
        # attention                                             #
        #########################################################

        # Excludes the situation where there is no solution.
        # read next comment for more explanations

        if (n_PARTIAL_ALIGNED > 2 or
            (n_PARTIAL_UNALIGNED >= 1 and
             n_EMPTY != self.number_of_children() - 1)):
            raise ValueError(impossible_msg)

        # From now on, there are at most two pq-trees which are partially filled
        # If there is one which is not aligned to the right, all the others are empty

        # First trivial case, no checking needed
        elif n_FULL == self.number_of_children():
            return FULL, True

        # Second trivial case, no checking needed
        elif n_EMPTY == self.number_of_children():
            return EMPTY, True

        # Third trivial case, no checking needed
        elif n_PARTIAL_UNALIGNED == 1:
            return (PARTIAL, UNALIGNED)

        # If there is just one partial element
        # and all the others are empty, we just reorder
        # the set to put it at the right end

        elif (n_PARTIAL_ALIGNED == 1 and
              n_EMPTY == self.number_of_children() - 1):

            if set_PARTIAL_ALIGNED[0] == self._children[-1]:
                return (PARTIAL, ALIGNED)

            else:
                return (PARTIAL, UNALIGNED)

        ##############################################################
        # 2/2                                                        #
        #                                                            #
        # We iteratively consider all the children, and check        #
        # that the elements containing v are indeed                  #
        # located on an interval.                                    #
        #                                                            #
        # We are also interested in knowing whether this interval is #
        # aligned to the right                                       #
        #                                                            #
        # Because of the previous tests, we can assume there are at  #
        # most two partial pq-trees and all of them are aligned to   #
        # their right                                                #
        ##############################################################

        else:

            new_children = []

            # Two variables to remember where we are
            # according to the interval

            seen_nonempty = False
            seen_right_end = False

            for i in self:

                type, aligned = f_seq[i]

                # We met an empty element
                if type == EMPTY:

                    # 2 possibilities :
                    #
                    #  * we have NOT met a non-empty element before
                    #    and it just means we are looking at the
                    #    leading empty elements
                    #
                    #  * we have met a non-empty element before and it
                    #    means we will never met another non-empty
                    #    element again => we have seen the right end
                    #    of the interval

                    new_children.append(i)

                    if seen_nonempty:
                        seen_right_end = True

                # We met a non-empty element
                else:
                    if seen_right_end:
                        raise ValueError(impossible_msg)

                    if type == PARTIAL:

                        # if we see an ALIGNED partial tree after
                        # having seen a nonempty element then the
                        # partial tree must be aligned to the left and
                        # so we have seen the right end

                        if seen_nonempty and aligned:
                            i.reverse()
                            seen_right_end = True

                            # right partial subtree
                            subtree = i
                            new_children.extend(subtree.simplify(v, left=True))

                        # If we see an UNALIGNED partial element after
                        # having met a nonempty element, there is no
                        # solution to the alignment problem

                        elif seen_nonempty and not aligned:
                            raise ValueError(impossible_msg)

                        # If we see an unaligned element but no non-empty
                        # element since the beginning, we are witnessing both the
                        # left and right end

                        elif not seen_nonempty and not aligned:
                            raise ValueError("Bon, ben ca arrive O_o")
                            seen_right_end = True

                        elif not seen_nonempty and aligned:

                            # left partial subtree
                            subtree = i

                            new_children.extend(subtree.simplify(v, right=True))

                    else:
                        new_children.append(i)

                    seen_nonempty = True

            # Setting the updated sequence of children
            self._children = new_children

            # Whether we achieved an alignment to the right is the
            # complement of whether we have seen the right end

            return (PARTIAL, not seen_right_end)

    def cardinality(self):
        r"""
        Return the number of orderings allowed by the structure.

        .. SEEALSO::

            :meth:`orderings` -- iterate over all admissible orderings

        EXAMPLES::

            sage: from sage.graphs.pq_trees import P, Q
            sage: q = Q([[0,3], [1,2], [2,3], [2,4], [4,0],[2,8], [2,9]])
            sage: q.cardinality()
            2
        """
        n = 1
        for c in self._children:
            if isinstance(c, PQ):
                n = n*c.cardinality()

        return n if (self.number_of_children() == 1) else 2*n

    def orderings(self):
        r"""
        Iterate over all orderings of the sets allowed by the structure.

        .. SEEALSO::

            :meth:`cardinality` -- return the number of orderings

        EXAMPLES::

            sage: from sage.graphs.pq_trees import P, Q
            sage: q = Q([[2,4], [1,2], [0,8], [0,5]])
            sage: for o in q.orderings():
            ....:    print(o)
            ({2, 4}, {1, 2}, {0, 8}, {0, 5})
            ({0, 5}, {0, 8}, {1, 2}, {2, 4})
        """
        if len(self._children) == 1:
            c = self._children[0]
            yield from (c.orderings() if isinstance(c, PQ) else [c])
        else:
            from itertools import product
            for o in product(*[x.orderings() if isinstance(x, PQ) else [x]
                               for x in self._children]):
                yield o
                yield o[::-1]
