# -*- coding: utf-8 -*-
# pylint: disable=all

# TODO: (high):
#
# check whether it makes sense to keep a list of solutions, and keep
# a global MILP up to date with this list

# TODO: (medium):
#
# can we somehow tweak gurobi so that
# minimal_subdistributions_iterator considers the minimal
# subdistributions with "smallest" (e.g., in the sorting order
# defined above) elements first?

r"""
A bijectionist's toolkit

AUTHORS:

- Alexander Grosz, Tobias Kietreiber, Stephan Pfannerer and Martin
  Rubey (2020): Initial version

Quick reference
===============

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~Bijectionist.set_intertwining_relations` | Set
    :meth:`~Bijectionist.set_constant_blocks` | Set
    :meth:`~Bijectionist.set_statistics` | Set
    :meth:`~Bijectionist.set_value_restrictions` | Set
    :meth:`~Bijectionist.set_distributions` | Set

    :meth:`~Bijectionist.statistics_table` | Return
    :meth:`~Bijectionist.statistics_fibers` | Return

    :meth:`~Bijectionist.constant_blocks` | Return
    :meth:`~Bijectionist.solutions_iterator` | Return
    :meth:`~Bijectionist.possible_values` | Return
    :meth:`~Bijectionist.minimal_subdistributions_iterator` | Return
    :meth:`~Bijectionist.minimal_subdistributions_blocks_iterator` | Return

A guided tour
=============

    EXAMPLES:

    We find a statistic `s` such that
    `(s, wex, fix) \sim (llis, des, adj)`::

        sage: N = 3
        sage: As = [list(Permutations(n)) for n in range(N+1)]
        sage: A = B = sum(As, [])
        sage: alpha1 = lambda p: len(p.weak_excedences())
        sage: alpha2 = lambda p: len(p.fixed_points())
        sage: beta1 = lambda p: len(p.descents(final_descent=True)) if p else 0
        sage: beta2 = lambda p: len([e for (e, f) in zip(p, p[1:]+[0]) if e == f+1])
        sage: tau = Permutation.longest_increasing_subsequence_length
        sage: def rotate_permutation(p):
        ....:     cycle = Permutation(tuple(range(1, len(p)+1)))
        ....:     return Permutation([cycle.inverse()(p(cycle(i))) for i in range(1, len(p)+1)])
        sage: bij = Bijectionist(A, B, tau)
        sage: bij.set_statistics((len, len), (alpha1, beta1), (alpha2, beta2))
        sage: a, b = bij.statistics_table()
        sage: table(a, header_row=True, frame=True)
        +-----------+--------+--------+--------+
        | a         | α_1(a) | α_2(a) | α_3(a) |
        +===========+========+========+========+
        | []        | 0      | 0      | 0      |
        +-----------+--------+--------+--------+
        | [1]       | 1      | 1      | 1      |
        +-----------+--------+--------+--------+
        | [1, 2]    | 2      | 2      | 2      |
        +-----------+--------+--------+--------+
        | [2, 1]    | 2      | 1      | 0      |
        +-----------+--------+--------+--------+
        | [1, 2, 3] | 3      | 3      | 3      |
        +-----------+--------+--------+--------+
        | [1, 3, 2] | 3      | 2      | 1      |
        +-----------+--------+--------+--------+
        | [2, 1, 3] | 3      | 2      | 1      |
        +-----------+--------+--------+--------+
        | [2, 3, 1] | 3      | 2      | 0      |
        +-----------+--------+--------+--------+
        | [3, 1, 2] | 3      | 1      | 0      |
        +-----------+--------+--------+--------+
        | [3, 2, 1] | 3      | 2      | 1      |
        +-----------+--------+--------+--------+

        sage: table(b, header_row=True, frame=True)
        +-----------+---+--------+--------+--------+
        | b         | τ | β_1(b) | β_2(b) | β_3(b) |
        +===========+===+========+========+========+
        | []        | 0 | 0      | 0      | 0      |
        +-----------+---+--------+--------+--------+
        | [1]       | 1 | 1      | 1      | 1      |
        +-----------+---+--------+--------+--------+
        | [1, 2]    | 2 | 2      | 1      | 0      |
        +-----------+---+--------+--------+--------+
        | [2, 1]    | 1 | 2      | 2      | 2      |
        +-----------+---+--------+--------+--------+
        | [1, 2, 3] | 3 | 3      | 1      | 0      |
        +-----------+---+--------+--------+--------+
        | [1, 3, 2] | 2 | 3      | 2      | 1      |
        +-----------+---+--------+--------+--------+
        | [2, 1, 3] | 2 | 3      | 2      | 1      |
        +-----------+---+--------+--------+--------+
        | [2, 3, 1] | 2 | 3      | 2      | 1      |
        +-----------+---+--------+--------+--------+
        | [3, 1, 2] | 2 | 3      | 2      | 0      |
        +-----------+---+--------+--------+--------+
        | [3, 2, 1] | 1 | 3      | 3      | 3      |
        +-----------+---+--------+--------+--------+

        sage: from sage.combinat.cyclic_sieving_phenomenon import orbit_decomposition
        sage: bij.set_constant_blocks(sum([orbit_decomposition(A, rotate_permutation) for A in As], []))
        sage: bij.constant_blocks()
        {{[1, 3, 2], [2, 1, 3], [3, 2, 1]}}
        sage: next(bij.solutions_iterator())
        {[]: 0,
         [1]: 1,
         [1, 2]: 1,
         [1, 2, 3]: 1,
         [1, 3, 2]: 2,
         [2, 1]: 2,
         [2, 1, 3]: 2,
         [2, 3, 1]: 2,
         [3, 1, 2]: 3,
         [3, 2, 1]: 2}

    There is no rotation invariant statistic on non crossing set partitions which is equidistributed
    with the Strahler number on ordered trees::

        sage: N=8; As = [[SetPartition(d.to_noncrossing_partition()) for d in DyckWords(n)] for n in range(N)]
        sage: A = sum(As, [])
        sage: B = sum([list(OrderedTrees(n)) for n in range(1, N+1)], [])
        sage: theta = lambda m: SetPartition([[i % m.size() + 1 for i in b] for b in m])

    The following code is equivalent to ``tau = findstat(397)``::

        sage: def tau(T):
        ....:     if len(T) == 0:
        ....:         return 1
        ....:     else:
        ....:         l = [tau(S) for S in T]
        ....:         m = max(l)
        ....:         if l.count(m) == 1:
        ....:             return m
        ....:         else:
        ....:             return m+1
        sage: bij = Bijectionist(A, B, tau)
        sage: bij.set_statistics((lambda a: a.size(), lambda b: b.node_number()-1))
        sage: from sage.combinat.cyclic_sieving_phenomenon import orbit_decomposition
        sage: bij.set_constant_blocks(sum([orbit_decomposition(A_n, theta) for A_n in As], []))
        sage: list(bij.solutions_iterator())
        []

    An example identifying `s` and `S`::

        sage: N = 4
        sage: A = [dyck_word for n in range(1, N) for dyck_word in DyckWords(n)]
        sage: B = [binary_tree for n in range(1, N) for binary_tree in BinaryTrees(n)]
        sage: concat_path = lambda D1, D2: DyckWord(list(D1) + list(D2))
        sage: concat_tree = lambda B1, B2: concat_path(B1.to_dyck_word(),
        ....:                                          B2.to_dyck_word()).to_binary_tree()
        sage: bij = Bijectionist(A, B)
        sage: bij.set_intertwining_relations((2, concat_path, concat_tree))
        sage: bij.set_statistics((lambda d: d.semilength(), lambda t: t.node_number()))
        sage: for D in bij.minimal_subdistributions_iterator():
        ....:     ascii_art(D)
        ( [ /\ ], [ o ] )
        (           [ o   ] )
        (           [  \  ] )
        ( [ /\/\ ], [   o ] )
        (           [   o ] )
        ( [  /\  ]  [  /  ] )
        ( [ /  \ ], [ o   ] )
        (             [ o     ] )
        (             [  \    ] )
        (             [   o   ] )
        (             [    \  ] )
        ( [ /\/\/\ ], [     o ] )
        (             [ o   ] )
        (             [  \  ] )
        (             [   o ] )
        ( [    /\  ]  [  /  ] )
        ( [ /\/  \ ], [ o   ] )
        (             [   o   ] )
        ( [  /\    ]  [  / \  ] )
        ( [ /  \/\ ], [ o   o ] )
        (                     [   o,     o ] )
        (                     [  /      /  ] )
        ( [           /\   ]  [ o      o   ] )
        ( [  /\/\    /  \  ]  [  \    /    ] )
        ( [ /    \, /    \ ], [   o  o     ] )

    TESTS:

    The following failed before commit c6d4d2e8804aa42afa08c72c887d50c725cc1a91::

        sage: N=4; A = B = [permutation for n in range(N) for permutation in Permutations(n)]
        sage: theta = lambda pi: Permutation([x+1 if x != len(pi) else 1 for x in pi[-1:]+pi[:-1]])
        sage: def tau(pi):
        ....:    n = len(pi)
        ....:    return sum([1 for i in range(1, n+1) for j in range(1, n+1)
        ....:                if i<j <= pi(i)<pi(j) or pi(i)<pi(j)<i<j])
        sage: bij = Bijectionist(A, B, tau)
        sage: bij.set_statistics((len, len))
        sage: from sage.combinat.cyclic_sieving_phenomenon import orbit_decomposition
        sage: bij.set_constant_blocks(orbit_decomposition(A, theta))
        sage: for solution in sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items()))):
        ....:     print(solution)
        {[]: 0, [1]: 0, [1, 2]: 0, [2, 1]: 0, [1, 2, 3]: 0, [1, 3, 2]: 0, [2, 1, 3]: 0, [3, 2, 1]: 0, [2, 3, 1]: 0, [3, 1, 2]: 1}
        {[]: 0, [1]: 0, [1, 2]: 0, [2, 1]: 0, [1, 2, 3]: 0, [1, 3, 2]: 0, [2, 1, 3]: 0, [3, 2, 1]: 0, [2, 3, 1]: 1, [3, 1, 2]: 0}
        {[]: 0, [1]: 0, [1, 2]: 0, [2, 1]: 0, [1, 2, 3]: 1, [1, 3, 2]: 0, [2, 1, 3]: 0, [3, 2, 1]: 0, [2, 3, 1]: 0, [3, 1, 2]: 0}

    A test including intertwining relations::

        sage: N = 2; A = B = [dyck_word for n in range(N+1) for dyck_word in DyckWords(n)]
        sage: alpha = lambda D: (D.area(), D.bounce())
        sage: beta = lambda D: (D.bounce(), D.area())
        sage: tau = lambda D: D.number_of_touch_points()

    The following looks correct::

        sage: bij = Bijectionist(A, B, tau)
        sage: bij.set_statistics((lambda d: d.semilength(), lambda d: d.semilength()))
        sage: for solution in sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items()))):
        ....:     print(solution)
        {[]: 0, [1, 0]: 1, [1, 0, 1, 0]: 1, [1, 1, 0, 0]: 2}
        {[]: 0, [1, 0]: 1, [1, 0, 1, 0]: 2, [1, 1, 0, 0]: 1}

    The following looks correct, because `alpha = beta \circ S`
    forces `S([1,0,1,0]) = [1,1,0,0]` and `s = tau \circ S` forces
    therefore `s([1,0,1,0]) = \tau(S([1,0,1,0])) =
    \tau([1,1,0,0]) = 1`::

        sage: bij = Bijectionist(A, B, tau)
        sage: bij.set_statistics((alpha, beta), (lambda d: d.semilength(), lambda d: d.semilength()))
        sage: for solution in bij.solutions_iterator():
        ....:     print(solution)
        {[]: 0, [1, 0]: 1, [1, 0, 1, 0]: 1, [1, 1, 0, 0]: 2}

    Now we introduce a intertwining relation::

        sage: concat_path = lambda D1, D2: DyckWord(list(D1) + list(D2))
        sage: pi_rho = (2, concat_path, lambda x, y: x+y)

    Without `\alpha` and `\beta` but with `\pi` and `\rho` the other
    values are forced because `s([1,0,1,0]) = s(\pi([1,0], [1,0])) =
    \rho(s([1,0]), s([1,0])) = 2`::

        sage: bij = Bijectionist(A, B, tau)
        sage: bij.set_intertwining_relations(pi_rho)
        sage: bij.set_statistics((lambda d: d.semilength(), lambda d: d.semilength()))
        sage: for solution in bij.solutions_iterator():
        ....:     print(solution)
        {[]: 0, [1, 0]: 1, [1, 0, 1, 0]: 2, [1, 1, 0, 0]: 1}

    Thus the combination of both constraints should be infeasible::

        sage: bij = Bijectionist(A, B, tau)
        sage: bij.set_statistics((alpha, beta), (lambda d: d.semilength(), lambda d: d.semilength()))
        sage: bij.set_intertwining_relations(pi_rho)
        sage: list(bij.solutions_iterator())
        []

    Repeating some tests, but using the constructor instead of set_XXX() methods:

        sage: N = 2; A = B = [dyck_word for n in range(N+1) for dyck_word in DyckWords(n)]
        sage: alpha = lambda D: (D.area(), D.bounce())
        sage: beta = lambda D: (D.bounce(), D.area())
        sage: tau = lambda D: D.number_of_touch_points()

        sage: bij = Bijectionist(A, B, tau, alpha_beta=((lambda d: d.semilength(), lambda d: d.semilength()),))
        sage: for solution in sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items()))):
        ....:     print(solution)
        {[]: 0, [1, 0]: 1, [1, 0, 1, 0]: 1, [1, 1, 0, 0]: 2}
        {[]: 0, [1, 0]: 1, [1, 0, 1, 0]: 2, [1, 1, 0, 0]: 1}


    Constant blocks::

        sage: A = B = list('abcd')
        sage: pi = lambda p1, p2: 'abcdefgh'[A.index(p1) + A.index(p2)]
        sage: rho = lambda s1, s2: (s1 + s2) % 2
        sage: bij = Bijectionist(A, B, lambda x: B.index(x) % 2, P=[['a', 'c']], pi_rho=((2, pi, rho),))
        sage: list(bij.solutions_iterator())
        [{'a': 0, 'b': 1, 'c': 0, 'd': 1}]
        sage: bij.constant_blocks()
        {{'a', 'c'}, {'b', 'd'}}


    Distributions::

        sage: A = B = [permutation for n in range(4) for permutation in Permutations(n)]
        sage: tau = Permutation.longest_increasing_subsequence_length
        sage: bij = Bijectionist(A, B, tau, alpha_beta=((len, len),), elements_distributions=(([Permutation([1, 2, 3]), Permutation([1, 3, 2])], [1, 3]),))
        sage: for sol in sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items()))):
        ....:     print(sol)
        {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 1, [1, 3, 2]: 3, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 2}
        {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 3, [1, 3, 2]: 1, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 2}
        {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 1, [1, 2, 3]: 1, [1, 3, 2]: 3, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 2}
        {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 1, [1, 2, 3]: 3, [1, 3, 2]: 1, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 2}


    Intertwining relations::

        sage: concat = lambda p1, p2: Permutation(p1 + [i + len(p1) for i in p2])

        sage: A = B = [permutation for n in range(4) for permutation in Permutations(n)]
        sage: bij = Bijectionist(A, B, Permutation.number_of_fixed_points, alpha_beta=((len, len),), pi_rho=((2, concat, lambda x, y: x + y),))
        sage: for solution in sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items()))):
        ....:     print(solution)
        {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 0, [1, 2, 3]: 3, [1, 3, 2]: 1, [2, 1, 3]: 1, [2, 3, 1]: 0, [3, 1, 2]: 0, [3, 2, 1]: 1}
        {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 0, [1, 2, 3]: 3, [1, 3, 2]: 1, [2, 1, 3]: 1, [2, 3, 1]: 0, [3, 1, 2]: 1, [3, 2, 1]: 0}
        {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 0, [1, 2, 3]: 3, [1, 3, 2]: 1, [2, 1, 3]: 1, [2, 3, 1]: 1, [3, 1, 2]: 0, [3, 2, 1]: 0}


    Statistics::

        sage: N = 4; A = B = [permutation for n in range(N) for permutation in Permutations(n)]
        sage: wex = lambda p: len(p.weak_excedences())
        sage: fix = lambda p: len(p.fixed_points())
        sage: des = lambda p: len(p.descents(final_descent=True)) if p else 0
        sage: bij = Bijectionist(A, B, fix, alpha_beta=((wex, des), (len, len)))
        sage: for solution in sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items()))):
        ....:     print(solution)
        {[]: 0, [1]: 1, [1, 2]: 0, [2, 1]: 2, [1, 2, 3]: 1, [1, 3, 2]: 0, [2, 1, 3]: 0, [2, 3, 1]: 1, [3, 1, 2]: 3, [3, 2, 1]: 1}
        {[]: 0, [1]: 1, [1, 2]: 0, [2, 1]: 2, [1, 2, 3]: 1, [1, 3, 2]: 0, [2, 1, 3]: 1, [2, 3, 1]: 0, [3, 1, 2]: 3, [3, 2, 1]: 1}
        {[]: 0, [1]: 1, [1, 2]: 0, [2, 1]: 2, [1, 2, 3]: 1, [1, 3, 2]: 0, [2, 1, 3]: 1, [2, 3, 1]: 1, [3, 1, 2]: 3, [3, 2, 1]: 0}
        {[]: 0, [1]: 1, [1, 2]: 0, [2, 1]: 2, [1, 2, 3]: 1, [1, 3, 2]: 1, [2, 1, 3]: 0, [2, 3, 1]: 0, [3, 1, 2]: 3, [3, 2, 1]: 1}
        {[]: 0, [1]: 1, [1, 2]: 0, [2, 1]: 2, [1, 2, 3]: 1, [1, 3, 2]: 1, [2, 1, 3]: 0, [2, 3, 1]: 1, [3, 1, 2]: 3, [3, 2, 1]: 0}
        {[]: 0, [1]: 1, [1, 2]: 0, [2, 1]: 2, [1, 2, 3]: 1, [1, 3, 2]: 1, [2, 1, 3]: 1, [2, 3, 1]: 0, [3, 1, 2]: 3, [3, 2, 1]: 0}

    Value restrictions::
        sage: A = B = [permutation for n in range(4) for permutation in Permutations(n)]
        sage: tau = Permutation.longest_increasing_subsequence_length
        sage: bij = Bijectionist(A, B, tau, alpha_beta=((len, len),), a_values=((Permutation([1, 2]), [1]),
        ....:                            (Permutation([3, 2, 1]), [2, 3, 4]),))
        sage: for sol in sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items()))):
        ....:     print(sol)
        {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 1, [1, 3, 2]: 2, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 3}
        ...
        {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 3, [1, 3, 2]: 2, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 1, [3, 2, 1]: 2}

        sage: A = B = [permutation for n in range(4) for permutation in Permutations(n)]
        sage: tau = Permutation.longest_increasing_subsequence_length
        sage: bij = Bijectionist(A, B, tau, a_values=((Permutation([1, 2]), [4, 5]),))
        Traceback (most recent call last):
        ...
        ValueError: No possible values found for singleton block [[1, 2]]

"""
# ****************************************************************************
#       Copyright (C) 2020 Martin Rubey <martin.rubey at tuwien.ac.at>
#                          Stephan Pfannerer
#                          Tobias Kietreiber
#                          Alexander Grosz
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ***************************************************************************
import itertools
from collections import namedtuple
from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException
from sage.rings.integer_ring import ZZ
from sage.combinat.set_partition import SetPartition
from sage.sets.disjoint_set import DisjointSet
from sage.structure.sage_object import SageObject
from copy import copy, deepcopy
from sage.misc.verbose import get_verbose

# TODO: (low) we frequently need variable names for subsets of A, B,
# Z.  In LaTeX, we mostly call them \tilde A, \tilde Z, etc. now.  It
# would be good to have a standard name in code, too.

# TODO: (medium) whenever possible, doctests of a method should only
# test this method.  Currently we have very many system tests, which
# is inconvenient when modifying the design substantially.


class Bijectionist(SageObject):
    r"""Solver class for bijection-statistic problems.

    INPUT:

    - ``A``, ``B`` -- sets of equal size, given as a list

    - ``tau`` (optional, default: ``None``) -- a function from ``B``
      to ``Z``, in case of ``None``, the identity map ``lambda x: x``
      is used

    - ``alpha`` (optional) -- a statistic from ``A`` to ``W``

    - ``beta`` (optional) -- a statistic from ``B`` to ``W``

    - ``P`` (optional) -- a partition of ``A``

    - ``pi_rho`` (optional) -- a triple ``(k, pi, rho)`` where

        - ``pi`` is a ``k``-ary operation composing objects in ``A``
          and

        - ``rho`` is a ``k``-ary function composing statistic values
          in `Z`

    ``W`` and ``Z`` can be arbitrary sets.  As a natural example we
    may think of the natural numbers or tuples of integers.

    We are looking for a statistic `s: A\to Z` and a bijection `S:
    A\to B` such that

    - `s = \tau \circ S`: the statistics `s` and `\tau` are
      equidistributed and `S` is an intertwining bijection.

    - `\alpha = \beta \circ S`: the statistics `\alpha` and `\beta`
      are equidistributed and `S` is an intertwining bijection.

    - `s` is constant on the blocks of `P`.

    - `s(\pi(a_1,\dots, a_k)) = \rho(s(a_1),\dots, s(a_k))`.

    Additionally, we may require that

    - `s(a)\in Z_a` for specified sets `Z_a\subseteq Z`, and

    - `s|_{\tilde A}` has a specified distribution for specified sets
      `\tilde A \subset A`.

    If `\tau` is the identity, the two unknown functions `s` and `S`
    coincide.  Although we do not exclude other bijective choices for
    `\tau`, they probably do not make sense.

    If we want that `S` is graded, i.e. if elements of `A` and `B`
    have a notion of size and `S` should preserve this size, we can
    add grading statistics as `\alpha` and `\beta`.  Since `\alpha`
    and `\beta` will be equidistributed with `S` as an intertwining
    bijection, `S` will then also be graded.

    In summary, we have the following two commutative diagrams, where
    `s` and `S` are unknown functions.

    .. MATH::

        \begin{array}{rrl}
                                          & A \\
            {\scriptstyle\alpha}\swarrow  & {\scriptstyle S}\downarrow & \searrow{\scriptstyle s}\\
            W \overset{\beta}{\leftarrow} & B                          & \overset{\tau}{\rightarrow} Z
        \end{array}
        \qquad
        \begin{array}{lcl}
            A^k                          &\overset{\pi}{\rightarrow} & A\\
            \downarrow{\scriptstyle s^k} &                           & \downarrow{\scriptstyle s}\\
            Z^k                          &\overset{\rho}{\rightarrow} & Z\\
        \end{array}

    .. NOTE::

        If `\tau` is the identity map, the partition `P` of `A`
        necessarily consists only of singletons.

    .. NOTE::

        The order of invocation of the methods with prefix ``set``,
        i.e., :meth:`set_statistics`,
        :meth:`set_intertwining_relations`,
        :meth:`set_constant_blocks`, etc., is irrelevant.  Calling
        any of these methods a second time overrides the previous
        specification.

    """
    def __init__(self, A, B, tau=None, alpha_beta=tuple(), P=[], pi_rho=tuple(), elements_distributions=tuple(), a_values=tuple(), solver=None, key=None):
        """
        Initialize the bijectionist.

        TESTS:

        Check that large input sets are handled well::

            sage: A = B = list(range(7000))
            sage: bij = Bijectionist(A, B)
        """
        # glossary of standard letters:
        # A, B, Z, W ... finite sets
        # ???? tilde_A, tilde_Z, ..., subsets?
        # P ... set partition of A
        # a in A, b in B, p in P
        # S: A -> B
        # alpha: A -> W, beta: B -> W
        # s: A -> Z, tau: B -> Z
        # k arity of pi and rho
        # pi: A^k -> A, rho: Z^k -> Z
        # a_tuple in A^k

        assert len(A) == len(set(A)), "A must have distinct items"
        assert len(B) == len(set(B)), "B must have distinct items"
        self._A = A
        self._B = B
        self._sorter = {}
        self._sorter["A"] = lambda x: sorted(x, key=self._A.index)
        self._sorter["B"] = lambda x: sorted(x, key=self._B.index)

        if tau is None:
            self._tau = {b: b for b in self._B}
        else:
            self._tau = {b: tau(b) for b in self._B}
        self._Z = set(self._tau.values())
        if key is not None and "Z" in key:
            self._sorter["Z"] = lambda x: sorted(x, key=key["Z"])
            self._Z = self._sorter["Z"](self._Z)
        else:
            try:
                self._Z = sorted(self._Z)
                self._sorter["Z"] = lambda x: sorted(x)
            except TypeError:
                self._sorter["Z"] = lambda x: list(x)
                self._Z = list(self._Z)

        # set optional inputs
        self.set_statistics(*alpha_beta)
        self.set_value_restrictions(*a_values)
        self.set_distributions(*elements_distributions)
        self.set_intertwining_relations(*pi_rho)
        self.set_constant_blocks(P)

        self._solver = solver

    def set_constant_blocks(self, P):
        r"""
        Declare that `s: A\to Z` is constant on each block of `P`.

        .. WARNING::

             Any restriction imposed by a previous invocation of
             :meth:`set_constant_blocks` will be overwritten,
             including restrictions discovered by
             :meth:`set_intertwining_relations` and
             :meth:`solutions_iterator`!

        A common example is to use the orbits of a bijection acting
        on `A`.  This can be achieved using the function
        :meth:`~sage.combinat.cyclic_sieving_phenomenon.orbit_decomposition`.

        INPUT:

        - ``P`` -- a set partition of `A`, singletons may be omitted

        EXAMPLES:

        Initially the partitions are set to singleton blocks.  The
        current partition can be reviewed using
        :meth:`constant_blocks`::

            sage: A = B = list('abcd')
            sage: bij = Bijectionist(A, B, lambda x: B.index(x) % 2)
            sage: bij.constant_blocks()
            {}

            sage: bij.set_constant_blocks([['a', 'c']])
            sage: bij.constant_blocks()
            {{'a', 'c'}}

        We now add a map that combines some blocks::

            sage: pi = lambda p1, p2: 'abcdefgh'[A.index(p1) + A.index(p2)]
            sage: rho = lambda s1, s2: (s1 + s2) % 2
            sage: bij.set_intertwining_relations((2, pi, rho))
            sage: list(bij.solutions_iterator())
            [{'a': 0, 'b': 1, 'c': 0, 'd': 1}]
            sage: bij.constant_blocks()
            {{'a', 'c'}, {'b', 'd'}}

        Setting constant blocks overrides any previous assignment::

            sage: bij.set_constant_blocks([['a', 'b']])
            sage: bij.constant_blocks()
            {{'a', 'b'}}

        If there is no solution, and the coarsest partition is
        requested, an error is raised::

            sage: bij.constant_blocks(optimal=True)
            Traceback (most recent call last):
            ...
            MIPSolverException: ...

        """
        self._P = DisjointSet(self._A)
        P = sorted(self._sorter["A"](p) for p in P)
        for p in P:
            for a in p:
                self._P.union(p[0], a)

        self._compute_possible_block_values()

    def constant_blocks(self, singletons=False, optimal=False):
        r"""
        Return the set partition `P` of `A` such that `s: A\to Z` is
        known to be constant on the blocks of `P`.

        INPUT:

        - ``singletons`` (optional, default: ``False``) -- whether or
          not to include singleton blocks in the output

        - ``optimal`` (optional, default: ``False``) -- whether or
          not to compute the coarsest possible partition

        .. NOTE::

            computing the coarsest possible partition may be
            computationally expensive, but may speed up generating
            solutions.

        EXAMPLES::

            sage: A = B = ["a", "b", "c"]
            sage: bij = Bijectionist(A, B, lambda x: 0)
            sage: bij.set_constant_blocks([["a", "b"]])
            sage: bij.constant_blocks()
            {{'a', 'b'}}

            sage: bij.constant_blocks(singletons=True)
            {{'a', 'b'}, {'c'}}

        """
        if optimal:
            self._forced_constant_blocks()
        if singletons:
            return SetPartition(self._P)
        return SetPartition(p for p in self._P if len(p) > 1)

    def set_statistics(self, *alpha_beta):
        r"""
        Set constraints of the form `\alpha = \beta\circ S`.

        .. WARNING::

             Any restriction imposed by a previous invocation of
             :meth:`set_statistics` will be overwritten!

        INPUT:

        - ``alpha_beta`` -- one or more pairs `(\alpha: A\to W,
          \beta: B\to W)`

        If the statistics `\alpha` and `\beta` are not
        equidistributed, an error is raised.

        EXAMPLES:

        We look for bijections `S` on permutations such that the
        number of weak exceedences of `S(\pi)` equals the number of
        descents of `\pi`, and statistics `s`, such that the number
        of fixed points of `S(\pi)` equals `s(\pi)`::

            sage: N = 4; A = B = [permutation for n in range(N) for permutation in Permutations(n)]
            sage: wex = lambda p: len(p.weak_excedences())
            sage: fix = lambda p: len(p.fixed_points())
            sage: des = lambda p: len(p.descents(final_descent=True)) if p else 0
            sage: adj = lambda p: len([e for (e, f) in zip(p, p[1:]+[0]) if e == f+1])
            sage: bij = Bijectionist(A, B, fix)
            sage: bij.set_statistics((wex, des), (len, len))
            sage: for solution in sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items()))):
            ....:     print(solution)
            {[]: 0, [1]: 1, [1, 2]: 0, [2, 1]: 2, [1, 2, 3]: 1, [1, 3, 2]: 0, [2, 1, 3]: 0, [2, 3, 1]: 1, [3, 1, 2]: 3, [3, 2, 1]: 1}
            {[]: 0, [1]: 1, [1, 2]: 0, [2, 1]: 2, [1, 2, 3]: 1, [1, 3, 2]: 0, [2, 1, 3]: 1, [2, 3, 1]: 0, [3, 1, 2]: 3, [3, 2, 1]: 1}
            {[]: 0, [1]: 1, [1, 2]: 0, [2, 1]: 2, [1, 2, 3]: 1, [1, 3, 2]: 0, [2, 1, 3]: 1, [2, 3, 1]: 1, [3, 1, 2]: 3, [3, 2, 1]: 0}
            {[]: 0, [1]: 1, [1, 2]: 0, [2, 1]: 2, [1, 2, 3]: 1, [1, 3, 2]: 1, [2, 1, 3]: 0, [2, 3, 1]: 0, [3, 1, 2]: 3, [3, 2, 1]: 1}
            {[]: 0, [1]: 1, [1, 2]: 0, [2, 1]: 2, [1, 2, 3]: 1, [1, 3, 2]: 1, [2, 1, 3]: 0, [2, 3, 1]: 1, [3, 1, 2]: 3, [3, 2, 1]: 0}
            {[]: 0, [1]: 1, [1, 2]: 0, [2, 1]: 2, [1, 2, 3]: 1, [1, 3, 2]: 1, [2, 1, 3]: 1, [2, 3, 1]: 0, [3, 1, 2]: 3, [3, 2, 1]: 0}

            sage: bij = Bijectionist(A, B, fix)
            sage: bij.set_statistics((wex, des), (fix, adj), (len, len))
            sage: for solution in sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items()))):
            ....:     print(solution)
            {[]: 0, [1]: 1, [1, 2]: 0, [2, 1]: 2, [1, 2, 3]: 1, [1, 3, 2]: 0, [2, 1, 3]: 1, [2, 3, 1]: 0, [3, 1, 2]: 3, [3, 2, 1]: 1}
            {[]: 0, [1]: 1, [1, 2]: 0, [2, 1]: 2, [1, 2, 3]: 1, [1, 3, 2]: 1, [2, 1, 3]: 0, [2, 3, 1]: 0, [3, 1, 2]: 3, [3, 2, 1]: 1}
            {[]: 0, [1]: 1, [1, 2]: 0, [2, 1]: 2, [1, 2, 3]: 1, [1, 3, 2]: 1, [2, 1, 3]: 1, [2, 3, 1]: 0, [3, 1, 2]: 3, [3, 2, 1]: 0}

        Calling this with non-equidistributed statistics yields an error::

            sage: bij = Bijectionist(A, B, fix)
            sage: bij.set_statistics((wex, fix))
            Traceback (most recent call last):
            ...
            ValueError: Statistics alpha and beta are not equidistributed!

        TESTS:

        Calling ``set_statistics`` without arguments should restore the previous state.::

            sage: N = 3; A = B = [permutation for n in range(N) for permutation in Permutations(n)]
            sage: wex = lambda p: len(p.weak_excedences())
            sage: fix = lambda p: len(p.fixed_points())
            sage: des = lambda p: len(p.descents(final_descent=True)) if p else 0
            sage: bij = Bijectionist(A, B, fix)
            sage: bij.set_statistics((wex, des), (len, len))
            sage: for solution in bij.solutions_iterator():
            ....:     print(solution)
            {[]: 0, [1]: 1, [1, 2]: 0, [2, 1]: 2}
            sage: bij.set_statistics()
            sage: for solution in sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items()))):
            ....:     print(solution)
            {[]: 0, [1]: 0, [1, 2]: 1, [2, 1]: 2}
            {[]: 0, [1]: 0, [1, 2]: 2, [2, 1]: 1}
            {[]: 0, [1]: 1, [1, 2]: 0, [2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 0}
            {[]: 0, [1]: 2, [1, 2]: 0, [2, 1]: 1}
            {[]: 0, [1]: 2, [1, 2]: 1, [2, 1]: 0}
            {[]: 1, [1]: 0, [1, 2]: 0, [2, 1]: 2}
            {[]: 1, [1]: 0, [1, 2]: 2, [2, 1]: 0}
            {[]: 1, [1]: 2, [1, 2]: 0, [2, 1]: 0}
            {[]: 2, [1]: 0, [1, 2]: 0, [2, 1]: 1}
            {[]: 2, [1]: 0, [1, 2]: 1, [2, 1]: 0}
            {[]: 2, [1]: 1, [1, 2]: 0, [2, 1]: 0}

        """
        self._n_statistics = len(alpha_beta)
        # TODO: (low) do we really want to recompute statistics every time?
        self._alpha = lambda p: tuple(arg[0](p) for arg in alpha_beta)
        self._beta = lambda p: tuple(arg[1](p) for arg in alpha_beta)

        # generate fibers
        self._statistics_fibers = {}
        for a in self._A:
            v = self._alpha(a)
            if v not in self._statistics_fibers:
                self._statistics_fibers[v] = ([], [])
            self._statistics_fibers[v][0].append(a)

        for b in self._B:
            v = self._beta(b)
            if v not in self._statistics_fibers:
                raise ValueError(f"Statistics alpha and beta do not have the same image, {v} is not a value of alpha, but of beta!")
            self._statistics_fibers[v][1].append(b)

        # check compatibility
        if not all(len(fiber[0]) == len(fiber[1])
                   for fiber in self._statistics_fibers.values()):
            raise ValueError("Statistics alpha and beta are not equidistributed!")

        self._W = list(self._statistics_fibers)

        # the possible values of s(a) are tau(beta^{-1}(alpha(a)))
        tau_beta_inverse = {}
        self._statistics_possible_values = {}
        for a in self._A:
            v = self._alpha(a)
            if v not in tau_beta_inverse:
                tau_beta_inverse[v] = set(self._tau[b]
                                          for b in self._statistics_fibers[v][1])
            self._statistics_possible_values[a] = tau_beta_inverse[v]

    def statistics_fibers(self):
        r"""
        Return a dictionary mapping statistic values in `W` to their
        preimages in `A` and `B`.

        This is a (computationally) fast way to obtain a first
        impression which objects in `A` should be mapped to which
        objects in `B`.

        EXAMPLES::

            sage: A = B = [permutation for n in range(4) for permutation in Permutations(n)]
            sage: tau = Permutation.longest_increasing_subsequence_length
            sage: wex = lambda p: len(p.weak_excedences())
            sage: fix = lambda p: len(p.fixed_points())
            sage: des = lambda p: len(p.descents(final_descent=True)) if p else 0
            sage: adj = lambda p: len([e for (e, f) in zip(p, p[1:]+[0]) if e == f+1])
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_statistics((len, len), (wex, des), (fix, adj))
            sage: table([[key, AB[0], AB[1]] for key, AB in bij.statistics_fibers().items()])
              (0, 0, 0)   [[]]                                [[]]
              (1, 1, 1)   [[1]]                               [[1]]
              (2, 2, 2)   [[1, 2]]                            [[2, 1]]
              (2, 1, 0)   [[2, 1]]                            [[1, 2]]
              (3, 3, 3)   [[1, 2, 3]]                         [[3, 2, 1]]
              (3, 2, 1)   [[1, 3, 2], [2, 1, 3], [3, 2, 1]]   [[1, 3, 2], [2, 1, 3], [2, 3, 1]]
              (3, 2, 0)   [[2, 3, 1]]                         [[3, 1, 2]]
              (3, 1, 0)   [[3, 1, 2]]                         [[1, 2, 3]]

        """
        return self._statistics_fibers

    def statistics_table(self, header=True):
        r"""
        Provide information about all elements of `A` with corresponding
        `\alpha` values and all elements of `B` with corresponding
        `\beta` and `\tau` values.

        INPUT:

        - ``header`` (optional, default: ``True``) -- whether to
          include a header with the standard greek letters.

        OUTPUT:

        A pair of lists suitable for :class:`~sage.misc.table.table`,
        where

        - the first contains the elements of `A` together with the
          values of `\alpha`

        - the second contains the elements of `B` together with the
          values of `\tau` and `\beta`

        EXAMPLES::

            sage: A = B = [permutation for n in range(4) for permutation in Permutations(n)]
            sage: tau = Permutation.longest_increasing_subsequence_length
            sage: wex = lambda p: len(p.weak_excedences())
            sage: fix = lambda p: len(p.fixed_points())
            sage: des = lambda p: len(p.descents(final_descent=True)) if p else 0
            sage: adj = lambda p: len([e for (e, f) in zip(p, p[1:]+[0]) if e == f+1])
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_statistics((wex, des), (fix, adj))
            sage: a, b = bij.statistics_table()
            sage: table(a, header_row=True, frame=True)
            +-----------+--------+--------+
            | a         | α_1(a) | α_2(a) |
            +===========+========+========+
            | []        | 0      | 0      |
            +-----------+--------+--------+
            | [1]       | 1      | 1      |
            +-----------+--------+--------+
            | [1, 2]    | 2      | 2      |
            +-----------+--------+--------+
            | [2, 1]    | 1      | 0      |
            +-----------+--------+--------+
            | [1, 2, 3] | 3      | 3      |
            +-----------+--------+--------+
            | [1, 3, 2] | 2      | 1      |
            +-----------+--------+--------+
            | [2, 1, 3] | 2      | 1      |
            +-----------+--------+--------+
            | [2, 3, 1] | 2      | 0      |
            +-----------+--------+--------+
            | [3, 1, 2] | 1      | 0      |
            +-----------+--------+--------+
            | [3, 2, 1] | 2      | 1      |
            +-----------+--------+--------+
            sage: table(b, header_row=True, frame=True)
            +-----------+---+--------+--------+
            | b         | τ | β_1(b) | β_2(b) |
            +===========+===+========+========+
            | []        | 0 | 0      | 0      |
            +-----------+---+--------+--------+
            | [1]       | 1 | 1      | 1      |
            +-----------+---+--------+--------+
            | [1, 2]    | 2 | 1      | 0      |
            +-----------+---+--------+--------+
            | [2, 1]    | 1 | 2      | 2      |
            +-----------+---+--------+--------+
            | [1, 2, 3] | 3 | 1      | 0      |
            +-----------+---+--------+--------+
            | [1, 3, 2] | 2 | 2      | 1      |
            +-----------+---+--------+--------+
            | [2, 1, 3] | 2 | 2      | 1      |
            +-----------+---+--------+--------+
            | [2, 3, 1] | 2 | 2      | 1      |
            +-----------+---+--------+--------+
            | [3, 1, 2] | 2 | 2      | 0      |
            +-----------+---+--------+--------+
            | [3, 2, 1] | 1 | 3      | 3      |
            +-----------+---+--------+--------+

        TESTS:

        If no statistics are given, the table should still be able to be generated::

            sage: A = B = [permutation for n in range(3) for permutation in Permutations(n)]
            sage: tau = Permutation.longest_increasing_subsequence_length
            sage: bij = Bijectionist(A, B, tau)
            sage: a, b = bij.statistics_table()
            sage: table(a, header_row=True, frame=True)
            +--------+
            | a      |
            +========+
            | []     |
            +--------+
            | [1]    |
            +--------+
            | [1, 2] |
            +--------+
            | [2, 1] |
            +--------+
            sage: table(b, header_row=True, frame=True)
            +--------+---+
            | b      | τ |
            +========+===+
            | []     | 0 |
            +--------+---+
            | [1]    | 1 |
            +--------+---+
            | [1, 2] | 2 |
            +--------+---+
            | [2, 1] | 1 |
            +--------+---+

        We can omit the header::

            sage: bij.statistics_table(header=True)[1]
            [['b', 'τ'], [[], 0], [[1], 1], [[1, 2], 2], [[2, 1], 1]]
            sage: bij.statistics_table(header=False)[1]
            [[[], 0], [[1], 1], [[1, 2], 2], [[2, 1], 1]]

        """
        # table for alpha
        n_statistics = self._n_statistics
        if header:
            output_alphas = [["a"] + ["\u03b1_" + str(i) + "(a)"
                                      for i in range(1, n_statistics + 1)]]
        else:
            output_alphas = []

        for a in self._A:
            if n_statistics > 0:
                output_alphas.append([a] + list(self._alpha(a)))
            else:
                output_alphas.append([a])

        # table for beta and tau
        if header:
            output_tau_betas = [["b", "\u03c4"] + ["\u03b2_" + str(i) + "(b)"
                                                   for i in range(1, n_statistics + 1)]]
        else:
            output_tau_betas = []
        for b in self._B:
            if n_statistics > 0:
                output_tau_betas.append([b, self._tau[b]] + list(self._beta(b)))
            else:
                output_tau_betas.append([b, self._tau[b]])

        return output_alphas, output_tau_betas

    def set_value_restrictions(self, *a_values):
        r"""
        Restrict the set of possible values `s(a)` for a given element
        `a`.

        .. WARNING::

             Any restriction imposed by a previous invocation of
             :meth:`set_value_restrictions` will be overwritten!

        INPUT:

        - ``a_values`` -- one or more pairs `(a\in A, \tilde
          Z\subseteq Z)`

        EXAMPLES:

        We may want to restrict the value of a given element to a
        single or multiple values.  We do not require that the
        specified values are in the image of `\tau`.  In some
        cases, the restriction may not be able to provide a better
        solution, as for size 3 in the following example. ::

            sage: A = B = [permutation for n in range(4) for permutation in Permutations(n)]
            sage: tau = Permutation.longest_increasing_subsequence_length
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_statistics((len, len))
            sage: bij.set_value_restrictions((Permutation([1, 2]), [1]),
            ....:                            (Permutation([3, 2, 1]), [2, 3, 4]))
            sage: for sol in sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items()))):
            ....:     print(sol)
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 1, [1, 3, 2]: 2, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 3}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 1, [1, 3, 2]: 2, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 3, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 1, [1, 3, 2]: 2, [2, 1, 3]: 2, [2, 3, 1]: 3, [3, 1, 2]: 2, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 1, [1, 3, 2]: 2, [2, 1, 3]: 3, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 1, [1, 3, 2]: 3, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 2, [1, 3, 2]: 1, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 3}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 2, [1, 3, 2]: 1, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 3, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 2, [1, 3, 2]: 1, [2, 1, 3]: 2, [2, 3, 1]: 3, [3, 1, 2]: 2, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 2, [1, 3, 2]: 1, [2, 1, 3]: 3, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 2, [1, 3, 2]: 2, [2, 1, 3]: 1, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 3}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 2, [1, 3, 2]: 2, [2, 1, 3]: 1, [2, 3, 1]: 2, [3, 1, 2]: 3, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 2, [1, 3, 2]: 2, [2, 1, 3]: 1, [2, 3, 1]: 3, [3, 1, 2]: 2, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 2, [1, 3, 2]: 2, [2, 1, 3]: 2, [2, 3, 1]: 1, [3, 1, 2]: 2, [3, 2, 1]: 3}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 2, [1, 3, 2]: 2, [2, 1, 3]: 2, [2, 3, 1]: 1, [3, 1, 2]: 3, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 2, [1, 3, 2]: 2, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 1, [3, 2, 1]: 3}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 2, [1, 3, 2]: 2, [2, 1, 3]: 2, [2, 3, 1]: 3, [3, 1, 2]: 1, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 2, [1, 3, 2]: 2, [2, 1, 3]: 3, [2, 3, 1]: 1, [3, 1, 2]: 2, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 2, [1, 3, 2]: 2, [2, 1, 3]: 3, [2, 3, 1]: 2, [3, 1, 2]: 1, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 2, [1, 3, 2]: 3, [2, 1, 3]: 1, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 2, [1, 3, 2]: 3, [2, 1, 3]: 2, [2, 3, 1]: 1, [3, 1, 2]: 2, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 2, [1, 3, 2]: 3, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 1, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 3, [1, 3, 2]: 1, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 3, [1, 3, 2]: 2, [2, 1, 3]: 1, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 3, [1, 3, 2]: 2, [2, 1, 3]: 2, [2, 3, 1]: 1, [3, 1, 2]: 2, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 3, [1, 3, 2]: 2, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 1, [3, 2, 1]: 2}

        However, an error occurs if the set of possible values is
        empty.  In this example, the image of `\tau` under any
        legal bijection is disjoint to the specified values.

        .. TODO::

            we now have to call
            :meth:`_compute_possible_block_values` for the error
            message. Is this intended behaviour?

            sage: A = B = [permutation for n in range(4) for permutation in Permutations(n)]
            sage: tau = Permutation.longest_increasing_subsequence_length
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_value_restrictions((Permutation([1, 2]), [4, 5]))
            sage: bij._compute_possible_block_values()
            Traceback (most recent call last):
            ...
            ValueError: No possible values found for singleton block [[1, 2]]

        TESTS::

            sage: A = B = [permutation for n in range(4) for permutation in Permutations(n)]
            sage: tau = Permutation.longest_increasing_subsequence_length
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_constant_blocks([[permutation for permutation in Permutations(n)] for n in range(4)])
            sage: bij.set_value_restrictions((Permutation([1, 2]), [4, 5]))
            sage: bij._compute_possible_block_values()
            Traceback (most recent call last):
            ...
            ValueError: No possible values found for block [[1, 2], [2, 1]]

            sage: A = B = [permutation for n in range(4) for permutation in Permutations(n)]
            sage: tau = Permutation.longest_increasing_subsequence_length
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_value_restrictions(((1, 2), [4, 5, 6]))
            Traceback (most recent call last):
            ...
            AssertionError: Element (1, 2) was not found in A

        """
        # reset values
        set_Z = set(self._Z)
        self._restrictions_possible_values = {a: set_Z for a in self._A}
        for a, values in a_values:
            assert a in self._A, f"Element {a} was not found in A"
            self._restrictions_possible_values[a] = self._restrictions_possible_values[a].intersection(values)

    def _compute_possible_block_values(self):
        r"""
        Update the dictionary of possible values of each block.
        """
        self._possible_block_values = {}  # P -> Power(Z)
        for p, block in self._P.root_to_elements_dict().items():
            self._possible_block_values[p] = set.intersection(*[self._restrictions_possible_values[a] for a in block],
                                                              *[self._statistics_possible_values[a] for a in block])
            if not self._possible_block_values[p]:
                if len(block) == 1:
                    raise ValueError(f"No possible values found for singleton block {block}")
                else:
                    raise ValueError(f"No possible values found for block {block}")

    def set_distributions(self, *elements_distributions):
        r"""
        Specify the distribution of `s` for a subset of elements.

        .. WARNING::

             Any restriction imposed by a previous invocation of
             :meth:`set_distributions` will be overwritten!

        INPUT:

        - one or more pairs of `(\tilde A\subseteq A, \tilde Z)`,
          where `\tilde Z` is a list of values in `Z` of the same
          size as `\tilde A`

        This method specifies that `\{s(a) | a\in\tilde A\}` equals
        ``\tilde Z`` as a multiset for each of the pairs.

        When specifying several distributions, the subsets of `A` do
        not have to be disjoint.

        EXAMPLES::

            sage: A = B = [permutation for n in range(4) for permutation in Permutations(n)]
            sage: tau = Permutation.longest_increasing_subsequence_length
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_statistics((len, len))
            sage: bij.set_distributions(([Permutation([1, 2, 3]), Permutation([1, 3, 2])], [1, 3]))
            sage: for sol in sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items()))):
            ....:     print(sol)
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 1, [1, 3, 2]: 3, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 3, [1, 3, 2]: 1, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 1, [1, 2, 3]: 1, [1, 3, 2]: 3, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 1, [1, 2, 3]: 3, [1, 3, 2]: 1, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 2}

            sage: bij.constant_blocks(optimal=True)
            {{[2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]}}
            sage: sorted(bij.minimal_subdistributions_blocks_iterator(), key=lambda d: (len(d[0]), d[0]))
            [([[]], [0]),
             ([[1]], [1]),
             ([[2, 1, 3]], [2]),
             ([[1, 2], [2, 1]], [1, 2]),
             ([[1, 2, 3], [1, 3, 2]], [1, 3])]

        We may also specify multiple, possibly overlapping distributions::

            sage: bij.set_distributions(([Permutation([1, 2, 3]), Permutation([1, 3, 2])], [1, 3]),
            ....:                       ([Permutation([1, 3, 2]), Permutation([3, 2, 1]),
            ....:                        Permutation([2, 1, 3])], [1, 2, 2]))
            sage: for sol in sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items()))):
            ....:     print(sol)
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 3, [1, 3, 2]: 1, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 1, [1, 2, 3]: 3, [1, 3, 2]: 1, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 2}

            sage: bij.constant_blocks(optimal=True)
            {{[1], [1, 3, 2]}, {[2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1]}}
            sage: sorted(bij.minimal_subdistributions_blocks_iterator(), key=lambda d: (len(d[0]), d[0]))
            [([[]], [0]),
             ([[1]], [1]),
             ([[1, 2, 3]], [3]),
             ([[2, 1, 3]], [2]),
             ([[1, 2], [2, 1]], [1, 2])]

        TESTS:

        Because of the current implementation of the output
        calculation, we do not improve our solution if we do not gain
        any unique solutions.::

            sage: A = B = [permutation for n in range(4) for permutation in Permutations(n)]
            sage: tau = Permutation.longest_increasing_subsequence_length
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_statistics((len, len))
            sage: bij.set_distributions(([Permutation([1, 2, 3]), Permutation([1, 3, 2])], [2, 3]))
            sage: for sol in sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items()))):
            ....:     print(sol)
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 2, [1, 3, 2]: 3, [2, 1, 3]: 1, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 2, [1, 3, 2]: 3, [2, 1, 3]: 2, [2, 3, 1]: 1, [3, 1, 2]: 2, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 2, [1, 3, 2]: 3, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 1, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 2, [1, 3, 2]: 3, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 1}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 3, [1, 3, 2]: 2, [2, 1, 3]: 1, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 3, [1, 3, 2]: 2, [2, 1, 3]: 2, [2, 3, 1]: 1, [3, 1, 2]: 2, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 3, [1, 3, 2]: 2, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 1, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 3, [1, 3, 2]: 2, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 1}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 1, [1, 2, 3]: 2, [1, 3, 2]: 3, [2, 1, 3]: 1, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 1, [1, 2, 3]: 2, [1, 3, 2]: 3, [2, 1, 3]: 2, [2, 3, 1]: 1, [3, 1, 2]: 2, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 1, [1, 2, 3]: 2, [1, 3, 2]: 3, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 1, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 1, [1, 2, 3]: 2, [1, 3, 2]: 3, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 1}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 1, [1, 2, 3]: 3, [1, 3, 2]: 2, [2, 1, 3]: 1, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 1, [1, 2, 3]: 3, [1, 3, 2]: 2, [2, 1, 3]: 2, [2, 3, 1]: 1, [3, 1, 2]: 2, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 1, [1, 2, 3]: 3, [1, 3, 2]: 2, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 1, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 1, [1, 2, 3]: 3, [1, 3, 2]: 2, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 1}

        Another example with statistics::

            sage: bij = Bijectionist(A, B, tau)
            sage: alpha = lambda p: p(1) if len(p) > 0 else 0
            sage: beta = lambda p: p(1) if len(p) > 0 else 0
            sage: bij.set_statistics((alpha, beta), (len, len))
            sage: for sol in sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items()))):
            ....:     print(sol)
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 1, [1, 2, 3]: 2, [1, 3, 2]: 3, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 1, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 1, [1, 2, 3]: 2, [1, 3, 2]: 3, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 1}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 1, [1, 2, 3]: 3, [1, 3, 2]: 2, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 1, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 1, [1, 2, 3]: 3, [1, 3, 2]: 2, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 1}

        The solution above is not unique.  We can add a feasible distribution to force uniqueness::

            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_statistics((alpha, beta), (len, len))
            sage: bij.set_distributions(([Permutation([1, 2, 3]), Permutation([3, 2, 1])], [1, 3]))
            sage: for sol in bij.solutions_iterator():
            ....:     print(sol)
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 1, [1, 2, 3]: 3, [1, 3, 2]: 2, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 1}

        Let us try to add a distribution that cannot be satisfied,
        because there is no solution where a permutation that starts
        with 1 is mapped onto 1::

            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_statistics((alpha, beta), (len, len))
            sage: bij.set_distributions(([Permutation([1, 2, 3]), Permutation([1, 3, 2])], [1, 3]))
            sage: list(bij.solutions_iterator())
            []

        The specified elements have to be in `A` and have to be of the same size::

            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_statistics((len, len))
            sage: bij.set_distributions(([Permutation([1, 2, 3, 4])], [1]))
            Traceback (most recent call last):
            ...
            ValueError: Element [1, 2, 3, 4] was not found in A!
            sage: bij.set_distributions(([Permutation([1, 2, 3])], [-1]))
            Traceback (most recent call last):
            ...
            ValueError: Value -1 was not found in tau(A)!

        Note that the same error occurs when an element that is not the first element of the list is
        not in `A`.

        """
        for elements, values in elements_distributions:
            assert len(elements) == len(values), f"{elements} and {values} are not of the same size!"
            for a, z in zip(elements, values):
                if a not in self._A:
                    raise ValueError(f"Element {a} was not found in A!")
                if z not in self._Z:
                    raise ValueError(f"Value {z} was not found in tau(A)!")
        self._elements_distributions = elements_distributions

    def set_intertwining_relations(self, *pi_rho):
        r"""
        Add restrictions of the form `s(\pi(a_1,\dots, a_k)) =
        \rho(s(a_1),\dots, s(a_k))`.

        .. WARNING::

             Any restriction imposed by a previous invocation of
             :meth:`set_intertwining_relations` will be overwritten!

        INPUT:

        - ``pi_rho`` -- one or more tuples `(k, \pi: A^k\to A, \rho:
          Z^k\to Z, \tilde A)` where `\tilde A` (optional) is an
          `k`-ary function that returns true if and only if an
          `k`-tuple of objects in `A` is in the domain of `\pi`

        EXAMPLES:

        We can concatenate two permutations, by increasing the values
        of the second permutation by the length of the first
        permutation::

            sage: concat = lambda p1, p2: Permutation(p1 + [i + len(p1) for i in p2])

        We may be interested in statistics on permutations which are
        equidistributed with the number of fixed points, such that
        concatenating permutations corresponds to adding statistic
        values::

            sage: A = B = [permutation for n in range(4) for permutation in Permutations(n)]
            sage: bij = Bijectionist(A, B, Permutation.number_of_fixed_points)
            sage: bij.set_statistics((len, len))
            sage: for solution in sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items()))):
            ....:     print(solution)
            ...
            {[]: 0, [1]: 1, [1, 2]: 0, [2, 1]: 2, [1, 2, 3]: 1, [1, 3, 2]: 0, [2, 1, 3]: 0, [2, 3, 1]: 1, [3, 1, 2]: 1, [3, 2, 1]: 3}
            ...
            {[]: 0, [1]: 1, [1, 2]: 0, [2, 1]: 2, [1, 2, 3]: 1, [1, 3, 2]: 1, [2, 1, 3]: 3, [2, 3, 1]: 0, [3, 1, 2]: 0, [3, 2, 1]: 1}
            ...

            sage: bij.set_intertwining_relations((2, concat, lambda x, y: x + y))
            sage: for solution in sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items()))):
            ....:     print(solution)
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 0, [1, 2, 3]: 3, [1, 3, 2]: 1, [2, 1, 3]: 1, [2, 3, 1]: 0, [3, 1, 2]: 0, [3, 2, 1]: 1}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 0, [1, 2, 3]: 3, [1, 3, 2]: 1, [2, 1, 3]: 1, [2, 3, 1]: 0, [3, 1, 2]: 1, [3, 2, 1]: 0}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 0, [1, 2, 3]: 3, [1, 3, 2]: 1, [2, 1, 3]: 1, [2, 3, 1]: 1, [3, 1, 2]: 0, [3, 2, 1]: 0}

        The domain of the composition may be restricted.  E.g., if we
        concatenate only permutations starting with a 1, we obtain
        fewer forced elements::

            sage: in_domain = lambda p1, p2: (not p1 or p1(1) == 1) and (not p2 or p2(1) == 1)
            sage: bij.set_intertwining_relations((2, concat, lambda x, y: x + y, in_domain))
            sage: for solution in sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items()))):
            ....:     print(solution)
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 0, [1, 2, 3]: 3, [1, 3, 2]: 0, [2, 1, 3]: 0, [2, 3, 1]: 1, [3, 1, 2]: 1, [3, 2, 1]: 1}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 0, [1, 2, 3]: 3, [1, 3, 2]: 0, [2, 1, 3]: 1, [2, 3, 1]: 0, [3, 1, 2]: 1, [3, 2, 1]: 1}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 0, [1, 2, 3]: 3, [1, 3, 2]: 0, [2, 1, 3]: 1, [2, 3, 1]: 1, [3, 1, 2]: 0, [3, 2, 1]: 1}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 0, [1, 2, 3]: 3, [1, 3, 2]: 0, [2, 1, 3]: 1, [2, 3, 1]: 1, [3, 1, 2]: 1, [3, 2, 1]: 0}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 0, [1, 2, 3]: 3, [1, 3, 2]: 1, [2, 1, 3]: 0, [2, 3, 1]: 0, [3, 1, 2]: 1, [3, 2, 1]: 1}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 0, [1, 2, 3]: 3, [1, 3, 2]: 1, [2, 1, 3]: 0, [2, 3, 1]: 1, [3, 1, 2]: 0, [3, 2, 1]: 1}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 0, [1, 2, 3]: 3, [1, 3, 2]: 1, [2, 1, 3]: 0, [2, 3, 1]: 1, [3, 1, 2]: 1, [3, 2, 1]: 0}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 0, [1, 2, 3]: 3, [1, 3, 2]: 1, [2, 1, 3]: 1, [2, 3, 1]: 0, [3, 1, 2]: 0, [3, 2, 1]: 1}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 0, [1, 2, 3]: 3, [1, 3, 2]: 1, [2, 1, 3]: 1, [2, 3, 1]: 0, [3, 1, 2]: 1, [3, 2, 1]: 0}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 0, [1, 2, 3]: 3, [1, 3, 2]: 1, [2, 1, 3]: 1, [2, 3, 1]: 1, [3, 1, 2]: 0, [3, 2, 1]: 0}

        We can also restrict according to several composition
        functions.  For example, we may additionally concatenate
        permutations by incrementing the elements of the first::

            sage: skew_concat = lambda p1, p2: Permutation([i + len(p2) for i in p1] + list(p2))
            sage: bij.set_intertwining_relations((2, skew_concat, lambda x, y: x + y))
            sage: for solution in sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items()))):
            ....:     print(solution)
            {[]: 0, [1]: 1, [1, 2]: 0, [2, 1]: 2, [1, 2, 3]: 0, [1, 3, 2]: 0, [2, 1, 3]: 1, [2, 3, 1]: 1, [3, 1, 2]: 1, [3, 2, 1]: 3}
            {[]: 0, [1]: 1, [1, 2]: 0, [2, 1]: 2, [1, 2, 3]: 0, [1, 3, 2]: 1, [2, 1, 3]: 0, [2, 3, 1]: 1, [3, 1, 2]: 1, [3, 2, 1]: 3}
            {[]: 0, [1]: 1, [1, 2]: 0, [2, 1]: 2, [1, 2, 3]: 1, [1, 3, 2]: 0, [2, 1, 3]: 0, [2, 3, 1]: 1, [3, 1, 2]: 1, [3, 2, 1]: 3}

        However, this yields no solution::

            sage: bij.set_intertwining_relations((2, concat, lambda x, y: x + y), (2, skew_concat, lambda x, y: x + y))
            sage: list(bij.solutions_iterator())
            []

        """
        Pi_Rho = namedtuple("Pi_Rho", "numargs pi rho domain")
        self._pi_rho = []

        for pi_rho_tuple in pi_rho:
            if len(pi_rho_tuple) == 3:
                k, pi, rho = pi_rho_tuple
                domain = None
            else:
                k, pi, rho, domain = pi_rho_tuple

            self._pi_rho.append(Pi_Rho(numargs=k, pi=pi, rho=rho, domain=domain))

    def _forced_constant_blocks(self):
        r"""
        Modify current partition into blocks to the coarsest possible
        one, meaning that after calling this function for every two
        distinct blocks `p_1`, `p_2` there exists a solution `s` with
        `s(p_1)\neq s(p_2)`.

        ALGORITHM:

        First we generate an initial solution.  For all blocks i, j
        that have the same value under this initial solution, we add
        the constraint `x[i, z] + x[j, z] <= 1` for all possible
        values `z\in Z`.  This constraint ensures that the `s` differs
        on the two blocks.  If this modified problem does not have a
        solution, we know that the two blocks always have the same
        value and join them.  Then we save all values of this new
        solution and continue looking at pairs of blocks that had the
        same value under all calculated solutions, until no blocks
        can be joined anymore.

        EXAMPLES:

        The easiest example is given by a constant `tau`, so everything
        is forced to be the same value:

            sage: A = B = [permutation for n in range(3) for permutation in Permutations(n)]
            sage: bij = Bijectionist(A, B, lambda x: 0)
            sage: bij.constant_blocks()
            {}
            sage: bij.constant_blocks(optimal=True)
            {{[], [1], [1, 2], [2, 1]}}

        In this other example we look at permutations with length 2 and 3::

            sage: N = 4
            sage: A = B = [permutation for n in range(2, N) for permutation in Permutations(n)]
            sage: tau = lambda p: p[0] if len(p) else 0
            sage: add_n = lambda p1: Permutation(p1 + [1 + len(p1)])
            sage: add_1 = lambda p1: Permutation([1] + [1 + i for i in p1])
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_intertwining_relations((1, add_n, lambda x: x + 1), (1, add_1, lambda x: x + 1))
            sage: bij.set_statistics((len, len))

            sage: bij.constant_blocks()
            {}
            sage: bij.constant_blocks(optimal=True)
            {{[1, 3, 2], [2, 1, 3]}}

        Indeed, ``[1,3,2]`` and ``[2,1,3]`` have the same value in
        all solutions, but different values are possible::

            sage: pi1 = Permutation([1,3,2]); pi2 = Permutation([2,1,3]);
            sage: set([(solution[pi1], solution[pi2]) for solution in bij.solutions_iterator()])
            {(2, 2), (3, 3)}

        Another example involving the cycle type of permutations::

            sage: A = B = [permutation for n in range(4) for permutation in Permutations(n)]
            sage: bij = Bijectionist(A, B, lambda x: x.cycle_type())

        Let us require that each permutation has the same value as its inverse::

            sage: from sage.combinat.cyclic_sieving_phenomenon import orbit_decomposition
            sage: P = orbit_decomposition([permutation for n in range(4) for permutation in Permutations(n)], Permutation.inverse)
            sage: bij.set_constant_blocks(P)
            sage: bij.constant_blocks()
            {{[2, 3, 1], [3, 1, 2]}}

            sage: concat = lambda p1, p2: Permutation(p1 + [i + len(p1) for i in p2])
            sage: union = lambda p1, p2: Partition(sorted(list(p1) + list(p2), reverse=True))
            sage: bij.set_intertwining_relations((2, concat, union))

        In this case we do not discover constant blocks by looking at the intertwining_relations only::

            sage: next(bij.solutions_iterator())
            ...
            sage: bij.constant_blocks()
            {{[2, 3, 1], [3, 1, 2]}}

            sage: bij.constant_blocks(optimal=True)
            {{[1, 3, 2], [2, 1, 3], [3, 2, 1]}, {[2, 3, 1], [3, 1, 2]}}

        TESTS::

            sage: N = 4
            sage: A = B = [permutation for n in range(N + 1) for permutation in Permutations(n)]
            sage: alpha1 = lambda p: len(p.weak_excedences())
            sage: alpha2 = lambda p: len(p.fixed_points())
            sage: beta1 = lambda p: len(p.descents(final_descent=True)) if p else 0
            sage: beta2 = lambda p: len([e for (e, f) in zip(p, p[1:] + [0]) if e == f + 1])
            sage: tau = Permutation.longest_increasing_subsequence_length
            sage: def rotate_permutation(p):
            ....:    cycle = Permutation(tuple(range(1, len(p) + 1)))
            ....:    return Permutation([cycle.inverse()(p(cycle(i))) for i in range(1, len(p) + 1)])
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_statistics((alpha1, beta1), (alpha2, beta2))
            sage: from sage.combinat.cyclic_sieving_phenomenon import orbit_decomposition
            sage: bij.set_constant_blocks(orbit_decomposition(A, rotate_permutation))
            sage: for p in bij.constant_blocks(): print(list(p))
            [[2, 1, 3, 4], [1, 2, 4, 3], [1, 3, 2, 4], [4, 2, 3, 1]]
            [[3, 2, 1], [1, 3, 2], [2, 1, 3]]
            [[2, 4, 3, 1], [3, 2, 4, 1], [2, 3, 1, 4], [1, 3, 4, 2]]
            [[1, 4, 2, 3], [3, 1, 2, 4], [4, 2, 1, 3], [4, 1, 3, 2]]
            [[1, 4, 3, 2], [3, 2, 1, 4]]
            [[2, 1, 4, 3], [4, 3, 2, 1]]
            [[2, 4, 1, 3], [3, 4, 2, 1], [4, 3, 1, 2], [3, 1, 4, 2]]

            sage: for p in bij.constant_blocks(optimal=True): sorted(p, key=len)
            [[1], [1, 2], [1, 2, 3], [1, 2, 3, 4]]
            [[1, 3, 2],
             [2, 1, 3],
             [3, 2, 1],
             [2, 3, 4, 1],
             [1, 3, 4, 2],
             [2, 1, 3, 4],
             [1, 3, 2, 4],
             [2, 3, 1, 4],
             [1, 2, 4, 3],
             [3, 2, 4, 1],
             [2, 1, 4, 3],
             [2, 4, 3, 1],
             [4, 2, 3, 1],
             [4, 3, 2, 1],
             [1, 4, 3, 2],
             [3, 2, 1, 4]]
            [[1, 4, 2, 3],
             [4, 2, 1, 3],
             [2, 4, 1, 3],
             [4, 3, 1, 2],
             [4, 1, 3, 2],
             [3, 4, 2, 1],
             [3, 1, 2, 4],
             [3, 1, 4, 2]]

        The permutation `[2, 1]` is in none of these blocks::

            sage: bij.set_constant_blocks(orbit_decomposition(A, rotate_permutation))
            sage: all(s[Permutation([2, 1])] == s[Permutation([1])] for s in bij.solutions_iterator())
            False

            sage: all(s[Permutation([2, 1])] == s[Permutation([1, 3, 2])] for s in bij.solutions_iterator())
            False

            sage: all(s[Permutation([2, 1])] == s[Permutation([1, 4, 2, 3])] for s in bij.solutions_iterator())
            False



            sage: A = B = ["a", "b", "c", "d", "e", "f"]
            sage: tau = {"a": 1, "b": 1, "c": 3, "d": 4, "e": 5, "f": 6}.get
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_distributions((["a", "b"], [1, 1]), (["c", "d", "e"], [3, 4, 5]))
            sage: bij.constant_blocks()
            {}
            sage: bij.constant_blocks(optimal=True)
            {{'a', 'b'}}

            sage: A = B = ["a", "b", "c", "d", "e", "f"]
            sage: tau = {"a": 1, "b": 1, "c": 5, "d": 4, "e": 4, "f": 6}.get
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_distributions((["a", "b"], [1, 1]), (["d", "e"], [4, 4]))
            sage: bij.constant_blocks(optimal=True)
            {{'a', 'b'}, {'d', 'e'}}

            sage: A = B = ["a", "b", "c", "d"]
            sage: tau = {"a": 1, "b": 1, "c": 2, "d": 2}.get
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.constant_blocks(optimal=True)
            {}
            sage: bij.set_constant_blocks([["a", "b"]])
            sage: bij.constant_blocks()
            {{'a', 'b'}}
            sage: bij.constant_blocks(optimal=True)
            {{'a', 'b'}, {'c', 'd'}}

        """
        bmilp = self._generate_and_solve_initial_bmilp()  # may throw Exception

        # generate blockwise preimage to determine which blocks have the same image
        solution = self._solution_by_blocks(bmilp)
        multiple_preimages = {(value,): preimages
                              for value, preimages in _invert_dict(solution).items()
                              if len(preimages) > 1}

        # check for each pair of blocks if a solution with different values on these block exists
        #     if yes, use the new solution to update the multiple_preimages dictionary, restart the check
        #     if no, the two blocks can be joined

        # _P has to be copied to not mess with the solution-process
        # since we do not want to regenerate the bmilp in each step, so blocks
        # have to stay consistent during the whole process
        tmp_P = deepcopy(self._P)
        updated_preimages = True
        while updated_preimages:
            updated_preimages = False
            for values in copy(multiple_preimages):  # copy to be able to modify dict
                if updated_preimages:
                    break
                for i, j in itertools.combinations(copy(multiple_preimages[values]), r=2):  # copy to be able to modify list
                    bmilp_veto = deepcopy(bmilp)  # adding constraints to a simple copy adds them to the original instance, too
                    try:
                        # veto the two blocks having the same value
                        for z in self._possible_block_values[i]:
                            if z in self._possible_block_values[j]:  # intersection
                                bmilp_veto.milp.add_constraint(bmilp_veto._x[i, z] + bmilp_veto._x[j, z] <= 1)
                        bmilp_veto.milp.solve()

                        # solution exists, update dictionary
                        solution = self._solution_by_blocks(bmilp_veto)
                        updated_multiple_preimages = {}
                        for values in multiple_preimages:
                            for p in multiple_preimages[values]:
                                solution_tuple = (*values, solution[p])  # tuple so actual solutions were equal in lookup
                                if solution_tuple not in updated_multiple_preimages:
                                    updated_multiple_preimages[solution_tuple] = []
                                updated_multiple_preimages[solution_tuple].append(p)
                        updated_preimages = True
                        multiple_preimages = updated_multiple_preimages
                        break
                    except MIPSolverException:
                        # no solution exists, join blocks
                        tmp_P.union(i, j)
                        if i in multiple_preimages[values] and j in multiple_preimages[values]:  # only one of the joined blocks should remain in the list
                            multiple_preimages[values].remove(j)
                    if len(multiple_preimages[values]) == 1:
                        del multiple_preimages[values]
                        break

        self.set_constant_blocks(tmp_P)

    def possible_values(self, p=None, optimal=False):
        r"""Return for each block the values of `s` compatible with the
        imposed restrictions.

        .. TODO::

            should this method update and return
            ``self._possible_block_values``?

        INPUT:

        - ``p`` (optional, default: ``None``) -- a block of `P`, or
          an element of a block of `P`, or a list of these

        - ``optimal`` (optional, default: ``False``) -- whether or
          not to compute the minimal possible set of statistic values,
          throws a MIPSolverException if no solution is found.

        .. NOTE::

            computing the minimal possible set of statistic values
            may be computationally expensive.

        TESTS::

            sage: A = B = ["a", "b", "c", "d"]
            sage: tau = {"a": 1, "b": 1, "c": 2, "d": 2}.get
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_constant_blocks([["a", "b"]])

        Test if all formats are really possible::

            sage: bij.possible_values(p="a")
            {'a': {1, 2}, 'b': {1, 2}}
            sage: bij.possible_values(p=["a", "b"])
            {'a': {1, 2}, 'b': {1, 2}}
            sage: bij.possible_values(p=[["a", "b"]])
            {'a': {1, 2}, 'b': {1, 2}}
            sage: bij.possible_values(p=[["a", "b"], ["c"]])
            {'a': {1, 2}, 'b': {1, 2}, 'c': {1, 2}}

        Test optimal::

            sage: bij.possible_values(p=["a", "c"], optimal=True)
            {'a': {1, 2}, 'b': {1, 2}, 'c': {1, 2}}

        Verify by listing all solutions::

            sage: sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items())))
            [{'a': 1, 'b': 1, 'c': 2, 'd': 2}, {'a': 2, 'b': 2, 'c': 1, 'd': 1}]

        Test if MIPSolverException is thrown::

            sage: A = B = list('ab')
            sage: bij = Bijectionist(A, B, lambda x: B.index(x) % 2)
            sage: bij.set_constant_blocks([['a', 'b']])
            sage: bij.possible_values(p="a")
            {'a': {0, 1}, 'b': {0, 1}}
            sage: bij.possible_values(p="a", optimal=True)
            Traceback (most recent call last):
            ...
            sage.numerical.mip.MIPSolverException: ...

        Another example::

            sage: A = B = [permutation for n in range(4) for permutation in Permutations(n)]
            sage: tau = Permutation.longest_increasing_subsequence_length
            sage: bij = Bijectionist(A, B, tau)
            sage: alpha = lambda p: p(1) if len(p) > 0 else 0
            sage: beta = lambda p: p(1) if len(p) > 0 else 0
            sage: bij.set_statistics((alpha, beta), (len, len))
            sage: for sol in sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items()))):
            ....:     print(sol)
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 1, [1, 2, 3]: 2, [1, 3, 2]: 3, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 1, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 1, [1, 2, 3]: 2, [1, 3, 2]: 3, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 1}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 1, [1, 2, 3]: 3, [1, 3, 2]: 2, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 1, [3, 2, 1]: 2}
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 1, [1, 2, 3]: 3, [1, 3, 2]: 2, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 1}
            sage: bij.possible_values(p=[Permutation([1]), Permutation([1, 2, 3]), Permutation([3, 1, 2])], optimal=True)
            {[1]: {1}, [1, 2, 3]: {2, 3}, [3, 1, 2]: {1, 2}}

        Another example::

            sage: N = 2; A = B = [dyck_word for n in range(N+1) for dyck_word in DyckWords(n)]
            sage: tau = lambda D: D.number_of_touch_points()
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_statistics((lambda d: d.semilength(), lambda d: d.semilength()))
            sage: for solution in sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items()))):
            ....:     print(solution)
            {[]: 0, [1, 0]: 1, [1, 0, 1, 0]: 1, [1, 1, 0, 0]: 2}
            {[]: 0, [1, 0]: 1, [1, 0, 1, 0]: 2, [1, 1, 0, 0]: 1}
            sage: bij.possible_values(p=[DyckWord([]), DyckWord([1, 0]), DyckWord([1, 0, 1, 0]), DyckWord([1, 1, 0, 0])], optimal=True)
            {[]: {0}, [1, 0]: {1}, [1, 0, 1, 0]: {1, 2}, [1, 1, 0, 0]: {1, 2}}

        .. TODO::

            test to show that the solution for all blocks is not more
            expensive than using :meth:`solutions_iterator`

        """
        # convert input to set of block representatives
        blocks = set()
        if p in self._A:
            blocks.add(self._P.find(p))
        elif type(p) is list:
            for p1 in p:
                if p1 in self._A:
                    blocks.add(self._P.find(p1))
                elif type(p1) is list:
                    for p2 in p1:
                        blocks.add(self._P.find(p2))

        if optimal:
            # function adding a solution to dict of solutions
            def add_solution(solutions, solution):
                for p, value in solution.items():
                    if p not in solutions:
                        solutions[p] = set()
                    solutions[p].add(value)

            # generate initial solution, solution dict and add solution
            bmilp = self._generate_and_solve_initial_bmilp()
            solution = self._solution(bmilp)
            solutions = {}
            add_solution(solutions, solution)

            # iterate through blocks and generate all values
            for p in blocks:
                veto_bmilp = deepcopy(bmilp)  # adding constraints to a simple copy adds them to the original instance, too
                for value in solutions[p]:
                    veto_bmilp.milp.add_constraint(veto_bmilp._x[p, value] == 0)
                while True:
                    try:
                        veto_bmilp.milp.solve()
                        # problem has a solution, so new value was found
                        solution = self._solution(veto_bmilp)
                        add_solution(solutions, solution)
                        # veto new value and try again
                        veto_bmilp.milp.add_constraint(veto_bmilp._x[p, solution[p]] == 0)
                    except MIPSolverException:
                        # no solution, so all possible values have been found
                        break

            # TODO: update possible block values if wanted

        # create dictionary to return
        possible_values = {}
        for p in blocks:
            for a in self._P.root_to_elements_dict()[p]:  # TODO: is this the format we want to return in or possible_values[block]?
                if optimal:
                    possible_values[a] = solutions[p]
                else:
                    possible_values[a] = self._possible_block_values[p]

        return possible_values

    def minimal_subdistributions_iterator(self, tA=None):
        r"""
        Return all minimal subsets `\tilde A` of `A` containing `tA`
        together with submultisets `\tilde Z` with `s(\tilde A) =
        \tilde Z` as multisets.

        .. TODO::

           should this method interact with ``self._elements_distributions``?

        INPUT:

        - ``tA`` (optional, default: ``None``) -- a subset of `A` TODO: add this

        If ``tA`` is not ``None``, return an iterator of the
        subdistributions containing ``tA``.

        TESTS::

            sage: A = B = [permutation for n in range(3) for permutation in Permutations(n)]
            sage: bij = Bijectionist(A, B, len)
            sage: bij.set_statistics((len, len))
            sage: for sol in bij.solutions_iterator():
            ....:     print(sol)
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 2}
            sage: sorted(bij.minimal_subdistributions_iterator())
            [([[]], [0]), ([[1]], [1]), ([[1, 2]], [2]), ([[2, 1]], [2])]

        Another example::

            sage: N = 2; A = B = [dyck_word for n in range(N+1) for dyck_word in DyckWords(n)]
            sage: tau = lambda D: D.number_of_touch_points()
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_statistics((lambda d: d.semilength(), lambda d: d.semilength()))
            sage: for solution in sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items()))):
            ....:     print(solution)
            {[]: 0, [1, 0]: 1, [1, 0, 1, 0]: 1, [1, 1, 0, 0]: 2}
            {[]: 0, [1, 0]: 1, [1, 0, 1, 0]: 2, [1, 1, 0, 0]: 1}
            sage: for subdistribution in bij.minimal_subdistributions_iterator():
            ....:     print(subdistribution)
            ([[]], [0])
            ([[1, 0]], [1])
            ([[1, 0, 1, 0], [1, 1, 0, 0]], [1, 2])

        An example with two elements of the same block in a subdistribution::

            sage: A = B = ["a", "b", "c", "d", "e"]
            sage: tau = {"a": 1, "b": 1, "c": 2, "d": 2, "e": 3}.get
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_constant_blocks([["a", "b"]])
            sage: bij.set_value_restrictions(("a", [1, 2]))
            sage: bij.constant_blocks(optimal=True)
            {{'a', 'b'}}
            sage: list(bij.minimal_subdistributions_iterator())
            [(['a', 'b', 'c', 'd', 'e'], [1, 1, 2, 2, 3])]
        """
        # see
        # https://mathoverflow.net/questions/406751/find-a-subdistribution/406975
        # and
        # https://gitlab.com/mantepse/bijection-tools/-/issues/29

        minimal_subdistribution = MixedIntegerLinearProgram(maximization=False, solver=self._solver)
        D = minimal_subdistribution.new_variable(binary=True)  # the subset of elements
        V = minimal_subdistribution.new_variable(integer=True)  # the subdistribution
        minimal_subdistribution.set_objective(sum(D[a] for a in self._A))
        minimal_subdistribution.add_constraint(sum(D[a] for a in self._A) >= 1)

        try:
            bmilp = self._generate_and_solve_initial_bmilp()
        except MIPSolverException:
            return
        s = self._solution(bmilp)
        while True:
            for v in self._Z:
                minimal_subdistribution.add_constraint(sum(D[a] for a in self._A if s[a] == v) == V[v])
            try:
                minimal_subdistribution.solve()
            except MIPSolverException:
                return
            d = minimal_subdistribution.get_values(D)  # a dict from A to {0, 1}
            new_s = self._find_counter_example(bmilp, s, d)
            if new_s is None:
                values = self._sorter["Z"](s[a] for a in self._A if d[a])
                yield ([a for a in self._A if d[a]], values)

                # get all variables with value 1
                active_vars = [D[a] for a in self._A
                               if minimal_subdistribution.get_values(D[a])]

                # add constraint that not all of these can be 1, thus vetoing
                # the current solution
                minimal_subdistribution.add_constraint(sum(active_vars) <= len(active_vars) - 1,
                                                       name="veto")
                # TODO: can we ignore that in the next step the same constraint is added again?
            else:
                s = new_s

    def _find_counter_example(self, bmilp, s0, d):
        r"""
        Return a solution `s` such that ``d`` is not a subdistribution of
        `s0`.

        TODO: better name

        INPUT:

        - ``bmilp``, the mixed linear integer program

        - ``s0``, a solution

        - ``d``, a subset of `A`, in the form of a dict from `A` to `\{0, 1\}`
        """
        for v in self._Z:
            v_in_d_count = sum(d[a] for a in self._A if s0[a] == v)
            if not v_in_d_count:
                continue

            veto_bmilp = deepcopy(bmilp)  # adding constraints to a simple copy adds them to the original instance, too
            # try to find a solution which has a different
            # subdistribution on d than s0
            v_in_d = sum(d[a] * veto_bmilp._x[self._P.find(a), v]
                         for a in self._A
                         if v in self._possible_block_values[self._P.find(a)])

            # it is sufficient to require that v occurs less often as
            # a value among {a | d[a] == 1} than it does in
            # v_in_d_count, because, if the distributions are
            # different, one such v must exist
            veto_bmilp.milp.add_constraint(v_in_d <= v_in_d_count - 1)
            try:
                veto_bmilp.milp.solve()
                return self._solution(veto_bmilp)
            except MIPSolverException:
                pass
        return

    def minimal_subdistributions_blocks_iterator(self, p=None):
        r"""Return all representatives of minimal subsets `\tilde P`
        of `P` containing `p` together with submultisets `\tilde Z`
        with `s(\tilde P) = \tilde Z` as multisets.

        .. WARNING::

            If there are several solutions with the same support
            (i.e., the sets of block representatives are the same),
            only one of these will be found, even if the
            distributions are different, see the doctest below.  To
            find all solutions, use
            :meth:`minimal_subdistributions_iterator`, which is,
            however, computationally more expensive.

        .. TODO::

            should this method interact with ``self._elements_distributions``?

        INPUT:

        - ``p`` (optional, default: ``None``) -- a subset of `P`

        If ``p`` is not ``None``, return an iterator of the
        subdistributions containing ``p``.

        .. TODO::

           the optional argument is not yet supported

        EXAMPLES::

            sage: A = B = [permutation for n in range(3) for permutation in Permutations(n)]
            sage: bij = Bijectionist(A, B, len)
            sage: bij.set_statistics((len, len))
            sage: for sol in bij.solutions_iterator():
            ....:     print(sol)
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 2}
            sage: sorted(bij.minimal_subdistributions_blocks_iterator())
            [([[]], [0]), ([[1]], [1]), ([[1, 2]], [2]), ([[2, 1]], [2])]

        Another example::

            sage: N = 2; A = B = [dyck_word for n in range(N+1) for dyck_word in DyckWords(n)]
            sage: tau = lambda D: D.number_of_touch_points()
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_statistics((lambda d: d.semilength(), lambda d: d.semilength()))
            sage: for solution in sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items()))):
            ....:     print(solution)
            {[]: 0, [1, 0]: 1, [1, 0, 1, 0]: 1, [1, 1, 0, 0]: 2}
            {[]: 0, [1, 0]: 1, [1, 0, 1, 0]: 2, [1, 1, 0, 0]: 1}
            sage: for subdistribution in bij.minimal_subdistributions_blocks_iterator():
            ....:     print(subdistribution)
            ([[]], [0])
            ([[1, 0]], [1])
            ([[1, 0, 1, 0], [1, 1, 0, 0]], [1, 2])

        An example with two elements of the same block in a subdistribution::

            sage: A = B = ["a", "b", "c", "d", "e"]
            sage: tau = {"a": 1, "b": 1, "c": 2, "d": 2, "e": 3}.get
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_constant_blocks([["a", "b"]])
            sage: bij.set_value_restrictions(("a", [1, 2]))
            sage: bij.constant_blocks(optimal=True)
            {{'a', 'b'}}
            sage: list(bij.minimal_subdistributions_blocks_iterator())
            [(['a', 'a', 'c', 'd', 'e'], [1, 1, 2, 2, 3])]

        An example with overlapping minimal subdistributions::

            sage: A = B = ["a", "b", "c", "d", "e"]
            sage: tau = {"a": 1, "b": 1, "c": 2, "d": 2, "e": 3}.get
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_distributions((["a", "b"], [1, 2]), (["a", "c", "d"], [1, 2, 3]))
            sage: sorted(bij.solutions_iterator(), key=lambda d: tuple(sorted(d.items())))
            [{'a': 1, 'b': 2, 'c': 2, 'd': 3, 'e': 1},
             {'a': 1, 'b': 2, 'c': 3, 'd': 2, 'e': 1},
             {'a': 2, 'b': 1, 'c': 1, 'd': 3, 'e': 2},
             {'a': 2, 'b': 1, 'c': 3, 'd': 1, 'e': 2}]
            sage: bij.constant_blocks(optimal=True)
            {{'a', 'e'}}
            sage: list(bij.minimal_subdistributions_blocks_iterator())
            [(['a', 'b'], [1, 2]), (['a', 'c', 'd'], [1, 2, 3])]

        Fedor Petrov's example from https://mathoverflow.net/q/424187::

            sage: A = B = ["a"+str(i) for i in range(1, 9)] + ["b"+str(i) for i in range(3, 9)] + ["d"]
            sage: tau = {b: 0 if i < 10 else 1 for i, b in enumerate(B)}.get
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_constant_blocks([["a"+str(i), "b"+str(i)] for i in range(1, 9) if "b"+str(i) in A])
            sage: d = [0]*8+[1]*4
            sage: bij.set_distributions((A[:8] + A[8+2:-1], d), (A[:8] + A[8:-3], d))
            sage: sorted([s[a] for a in A] for s in bij.solutions_iterator())
            [[0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1],
             [0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0],
             [0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0],
             [0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0],
             [0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0],
             [1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0],
             [1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0],
             [1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0],
             [1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0],
             [1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1],
             [1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1]]

            sage: sorted(bij.minimal_subdistributions_blocks_iterator())
            [(['a1', 'a2', 'a3', 'a4', 'a5', 'a5', 'a6', 'a6', 'a7', 'a7', 'a8', 'a8'],
              [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1]),
             (['a3', 'a4', 'd'], [0, 0, 1]),
             (['a7', 'a8', 'd'], [0, 0, 1])]

        The following solution is not found, because it happens to
        have the same support as the other::

            sage: D = set(A).difference(['b7', 'b8', 'd'])
            sage: sorted(a.replace("b", "a") for a in D)
            ['a1', 'a2', 'a3', 'a3', 'a4', 'a4', 'a5', 'a5', 'a6', 'a6', 'a7', 'a8']
            sage: set(tuple(sorted(s[a] for a in D)) for s in bij.solutions_iterator())
            {(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1)}

        But it is, by design, included here::

            sage: sorted(D) in [d for d, _ in bij.minimal_subdistributions_iterator()]
            True

        """
        # see
        # https://mathoverflow.net/questions/406751/find-a-subdistribution/406975
        # and
        # https://gitlab.com/mantepse/bijection-tools/-/issues/29
        # see https://mathoverflow.net/q/424187 for Fedor Petrov's example

        minimal_subdistribution = MixedIntegerLinearProgram(maximization=False, solver=self._solver)
        D = minimal_subdistribution.new_variable(integer=True, nonnegative=True)  # the submultiset of elements
        X = minimal_subdistribution.new_variable(binary=True)  # the support of D
        V = minimal_subdistribution.new_variable(integer=True, nonnegative=True)  # the subdistribution
        P = _disjoint_set_roots(self._P)
        minimal_subdistribution.set_objective(sum(D[p] for p in P))
        minimal_subdistribution.add_constraint(sum(D[p] for p in P) >= 1)
        for p in P:
            minimal_subdistribution.add_constraint(D[p] <= len(self._P.root_to_elements_dict()[p]))
            minimal_subdistribution.add_constraint(X[p]*len(self._P.root_to_elements_dict()[p]) >= D[p] >= X[p])

        def add_counter_example_constraint(s):
            for v in self._Z:
                minimal_subdistribution.add_constraint(sum(D[p] for p in P
                                                           if s[p] == v) == V[v])

        try:
            bmilp = self._generate_and_solve_initial_bmilp()
        except MIPSolverException:
            return
        s = self._solution_by_blocks(bmilp)
        add_counter_example_constraint(s)
        while True:
            try:
                minimal_subdistribution.solve()
            except MIPSolverException:
                return
            d = minimal_subdistribution.get_values(D)  # a dict from P to multiplicities
            new_s = self._find_counter_example2(bmilp, P, s, d)
            if new_s is None:
                yield ([p for p in P for _ in range(ZZ(d[p]))],
                       self._sorter["Z"](s[p]
                                         for p in P
                                         for _ in range(ZZ(d[p]))))

                support = [X[p] for p in P if d[p]]
                # add constraint that the support is different
                minimal_subdistribution.add_constraint(sum(support) <= len(support) - 1,
                                                       name="veto")
            else:
                s = new_s
                add_counter_example_constraint(s)

    def _find_counter_example2(self, bmilp, P, s0, d):
        r"""
        Return a solution `s` such that ``d`` is not a subdistribution of
        `s0`.

        .. TODO::

            find a better name - possibly not relevant if we
            implement the cache of solutions

        INPUT:

        - ``bmilp``, the mixed linear integer program

        - ``P``, the representatives of the blocks

        - ``s0``, a solution

        - ``d``, a subset of `A`, in the form of a dict from `A` to `\{0, 1\}`

        """
        for v in self._Z:
            v_in_d_count = sum(d[p] for p in P if s0[p] == v)
            if not v_in_d_count:
                continue

            veto_bmilp = deepcopy(bmilp)  # adding constraints to a simple copy adds them to the original instance, too
            # try to find a solution which has a different
            # subdistribution on d than s0
            v_in_d = sum(d[p] * veto_bmilp._x[p, v]
                         for p in P
                         if v in self._possible_block_values[p])

            # it is sufficient to require that v occurs less often as
            # a value among {a | d[a] == 1} than it does in
            # v_in_d_count, because, if the distributions are
            # different, one such v must exist
            veto_bmilp.milp.add_constraint(v_in_d <= v_in_d_count - 1)
            try:
                veto_bmilp.milp.solve()
                return self._solution_by_blocks(veto_bmilp)
            except MIPSolverException:
                pass
        return

    def _preprocess_intertwining_relations(self):
        r"""

        .. TODO::

            (medium) untangle side effect and return value if possible

        Make `self._P` be the finest set partition coarser than `self._P`
        such that composing elements preserves blocks.

        Suppose that `p_1`, `p_2` are blocks of `P`, and `a_1, a'_1
        \in p_1` and `a_2, a'_2\in p_2`.  Then,

        .. MATH:

            s(\pi(a_1, a_2))
            = \rho(s(a_1), s(a_2))
            = \rho(s(a'_1), s(a'_2))
            = s(\pi(a'_1, a'_2)).

        Therefore, `\pi(a_1, a_2)` and `\pi(a'_1, a'_2)` are in the
        same block.

        In other words, `s(\pi(a_1,\dots,a_k))` only depends on the
        blocks of `a_1,\dots,a_k`.

        .. TODO::

            create one test with one and one test with two
            intertwining_relations

        """
        images = {}  # A^k -> A, a_1,...,a_k to pi(a_1,...,a_k), for all pi
        origins_by_elements = []  # (pi/rho, pi(a_1,...,a_k), a_1,...,a_k)
        for composition_index, pi_rho in enumerate(self._pi_rho):
            for a_tuple in itertools.product(*([self._A]*pi_rho.numargs)):
                if pi_rho.domain is not None and not pi_rho.domain(*a_tuple):
                    continue
                a = pi_rho.pi(*a_tuple)
                if a in self._A:
                    if a in images:
                        # this happens if there are several pi's of the same arity
                        images[a_tuple].add(a)  # TODO: (low) wouldn't self._P.find(a) be more efficient here?
                    else:
                        images[a_tuple] = set((a,))  # TODO: (low) wouldn't self._P.find(a) be more efficient here?
                    origins_by_elements.append((composition_index, a, a_tuple))

        # merge blocks
        something_changed = True
        while something_changed:
            something_changed = False
            # collect (preimage, image) pairs by (representatives) of
            # the blocks of the elements of the preimage
            updated_images = {}  # (p_1,...,p_k) to {a_1,....}
            for a_tuple, image_set in images.items():
                representatives = tuple(self._P.find(a) for a in a_tuple)
                if representatives in updated_images:
                    updated_images[representatives].update(image_set)
                else:
                    updated_images[representatives] = image_set

            # merge blocks
            for a_tuple, image_set in updated_images.items():
                image = image_set.pop()
                while image_set:
                    self._P.union(image, image_set.pop())
                    something_changed = True
                # we keep a representative
                image_set.add(image)

            images = updated_images

        origins = set()
        for composition_index, image, preimage in origins_by_elements:
            origins.add((composition_index,
                         self._P.find(image),
                         tuple(self._P.find(a) for a in preimage)))
        return origins

    def solutions_iterator(self):
        r"""
        An iterator over all solutions of the problem.

        OUTPUT: An iterator over all possible mappings `s: A\to Z`

        ALGORITHM:

        We solve an integer linear program with a binary variable
        `x_{p, z}` for each partition block `p\in P` and each
        statistic value `z\in Z`:

        - `x_{p, z} = 1` if and only if `s(a) = z` for all `a\in p`.

        Then we add the constraint `\sum_{x\in V} x<|V|`, where `V`
        is the set containing all `x` with `x = 1`, that is, those
        indicator variables representing the current solution.
        Therefore, a solution of this new program must be different
        from all those previously obtained.

        INTEGER LINEAR PROGRAM:

        * Let `m_w(p)`, for a block `p` of `P`, be the multiplicity
          of the value `w` in `W` under `\alpha`, that is, the number
          of elements `a \in p` with `\alpha(a)=w`.

        * Let `n_w(z)` be the number of elements `b \in B` with
          `\beta(b)=w` and `\tau(b)=z` for `w \in W`, `z \in Z`.

        * Let `k` be the arity of a pair `(\pi, \rho)` in an
          intertwining relation.

        and the following constraints:

        * because every block is assigned precisely one value, for
          all `p\in P`,

        .. MATH::

            \sum_z x_{p, z} = 1.

        * because the statistics `s` and `\tau` and also `\alpha` and
          `\beta` are equidistributed, for all `w\in W` and `z\in Z`,

        .. MATH::

            \sum_p m_w(p) x_{p, z} = n_w(z).

        * for each intertwining relation `s(\pi(a_1,\dots, a_k)) =
          \rho(s(a_1),\dots, s(a_r))`, and for all `k`-combinations
          of blocks `p_i\in P` such that there exist `(a_1,\dots,
          a_k)\in p_1\times\dots\times p_k` with `\pi(a_1,\dots,
          a_k)\in W` and `z = \rho(z_1,\dots, z_k)`,

        .. MATH::

            x_{p, z} \geq 1-k + \sum_{i=1}^k x_{p_i, z_i}.

        * for each distribution restriction, i.e. a set of elements
          `e` and a distribution of values given by integers `d_z`
          representing the multiplicity of each `z \in Z`, and `r_p =
          |p \cap e|` indicating the relative size of block `p` in
          the set of elements of the distribution,

        .. MATH::

            \sum_p r_p x_{p, z} = d_z.

        EXAMPLES::

            sage: A = B = list('abc')
            sage: bij = Bijectionist(A, B, lambda x: B.index(x) % 2, solver="GLPK")
            sage: next(bij.solutions_iterator())
            {'a': 0, 'b': 1, 'c': 0}

            sage: list(bij.solutions_iterator())
            [{'a': 0, 'b': 1, 'c': 0},
             {'a': 1, 'b': 0, 'c': 0},
             {'a': 0, 'b': 0, 'c': 1}]

            sage: N = 4
            sage: A = B = [permutation for n in range(N) for permutation in Permutations(n)]

        Let `\tau` be the number of non-left-to-right-maxima of a
        permutation::

            sage: def tau(pi):
            ....:    pi = list(pi)
            ....:    i = count = 0
            ....:    for j in range(len(pi)):
            ....:        if pi[j] > i:
            ....:            i = pi[j]
            ....:        else:
            ....:            count += 1
            ....:    return count

        We look for a statistic which is constant on conjugacy classes::

            sage: P = [list(a) for n in range(N) for a in Permutations(n).conjugacy_classes()]

            sage: bij = Bijectionist(A, B, tau, solver="GLPK")
            sage: bij.set_statistics((len, len))
            sage: bij.set_constant_blocks(P)
            sage: for solution in bij.solutions_iterator():
            ....:     print(solution)
            {[]: 0, [1]: 0, [1, 2]: 1, [2, 1]: 0, [1, 2, 3]: 0, [1, 3, 2]: 1, [2, 1, 3]: 1, [3, 2, 1]: 1, [2, 3, 1]: 2, [3, 1, 2]: 2}
            {[]: 0, [1]: 0, [1, 2]: 0, [2, 1]: 1, [1, 2, 3]: 0, [1, 3, 2]: 1, [2, 1, 3]: 1, [3, 2, 1]: 1, [2, 3, 1]: 2, [3, 1, 2]: 2}

        Setting the verbosity prints the MILP which is solved::

            sage: set_verbose(2)
            sage: _ = list(bij.solutions_iterator())
            Constraints are:
                block []: 1 <= x_0 <= 1
                block [1]: 1 <= x_1 <= 1
                block [1, 2]: 1 <= x_2 + x_3 <= 1
                block [2, 1]: 1 <= x_4 + x_5 <= 1
                block [1, 2, 3]: 1 <= x_6 + x_7 + x_8 <= 1
                block [1, 3, 2]: 1 <= x_9 + x_10 + x_11 <= 1
                block [2, 3, 1]: 1 <= x_12 + x_13 + x_14 <= 1
                statistics: 1 <= x_0 <= 1
                statistics: 1 <= x_1 <= 1
                statistics: 1 <= x_2 + x_4 <= 1
                statistics: 1 <= x_3 + x_5 <= 1
                statistics: 1 <= x_6 + 3 x_9 + 2 x_12 <= 1
                statistics: 3 <= x_7 + 3 x_10 + 2 x_13 <= 3
                statistics: 2 <= x_8 + 3 x_11 + 2 x_14 <= 2
            Variables are:
                x_0: s([]) = 0
                x_1: s([1]) = 0
                x_2: s([1, 2]) = 0
                x_3: s([1, 2]) = 1
                x_4: s([2, 1]) = 0
                x_5: s([2, 1]) = 1
                x_6: s([1, 2, 3]) = 0
                x_7: s([1, 2, 3]) = 1
                x_8: s([1, 2, 3]) = 2
                x_9: s([1, 3, 2]) = s([2, 1, 3]) = s([3, 2, 1]) = 0
                x_10: s([1, 3, 2]) = s([2, 1, 3]) = s([3, 2, 1]) = 1
                x_11: s([1, 3, 2]) = s([2, 1, 3]) = s([3, 2, 1]) = 2
                x_12: s([2, 3, 1]) = s([3, 1, 2]) = 0
                x_13: s([2, 3, 1]) = s([3, 1, 2]) = 1
                x_14: s([2, 3, 1]) = s([3, 1, 2]) = 2
            after vetoing
            Constraints are:
                block []: 1 <= x_0 <= 1
                block [1]: 1 <= x_1 <= 1
                block [1, 2]: 1 <= x_2 + x_3 <= 1
                block [2, 1]: 1 <= x_4 + x_5 <= 1
                block [1, 2, 3]: 1 <= x_6 + x_7 + x_8 <= 1
                block [1, 3, 2]: 1 <= x_9 + x_10 + x_11 <= 1
                block [2, 3, 1]: 1 <= x_12 + x_13 + x_14 <= 1
                statistics: 1 <= x_0 <= 1
                statistics: 1 <= x_1 <= 1
                statistics: 1 <= x_2 + x_4 <= 1
                statistics: 1 <= x_3 + x_5 <= 1
                statistics: 1 <= x_6 + 3 x_9 + 2 x_12 <= 1
                statistics: 3 <= x_7 + 3 x_10 + 2 x_13 <= 3
                statistics: 2 <= x_8 + 3 x_11 + 2 x_14 <= 2
                veto: x_0 + x_1 + x_3 + x_4 + x_6 + x_10 + x_14 <= 6
            after vetoing
            Constraints are:
                block []: 1 <= x_0 <= 1
                block [1]: 1 <= x_1 <= 1
                block [1, 2]: 1 <= x_2 + x_3 <= 1
                block [2, 1]: 1 <= x_4 + x_5 <= 1
                block [1, 2, 3]: 1 <= x_6 + x_7 + x_8 <= 1
                block [1, 3, 2]: 1 <= x_9 + x_10 + x_11 <= 1
                block [2, 3, 1]: 1 <= x_12 + x_13 + x_14 <= 1
                statistics: 1 <= x_0 <= 1
                statistics: 1 <= x_1 <= 1
                statistics: 1 <= x_2 + x_4 <= 1
                statistics: 1 <= x_3 + x_5 <= 1
                statistics: 1 <= x_6 + 3 x_9 + 2 x_12 <= 1
                statistics: 3 <= x_7 + 3 x_10 + 2 x_13 <= 3
                statistics: 2 <= x_8 + 3 x_11 + 2 x_14 <= 2
                veto: x_0 + x_1 + x_3 + x_4 + x_6 + x_10 + x_14 <= 6
                veto: x_0 + x_1 + x_2 + x_5 + x_6 + x_10 + x_14 <= 6

            sage: set_verbose(0)

        TESTS:

        An unfeasible problem::

            sage: A = ["a", "b", "c", "d"]; B = [1, 2, 3, 4]
            sage: bij = Bijectionist(A, B)
            sage: bij.set_value_restrictions(("a", [1, 2]), ("b", [1, 2]), ("c", [1, 3]), ("d", [2, 3]))
            sage: list(bij.solutions_iterator())
            []

        """
        try:
            bmilp = self._generate_and_solve_initial_bmilp()
        except MIPSolverException:
            return
        while True:
            yield self._solution(bmilp)
            bmilp.veto_current_solution()
            if get_verbose() >= 2:
                print("after vetoing")
                self._show_bmilp(bmilp, variables=False)
            try:
                bmilp.milp.solve()
            except MIPSolverException:
                return

    def _solution(self, bmilp):
        """
        Return the bmilp solution as a dictionary from `A` to
        `Z`.

        """
        map = {}  # A -> Z, a +-> s(a)
        for p, block in self._P.root_to_elements_dict().items():
            for z in self._possible_block_values[p]:
                if bmilp.milp.get_values(bmilp._x[p, z]) == 1:
                    for a in block:
                        map[a] = z
                    break
        return map

    def _solution_by_blocks(self, bmilp):
        """
        Return the bmilp solution as a dictionary from block
        representatives of `P` to `Z`.

        """
        map = {}  # P -> Z, a +-> s(a)
        for p in _disjoint_set_roots(self._P):
            for z in self._possible_block_values[p]:
                if bmilp.milp.get_values(bmilp._x[p, z]) == 1:
                    map[p] = z
                    break
        return map

    def _show_bmilp(self, bmilp, variables=True):
        """
        Print the constraints and variables of the current MILP
        together with some explanations.

        """
        print("Constraints are:")
        b = bmilp.milp.get_backend()
        varid_name = {}
        for i in range(b.ncols()):
            s = b.col_name(i)
            default_name = str(bmilp.milp.linear_functions_parent()({i: 1}))
            if s and s != default_name:
                varid_name[i] = s
            else:
                varid_name[i] = default_name
        for i, (lb, (indices, values), ub) in enumerate(bmilp.milp.constraints()):
            if b.row_name(i):
                print("    "+b.row_name(i)+":", end=" ")
            if lb is not None:
                print(str(ZZ(lb))+" <=", end=" ")
            first = True
            for j, c in sorted(zip(indices, values)):
                c = ZZ(c)
                if c == 0:
                    continue
                print((("+ " if (not first and c > 0) else "") +
                       ("" if c == 1 else
                        ("- " if c == -1 else
                         (str(c) + " " if first and c < 0 else
                          ("- " + str(abs(c)) + " " if c < 0 else str(c) + " "))))
                       + varid_name[j]), end=" ")
                first = False
            # Upper bound
            print("<= "+str(ZZ(ub)) if ub is not None else "")

        if variables:
            print("Variables are:")
            for (p, z), v in bmilp._x.items():
                print(f"    {v}: " + "".join([f"s({a}) = "
                                              for a in self._P.root_to_elements_dict()[p]]) + f"{z}")

    def _generate_and_solve_initial_bmilp(self):
        r"""
        Generate a _BijectionistMILP, add all relevant constraints and call MILP.solve().
        """
        preimage_blocks = self._preprocess_intertwining_relations()
        self._compute_possible_block_values()

        bmilp = _BijectionistMILP(self)
        n = bmilp.milp.number_of_variables()
        bmilp.add_alpha_beta_constraints()
        bmilp.add_distribution_constraints()
        bmilp.add_interwining_relation_constraints(preimage_blocks)
        if get_verbose() >= 2:
            self._show_bmilp(bmilp)
        assert n == bmilp.milp.number_of_variables(), "The number of variables increased."
        bmilp.milp.solve()
        return bmilp


class _BijectionistMILP(SageObject):
    r"""
    Wrapper class for the MixedIntegerLinearProgram (MILP).  This class is used to manage the MILP,
    add constraints, solve the problem and check for uniqueness of solution values.
    """
    def __init__(self, bijectionist: Bijectionist):
        # TODO: it would be cleaner not to pass the full bijectionist
        # instance, but only those attributes we actually use:
        # _possible_block_values
        # _elements_distributions
        # _W, _Z, _A, _B, _P, _alpha, _beta, _tau, _pi_rho
        self.milp = MixedIntegerLinearProgram(solver=bijectionist._solver)
        self.milp.set_objective(None)
        self._x = self.milp.new_variable(binary=True)  # indexed by P x Z

        self._bijectionist = bijectionist

        for p in _disjoint_set_roots(bijectionist._P):
            name = f"block {p}"
            self.milp.add_constraint(sum(self._x[p, z]
                                         for z in bijectionist._possible_block_values[p]) == 1,
                                     name=name[:50])

    def add_alpha_beta_constraints(self):
        r"""
        Add constraints enforcing that `(alpha, s)` is equidistributed
        with `(beta, tau)` and `S` is the intertwining bijection.

        We do this by adding

        .. MATH::

            \sum_{a\in A, z\in Z} x_{p(a), z} s^z t^{\alpha(a)}
            = \sum_{b\in B} s^{\tau(b)} t(\beta(b))

        as a matrix equation.

        """
        W = self._bijectionist._W
        Z = self._bijectionist._Z
        AZ_matrix = [[ZZ(0)]*len(W) for _ in range(len(Z))]
        B_matrix = [[ZZ(0)]*len(W) for _ in range(len(Z))]

        W_dict = {w: i for i, w in enumerate(W)}
        Z_dict = {z: i for i, z in enumerate(Z)}

        for a in self._bijectionist._A:
            p = self._bijectionist._P.find(a)
            for z in self._bijectionist._possible_block_values[p]:
                w_index = W_dict[self._bijectionist._alpha(a)]
                z_index = Z_dict[z]
                AZ_matrix[z_index][w_index] += self._x[p, z]

        for b in self._bijectionist._B:
            w_index = W_dict[self._bijectionist._beta(b)]
            z_index = Z_dict[self._bijectionist._tau[b]]
            B_matrix[z_index][w_index] += 1

        # TODO: (low) I am not sure that this is the best way to
        # filter out empty conditions
        for w in range(len(W)):
            for z in range(len(Z)):
                c = AZ_matrix[z][w] - B_matrix[z][w]
                if c.is_zero():
                    continue
                if c in ZZ:
                    raise MIPSolverException
                self.milp.add_constraint(c == 0, name="statistics")

    def add_distribution_constraints(self):
        r"""
        Add constraints so the distributions given by
        :meth:`~Bijectionist.set_distributions` are fulfilled.

        To accomplish this we add

        .. MATH::

            \sum_{a\in elements} x_{p(a), z}t^z = \sum_{z\in values} t^z,

        where `p(a)` is the block containing `a`, for each given
        distribution as a vector equation.

        """
        Z = self._bijectionist._Z
        Z_dict = {z: i for i, z in enumerate(Z)}
        for elements, values in self._bijectionist._elements_distributions:
            elements_sum = [ZZ(0)]*len(Z_dict)
            values_sum = [ZZ(0)]*len(Z_dict)
            for a in elements:
                p = self._bijectionist._P.find(a)
                for z in self._bijectionist._possible_block_values[p]:
                    elements_sum[Z_dict[z]] += self._x[p, z]
            for z in values:
                values_sum[Z_dict[z]] += 1

            # TODO: (low) I am not sure that this is the best way to
            # filter out empty conditions
            for element, value in zip(elements_sum, values_sum):
                c = element - value
                if c.is_zero():
                    continue
                if c in ZZ:
                    raise MIPSolverException
                self.milp.add_constraint(c == 0, name=f"d: {element} == {value}")

    def add_interwining_relation_constraints(self, origins):
        r"""
        Add constraints corresponding to the given intertwining
        relations.

        This adds the constraints imposed by
        :meth:`~Bijectionist.set_intertwining_relations`.

        .. MATH::

            s(\pi(a_1,\dots, a_k)) = \rho(s(a_1),\dots, s(a_k))`

        for each pair `(\pi, \rho)`.  The relation implies
        immediately that `s(\pi(a_1,\dots, a_k))` only depends on the
        blocks of `a_1,\dots, a_k`.

        The MILP formulation is as follows.  Let `a_1,\dots,a_k \in
        A` and let `a = \pi(a_1,\dots,a_k)`.  Let `z_1,\dots,z_k \in
        Z` and let `z = \rho(z_1,\dots,z_k)`.  Suppose that `a_i\in
        p_i` for all `i` and that `a\in p`.

        We then want to model the implication

        .. MATH::

           x_{p_1, z_1} = 1,\dots, x_{p_k, z_k} = 1 \Rightarrow x_{p, z} = 1.

        We achieve this by requiring

        .. MATH::

            x_{p, z}\geq 1 - k + \sum_{i=1}^k x_{p_i, z_i}.

        Not that `z` must be a possible value of `p` and each `z_i`
        must be a possible value of `p_i`.

        INPUT:

        - origins, a list of triples `((\pi/\rho, p,
          (p_1,\dots,p_k))`, where `p` is the block of
          `\rho(s(a_1),\dots, s(a_k))`, for any `a_i\in p_i`.

        """
        for composition_index, image_block, preimage_blocks in origins:
            pi_rho = self._bijectionist._pi_rho[composition_index]
            # iterate over all possible value combinations of the origin blocks
            for z_tuple in itertools.product(*[self._bijectionist._possible_block_values[p]
                                               for p in preimage_blocks]):
                rhs = 1 - pi_rho.numargs + sum(self._x[p_i, z_i]
                                               for p_i, z_i in zip(preimage_blocks, z_tuple))
                z = pi_rho.rho(*z_tuple)
                if z in self._bijectionist._possible_block_values[image_block]:
                    c = self._x[image_block, z] - rhs
                    if c.is_zero():
                        continue
                    self.milp.add_constraint(c >= 0,
                                             name=f"pi/rho({composition_index})")
                else:
                    self.milp.add_constraint(rhs <= 0,
                                             name=f"pi/rho({composition_index})")

    def veto_current_solution(self):
        r"""
        Add a constraint vetoing the current solution.

        This adds a constraint such that the next call to
        :meth:`MixedIntegerLinearProgram.solve()` must return a
        solution different from the current one.

        We require that the MILP currently has a solution.

        .. WARNING::

            The underlying MILP will be modified!

        ALGORITHM:

        We add the constraint `\sum_{x\in V} x < |V|`` where `V` is
        the set of variables `x_{p, z}` with value 1, that is, the
        set of variables indicating the current solution.

        """
        # get all variables with value 1
        active_vars = [self._x[p, z]
                       for p in _disjoint_set_roots(self._bijectionist._P)
                       for z in self._bijectionist._possible_block_values[p]
                       if self.milp.get_values(self._x[p, z])]

        # add constraint that not all of these can be 1, thus vetoing
        # the current solution
        self.milp.add_constraint(sum(active_vars) <= len(active_vars) - 1,
                                 name="veto")


def _invert_dict(d):
    """
    Return the dictionary whose keys are the values of the input and
    whose values are the lists of preimages.
    """
    preimages = {}
    for k, v in d.items():
        preimages[v] = preimages.get(v, []) + [k]
    return preimages


def _disjoint_set_roots(d):
    """
    Return the representatives of the blocks of the disjoint set.
    """
    return d.root_to_elements_dict().keys()


"""
TESTS::

    sage: As = Bs = [[],
    ....:            [(1,i,j) for i in [-1,0,1] for j in [-1,1]],
    ....:            [(2,i,j) for i in [-1,0,1] for j in [-1,1]],
    ....:            [(3,i,j) for i in [-2,-1,0,1,2] for j in [-1,1]]]

    # adding [(2,-2,-1), (2,2,-1), (2,-2,1), (2,2,1)] makes it take (seemingly) forever

    sage: c1 = lambda a, b: (a[0]+b[0], a[1]*b[1], a[2]*b[2])
    sage: c2 = lambda a: (a[0], -a[1], a[2])

    sage: bij = Bijectionist(sum(As, []), sum(Bs, []))
    sage: bij.set_statistics((lambda x: x[0], lambda x: x[0]))
    sage: bij.set_intertwining_relations((2, c1, c1), (1, c2, c2))
    sage: l = list(bij.solutions_iterator()); len(l)
    64

A brute force check would be difficult::

    sage: prod([factorial(len(A)) for A in As])
    1881169920000

Let us try a smaller example::

    sage: As = Bs = [[],
    ....:            [(1,i,j) for i in [-1,0,1] for j in [-1,1]],
    ....:            [(2,i,j) for i in [-1,1] for j in [-1,1]],
    ....:            [(3,i,j) for i in [-1,1] for j in [-1,1]]]

    sage: bij = Bijectionist(sum(As, []), sum(Bs, []))
    sage: bij.set_statistics((lambda x: x[0], lambda x: x[0]))
    sage: bij.set_intertwining_relations((2, c1, c1), (1, c2, c2))
    sage: l1 = list(bij.solutions_iterator()); len(l1)
    16
    sage: prod([factorial(len(A)) for A in As])
    414720

    sage: pis = cartesian_product([Permutations(len(A)) for A in As])
    sage: it = ({a: Bs[n][pi[n][i]-1] for n, A in enumerate(As) for i, a in enumerate(A)} for pi in pis)
    sage: A = sum(As, [])
    sage: respects_c1 = lambda s: all(c1(a1, a2) not in A or s[c1(a1, a2)] == c1(s[a1], s[a2]) for a1 in A for a2 in A)
    sage: respects_c2 = lambda s: all(c2(a1) not in A or s[c2(a1)] == c2(s[a1]) for a1 in A)
    sage: l2 = [s for s in it if respects_c1(s) and respects_c2(s)]
    sage: sorted(l1, key=lambda s: tuple(s.items())) == l2
    True

"""


"""
Our benchmark example::

    sage: from sage.combinat.cyclic_sieving_phenomenon import orbit_decomposition
    sage: alpha1 = lambda p: len(p.weak_excedences())
    sage: alpha2 = lambda p: len(p.fixed_points())
    sage: beta1 = lambda p: len(p.descents(final_descent=True)) if p else 0
    sage: beta2 = lambda p: len([e for (e, f) in zip(p, p[1:]+[0]) if e == f+1])
    sage: gamma = Permutation.longest_increasing_subsequence_length
    sage: def rotate_permutation(p):
    ....:    cycle = Permutation(tuple(range(1, len(p)+1)))
    ....:    return Permutation([cycle.inverse()(p(cycle(i))) for i in range(1, len(p)+1)])

    sage: N=5
    sage: As = [list(Permutations(n)) for n in range(N+1)]
    sage: A = B = sum(As, [])
    sage: bij = Bijectionist(A, B, gamma)
    sage: bij.set_statistics((len, len), (alpha1, beta1), (alpha2, beta2))
    sage: bij.set_constant_blocks(sum([orbit_decomposition(A, rotate_permutation) for A in As], []))

    sage: P = bij.constant_blocks(optimal=True)
    sage: P = [sorted(p, key=lambda p: (len(p), p)) for p in P]
    sage: P = sorted(P, key=lambda p: (len(next(iter(p))), len(p)))
    sage: for p in P:
    ....:     print(p)
    [[1], [1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 2, 3, 4, 5]]
    [[2, 1], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 2, 1], ...
    [[3, 1, 2], [1, 4, 2, 3], [2, 4, 1, 3], [3, 1, 2, 4], [3, 1, 4, 2], ...
    [[4, 1, 2, 3], [1, 5, 2, 3, 4], [4, 1, 2, 3, 5], [4, 5, 1, 2, 3], [5, 1, 2, 4, 3], ...
    [[1, 3, 2, 5, 4], [2, 1, 3, 5, 4], [2, 1, 4, 3, 5], [5, 2, 4, 3, 1], [5, 3, 2, 4, 1]]
    [[1, 3, 5, 2, 4], [2, 4, 1, 3, 5], [3, 5, 2, 4, 1], [4, 1, 3, 5, 2], [5, 2, 4, 1, 3]]
    ...

    sage: for d in sorted(bij.minimal_subdistributions_blocks_iterator(), key=lambda d: (len(d[0]), d[0])):
    ....:     print(d)
    ([[]], [0])
    ([[1]], [1])
    ([[2, 1]], [2])
    ([[3, 1, 2]], [3])
    ([[4, 1, 2, 3]], [4])
    ([[5, 1, 2, 3, 4]], [5])
    ([[2, 1, 4, 5, 3], [2, 3, 5, 1, 4], [2, 4, 1, 5, 3], [2, 4, 5, 1, 3]], [2, 3, 3, 3])
    ([[2, 1, 5, 3, 4], [2, 5, 1, 3, 4], [3, 1, 5, 2, 4], [3, 5, 1, 2, 4]], [3, 3, 4, 4])
    ([[1, 3, 2, 5, 4], [1, 3, 5, 2, 4], [1, 4, 2, 5, 3], [1, 4, 5, 2, 3], [1, 4, 5, 3, 2], [1, 5, 4, 2, 3], [1, 5, 4, 3, 2]], [2, 2, 3, 3, 3, 3, 3])

    sage: l = list(bij.solutions_iterator()); len(l)  # long time
    504

    sage: for a, d in bij.minimal_subdistributions_iterator():  # long time
    ....:     print(sorted(a), sorted(d))
"""
