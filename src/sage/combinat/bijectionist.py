# sage.doctest: needs sage.numerical.mip
r"""
A bijectionist's toolkit

AUTHORS:

- Alexander Grosz, Tobias Kietreiber, Stephan Pfannerer and Martin
  Rubey (2020-2022): Initial version

Quick reference
===============

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~Bijectionist.set_statistics` | Declare statistics that are preserved by the bijection.
    :meth:`~Bijectionist.set_value_restrictions` | Restrict the values of the statistic on an element.
    :meth:`~Bijectionist.set_constant_blocks` | Declare that the statistic is constant on some sets.
    :meth:`~Bijectionist.set_distributions` | Restrict the distribution of values of the statistic on some elements.
    :meth:`~Bijectionist.set_intertwining_relations` | Declare that the statistic intertwines with other maps.
    :meth:`~Bijectionist.set_quadratic_relation` | Declare that the statistic satisfies a certain relation.
    :meth:`~Bijectionist.set_homomesic` | Declare that the statistic is homomesic with respect to a given set partition.


    :meth:`~Bijectionist.statistics_table` | Print a table collecting information on the given statistics.
    :meth:`~Bijectionist.statistics_fibers` | Collect elements with the same statistics.

    :meth:`~Bijectionist.constant_blocks` | Return the blocks on which the statistic is constant.
    :meth:`~Bijectionist.solutions_iterator` | Iterate over all possible solutions.
    :meth:`~Bijectionist.possible_values` | Return all possible values for a given element.
    :meth:`~Bijectionist.minimal_subdistributions_iterator` | Iterate over the minimal subdistributions.
    :meth:`~Bijectionist.minimal_subdistributions_blocks_iterator` | Iterate over the minimal subdistributions.

A guided tour
=============

Consider the following combinatorial statistics on a permutation:

    * `wex`, the number of weak excedences,
    * `fix`, the number of fixed points,
    * `des`, the number of descents (after appending `0`),
    * `adj`, the number of adjacencies (after appending `0`), and
    * `llis`, the length of a longest increasing subsequence.

Moreover, let `rot` be action of rotation on a permutation, i.e., the
conjugation with the long cycle.

It is known that there must exist a statistic `s` on permutations,
which is equidistributed with `llis` but additionally invariant under
`rot`.  However, at least very small cases do not contradict the
possibility that one can even find a statistic `s`, invariant under
`rot` and such that `(s, wex, fix) \sim (llis, des, adj)`.  Let us
check this for permutations of size at most `3`::

    sage: N = 3
    sage: A = B = [pi for n in range(N+1) for pi in Permutations(n)]
    sage: def alpha1(p): return len(p.weak_excedences())
    sage: def alpha2(p): return len(p.fixed_points())
    sage: def beta1(p): return len(p.descents(final_descent=True)) if p else 0
    sage: def beta2(p): return len([e for (e, f) in zip(p, p[1:]+[0]) if e == f+1])
    sage: tau = Permutation.longest_increasing_subsequence_length
    sage: def rotate_permutation(p):
    ....:     cycle = Permutation(tuple(range(1, len(p)+1)))
    ....:     return Permutation([cycle.inverse()(p(cycle(i))) for i in range(1, len(p)+1)])
    sage: bij = Bijectionist(A, B, tau)
    sage: bij.set_statistics((len, len), (alpha1, beta1), (alpha2, beta2))
    sage: a, b = bij.statistics_table()
    sage: table(a, header_row=True, frame=True)
    ┌───────────┬────────┬────────┬────────┐
    │ a         │ α_1(a) │ α_2(a) │ α_3(a) │
    ╞═══════════╪════════╪════════╪════════╡
    │ []        │ 0      │ 0      │ 0      │
    ├───────────┼────────┼────────┼────────┤
    │ [1]       │ 1      │ 1      │ 1      │
    ├───────────┼────────┼────────┼────────┤
    │ [1, 2]    │ 2      │ 2      │ 2      │
    ├───────────┼────────┼────────┼────────┤
    │ [2, 1]    │ 2      │ 1      │ 0      │
    ├───────────┼────────┼────────┼────────┤
    │ [1, 2, 3] │ 3      │ 3      │ 3      │
    ├───────────┼────────┼────────┼────────┤
    │ [1, 3, 2] │ 3      │ 2      │ 1      │
    ├───────────┼────────┼────────┼────────┤
    │ [2, 1, 3] │ 3      │ 2      │ 1      │
    ├───────────┼────────┼────────┼────────┤
    │ [2, 3, 1] │ 3      │ 2      │ 0      │
    ├───────────┼────────┼────────┼────────┤
    │ [3, 1, 2] │ 3      │ 1      │ 0      │
    ├───────────┼────────┼────────┼────────┤
    │ [3, 2, 1] │ 3      │ 2      │ 1      │
    └───────────┴────────┴────────┴────────┘

    sage: table(b, header_row=True, frame=True)
    ┌───────────┬───┬────────┬────────┬────────┐
    │ b         │ τ │ β_1(b) │ β_2(b) │ β_3(b) │
    ╞═══════════╪═══╪════════╪════════╪════════╡
    │ []        │ 0 │ 0      │ 0      │ 0      │
    ├───────────┼───┼────────┼────────┼────────┤
    │ [1]       │ 1 │ 1      │ 1      │ 1      │
    ├───────────┼───┼────────┼────────┼────────┤
    │ [1, 2]    │ 2 │ 2      │ 1      │ 0      │
    ├───────────┼───┼────────┼────────┼────────┤
    │ [2, 1]    │ 1 │ 2      │ 2      │ 2      │
    ├───────────┼───┼────────┼────────┼────────┤
    │ [1, 2, 3] │ 3 │ 3      │ 1      │ 0      │
    ├───────────┼───┼────────┼────────┼────────┤
    │ [1, 3, 2] │ 2 │ 3      │ 2      │ 1      │
    ├───────────┼───┼────────┼────────┼────────┤
    │ [2, 1, 3] │ 2 │ 3      │ 2      │ 1      │
    ├───────────┼───┼────────┼────────┼────────┤
    │ [2, 3, 1] │ 2 │ 3      │ 2      │ 1      │
    ├───────────┼───┼────────┼────────┼────────┤
    │ [3, 1, 2] │ 2 │ 3      │ 2      │ 0      │
    ├───────────┼───┼────────┼────────┼────────┤
    │ [3, 2, 1] │ 1 │ 3      │ 3      │ 3      │
    └───────────┴───┴────────┴────────┴────────┘

    sage: from sage.combinat.cyclic_sieving_phenomenon import orbit_decomposition
    sage: bij.set_constant_blocks(orbit_decomposition(A, rotate_permutation))
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

On the other hand, we can check that there is no rotation invariant
statistic on non-crossing set partitions which is equidistributed
with the Strahler number on ordered trees::

    sage: N = 8
    sage: A = [SetPartition(d.to_noncrossing_partition()) for n in range(N) for d in DyckWords(n)]
    sage: B = [t for n in range(1, N+1) for t in OrderedTrees(n)]
    sage: def theta(m): return SetPartition([[i % m.size() + 1 for i in b] for b in m])

Code for the Strahler number can be obtained from FindStat.  The
following code is equivalent to ``tau = findstat(397)``::

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
    sage: bij.set_constant_blocks(orbit_decomposition(A, theta))
    sage: list(bij.solutions_iterator())
    []

Next we demonstrate how to search for a bijection, instead An example
identifying `s` and `S`::

    sage: N = 4
    sage: A = [dyck_word for n in range(1, N) for dyck_word in DyckWords(n)]
    sage: B = [binary_tree for n in range(1, N) for binary_tree in BinaryTrees(n)]
    sage: concat_path = lambda D1, D2: DyckWord(list(D1) + list(D2))
    sage: concat_tree = lambda B1, B2: concat_path(B1.to_dyck_word(),
    ....:                                          B2.to_dyck_word()).to_binary_tree()
    sage: bij = Bijectionist(A, B)
    sage: bij.set_intertwining_relations((2, concat_path, concat_tree))
    sage: bij.set_statistics((lambda d: d.semilength(), lambda t: t.node_number()))
    sage: for D in sorted(bij.minimal_subdistributions_iterator(), key=lambda x: (len(x[0][0]), x)):
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

The output is in a form suitable for FindStat::

    sage: findmap(list(bij.minimal_subdistributions_iterator()))            # optional -- internet
    0: Mp00034 (quality [100])
    1: Mp00061oMp00023 (quality [100])
    2: Mp00018oMp00140 (quality [100])

TESTS::

    sage: N = 4; A = B = [permutation for n in range(N) for permutation in Permutations(n)]
    sage: def theta(pi): return Permutation([x+1 if x != len(pi) else 1 for x in pi[-1:]+pi[:-1]])
    sage: def tau(pi):
    ....:    n = len(pi)
    ....:    return sum([1 for i in range(1, n+1) for j in range(1, n+1)
    ....:                if i<j <= pi(i)<pi(j) or pi(i)<pi(j)<i<j])
    sage: bij = Bijectionist(A, B, tau)
    sage: bij.set_statistics((len, len))
    sage: from sage.combinat.cyclic_sieving_phenomenon import orbit_decomposition
    sage: bij.set_constant_blocks(orbit_decomposition(A, theta))
    sage: for solution in sorted(bij.solutions_iterator(), key=lambda d: sorted(d.items())):
    ....:     print(solution)
    {[]: 0, [1]: 0, [1, 2]: 0, [2, 1]: 0, [1, 2, 3]: 0, [1, 3, 2]: 0, [2, 1, 3]: 0, [3, 2, 1]: 0, [2, 3, 1]: 0, [3, 1, 2]: 1}
    {[]: 0, [1]: 0, [1, 2]: 0, [2, 1]: 0, [1, 2, 3]: 0, [1, 3, 2]: 0, [2, 1, 3]: 0, [3, 2, 1]: 0, [2, 3, 1]: 1, [3, 1, 2]: 0}
    {[]: 0, [1]: 0, [1, 2]: 0, [2, 1]: 0, [1, 2, 3]: 1, [1, 3, 2]: 0, [2, 1, 3]: 0, [3, 2, 1]: 0, [2, 3, 1]: 0, [3, 1, 2]: 0}

A test including intertwining relations::

    sage: N = 2; A = B = [dyck_word for n in range(N+1) for dyck_word in DyckWords(n)]
    sage: def alpha(D): return (D.area(), D.bounce())
    sage: def beta(D): return (D.bounce(), D.area())
    sage: def tau(D): return D.number_of_touch_points()

The following looks correct::

    sage: bij = Bijectionist(A, B, tau)
    sage: bij.set_statistics((lambda d: d.semilength(), lambda d: d.semilength()))
    sage: for solution in sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items()))):
    ....:     print(solution)
    {[]: 0, [1, 0]: 1, [1, 0, 1, 0]: 1, [1, 1, 0, 0]: 2}
    {[]: 0, [1, 0]: 1, [1, 0, 1, 0]: 2, [1, 1, 0, 0]: 1}

The following looks correct, because `alpha = beta \circ S` forces
`S([1,0,1,0]) = [1,1,0,0]` and `s = tau \circ S` forces therefore `s([1,0,1,0])
= \tau(S([1,0,1,0])) = \tau([1,1,0,0]) = 1`::

    sage: bij = Bijectionist(A, B, tau)
    sage: bij.set_statistics((alpha, beta), (lambda d: d.semilength(), lambda d: d.semilength()))
    sage: for solution in bij.solutions_iterator():
    ....:     print(solution)
    {[]: 0, [1, 0]: 1, [1, 0, 1, 0]: 1, [1, 1, 0, 0]: 2}

Now we introduce a intertwining relation::

    sage: concat_path = lambda D1, D2: DyckWord(list(D1) + list(D2))
    sage: pi_rho = (2, concat_path, lambda x, y: x+y)

Without `\alpha` and `\beta` but with `\pi` and `\rho` the other values are
forced because `s([1,0,1,0]) = s(\pi([1,0], [1,0])) = \rho(s([1,0]), s([1,0]))
= 2`::

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
    sage: def alpha(D): return (D.area(), D.bounce())
    sage: def beta(D): return (D.bounce(), D.area())
    sage: def tau(D): return D.number_of_touch_points()

    sage: bij = Bijectionist(A, B, tau, alpha_beta=((lambda d: d.semilength(), lambda d: d.semilength()),))
    sage: for solution in sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items()))):
    ....:     print(solution)
    {[]: 0, [1, 0]: 1, [1, 0, 1, 0]: 1, [1, 1, 0, 0]: 2}
    {[]: 0, [1, 0]: 1, [1, 0, 1, 0]: 2, [1, 1, 0, 0]: 1}

Constant blocks::

    sage: A = B = 'abcd'
    sage: def pi(p1, p2): return 'abcdefgh'[A.index(p1) + A.index(p2)]
    sage: def rho(s1, s2): return (s1 + s2) % 2
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

    sage: def concat(p1, p2): return Permutation(p1 + [i + len(p1) for i in p2])

    sage: A = B = [permutation for n in range(4) for permutation in Permutations(n)]
    sage: bij = Bijectionist(A, B, Permutation.number_of_fixed_points, alpha_beta=((len, len),), pi_rho=((2, concat, lambda x, y: x + y),))
    sage: for solution in sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items()))):
    ....:     print(solution)
    {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 0, [1, 2, 3]: 3, [1, 3, 2]: 1, [2, 1, 3]: 1, [2, 3, 1]: 0, [3, 1, 2]: 0, [3, 2, 1]: 1}
    {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 0, [1, 2, 3]: 3, [1, 3, 2]: 1, [2, 1, 3]: 1, [2, 3, 1]: 0, [3, 1, 2]: 1, [3, 2, 1]: 0}
    {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 0, [1, 2, 3]: 3, [1, 3, 2]: 1, [2, 1, 3]: 1, [2, 3, 1]: 1, [3, 1, 2]: 0, [3, 2, 1]: 0}

Statistics::

    sage: N = 4; A = B = [permutation for n in range(N) for permutation in Permutations(n)]
    sage: def wex(p): return len(p.weak_excedences())
    sage: def fix(p): return len(p.fixed_points())
    sage: def des(p): return len(p.descents(final_descent=True)) if p else 0
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
    sage: alpha_beta = [(len, len)]
    sage: value_restrictions = [(Permutation([1, 2]), [1]), (Permutation([3, 2, 1]), [2, 3, 4])]
    sage: bij = Bijectionist(A, B, tau, alpha_beta=alpha_beta, value_restrictions=value_restrictions)
    sage: for sol in sorted(list(bij.solutions_iterator()), key=lambda d: tuple(sorted(d.items()))):
    ....:     print(sol)
    {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 1, [1, 3, 2]: 2, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 2, [3, 2, 1]: 3}
    ...
    {[]: 0, [1]: 1, [1, 2]: 1, [2, 1]: 2, [1, 2, 3]: 3, [1, 3, 2]: 2, [2, 1, 3]: 2, [2, 3, 1]: 2, [3, 1, 2]: 1, [3, 2, 1]: 2}

    sage: A = B = [permutation for n in range(4) for permutation in Permutations(n)]
    sage: tau = Permutation.longest_increasing_subsequence_length
    sage: bij = Bijectionist(A, B, tau, value_restrictions=((Permutation([1, 2]), [4, 5]),))
    sage: next(bij.solutions_iterator())
    Traceback (most recent call last):
    ...
    ValueError: no possible values found for singleton block [[1, 2]]
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
from collections import namedtuple, defaultdict
from sage.numerical.mip import MixedIntegerLinearProgram, MIPSolverException
from sage.rings.integer_ring import ZZ
from sage.combinat.set_partition import SetPartition, SetPartitions
from sage.sets.disjoint_set import DisjointSet
from sage.structure.sage_object import SageObject
from copy import copy
from sage.misc.verbose import get_verbose


class Bijectionist(SageObject):
    r"""
    A toolbox to list all possible bijections between two finite sets
    under various constraints.

    INPUT:

    - ``A``, ``B`` -- sets of equal size, given as a list

    - ``tau`` -- (optional) a function from ``B`` to ``Z``, in case of
      ``None``, the identity map ``lambda x: x`` is used

    - ``alpha_beta`` -- (optional) a list of pairs of statistics ``alpha`` from
      ``A`` to ``W`` and ``beta`` from ``B`` to ``W``

    - ``P`` -- (optional) a partition of ``A``

    - ``pi_rho`` -- (optional) a list of triples ``(k, pi, rho)``, where

      * ``pi`` -- a ``k``-ary operation composing objects in ``A`` and
      * ``rho`` -- a ``k``-ary function composing statistic values in ``Z``

    - ``elements_distributions`` -- (optional) a list of pairs ``(tA, tZ)``,
      specifying the distributions of ``tA``

    - ``value_restrictions`` -- (optional) a list of pairs ``(a, tZ)``,
      restricting the possible values of ``a``

    - ``solver`` -- (optional) the backend used to solve the mixed integer
      linear programs

    ``W`` and ``Z`` can be arbitrary sets.  As a natural example we may think
    of the natural numbers or tuples of integers.

    We are looking for a statistic `s: A\to Z` and a bijection `S: A\to B` such
    that

    - `s = \tau \circ S`: the statistics `s` and `\tau` are equidistributed and
      `S` is an intertwining bijection.

    - `\alpha = \beta \circ S`: the statistics `\alpha` and `\beta` are
      equidistributed and `S` is an intertwining bijection.

    - `s` is constant on the blocks of `P`.

    - `s(\pi(a_1,\dots, a_k)) = \rho(s(a_1),\dots, s(a_k))`.

    Additionally, we may require that

    - `s(a)\in Z_a` for specified sets `Z_a\subseteq Z`, and

    - `s|_{\tilde A}` has a specified distribution for specified sets `\tilde A
      \subset A`.

    If `\tau` is the identity, the two unknown functions `s` and `S` coincide.
    Although we do not exclude other bijective choices for `\tau`, they
    probably do not make sense.

    If we want that `S` is graded, i.e. if elements of `A` and `B` have a
    notion of size and `S` should preserve this size, we can add grading
    statistics as `\alpha` and `\beta`.  Since `\alpha` and `\beta` will be
    equidistributed with `S` as an intertwining bijection, `S` will then also
    be graded.

    In summary, we have the following two commutative diagrams, where `s` and
    `S` are unknown functions.

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

        If `\tau` is the identity map, the partition `P` of `A` necessarily
        consists only of singletons.

    .. NOTE::

        The order of invocation of the methods with prefix ``set``, i.e.,
        :meth:`set_statistics`, :meth:`set_intertwining_relations`,
        :meth:`set_constant_blocks`, etc., is irrelevant.  Calling any of these
        methods a second time overrides the previous specification.
    """
    def __init__(self, A, B, tau=None, alpha_beta=tuple(), P=None,
                 pi_rho=tuple(), phi_psi=tuple(), Q=None,
                 elements_distributions=tuple(),
                 value_restrictions=tuple(), solver=None, key=None):
        """
        Initialize the bijectionist.

        TESTS:

        Check that large input sets are handled well::

            sage: A = B = range(20000)
            sage: bij = Bijectionist(A, B)
        """
        # glossary of standard letters:
        # A, B, Z, W ... finite sets
        # P ... set partition of A
        # tA, tB, tZ, tP ... subsets
        # a in A, b in B, p in P
        # S: A -> B
        # alpha: A -> W, beta: B -> W
        # s: A -> Z, tau: B -> Z
        # k arity of pi and rho
        # pi: A^k -> A, rho: Z^k -> Z
        # a_tuple in A^k
        self._A = list(A)
        self._B = list(B)
        assert len(self._A) == len(set(self._A)), "A must have distinct items"
        assert len(self._B) == len(set(self._B)), "B must have distinct items"
        self._bmilp = None
        self._sorter = {}
        self._sorter["A"] = lambda x: sorted(x, key=self._A.index)
        self._sorter["B"] = lambda x: sorted(x, key=self._B.index)

        if tau is None:
            self._tau = {b: b for b in self._B}
        else:
            self._tau = {b: tau(b) for b in self._B}
        # we store Z as a list to keep an order
        self._Z = set(self._tau.values())
        if key is not None and "Z" in key:
            self._sorter["Z"] = lambda x: sorted(x, key=key["Z"])
            self._Z = self._sorter["Z"](self._Z)
        else:
            try:
                self._Z = sorted(self._Z)
                self._sorter["Z"] = lambda x: sorted(x)
            except TypeError:
                self._Z = list(self._Z)
                self._sorter["Z"] = lambda x: list(x)
        if P is None:
            P = []

        # set optional inputs
        self.set_statistics(*alpha_beta)
        self.set_value_restrictions(*value_restrictions)
        self.set_distributions(*elements_distributions)
        self.set_quadratic_relation(*phi_psi)
        self.set_homomesic(Q)
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

        - ``P`` -- set partition of `A`, singletons may be omitted

        EXAMPLES:

        Initially the partitions are set to singleton blocks.  The
        current partition can be reviewed using
        :meth:`constant_blocks`::

            sage: A = B = 'abcd'
            sage: bij = Bijectionist(A, B, lambda x: B.index(x) % 2)
            sage: bij.constant_blocks()
            {}

            sage: bij.set_constant_blocks([['a', 'c']])
            sage: bij.constant_blocks()
            {{'a', 'c'}}

        We now add a map that combines some blocks::

            sage: def pi(p1, p2): return 'abcdefgh'[A.index(p1) + A.index(p2)]
            sage: def rho(s1, s2): return (s1 + s2) % 2
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
            StopIteration
        """
        self._bmilp = None
        self._P = DisjointSet(self._A)
        P = sorted(self._sorter["A"](p) for p in P)
        for p in P:
            for a in p:
                self._P.union(p[0], a)

        self._possible_block_values = None

    def constant_blocks(self, singletons=False, optimal=False):
        r"""
        Return the set partition `P` of `A` such that `s: A\to Z` is
        known to be constant on the blocks of `P`.

        INPUT:

        - ``singletons`` -- boolean (default: ``False``); whether or not to
          include singleton blocks in the output

        - ``optimal`` -- boolean (default: ``False``); whether or not to
          compute the coarsest possible partition

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

        ALGORITHM:

        We add

        .. MATH::

            \sum_{a\in A, z\in Z} x_{p(a), z} s^z t^{\alpha(a)}
            = \sum_{b\in B} s^{\tau(b)} t(\beta(b))

        as a matrix equation to the MILP.

        EXAMPLES:

        We look for bijections `S` on permutations such that the
        number of weak exceedences of `S(\pi)` equals the number of
        descents of `\pi`, and statistics `s`, such that the number
        of fixed points of `S(\pi)` equals `s(\pi)`::

            sage: N = 4; A = B = [permutation for n in range(N) for permutation in Permutations(n)]
            sage: def wex(p): return len(p.weak_excedences())
            sage: def fix(p): return len(p.fixed_points())
            sage: def des(p): return len(p.descents(final_descent=True)) if p else 0
            sage: def adj(p): return len([e for (e, f) in zip(p, p[1:]+[0]) if e == f+1])
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
            ValueError: statistics alpha and beta are not equidistributed

        TESTS:

        Calling ``set_statistics`` without arguments should restore the previous state::

            sage: N = 3; A = B = [permutation for n in range(N) for permutation in Permutations(n)]
            sage: def wex(p): return len(p.weak_excedences())
            sage: def fix(p): return len(p.fixed_points())
            sage: def des(p): return len(p.descents(final_descent=True)) if p else 0
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
        self._bmilp = None
        self._n_statistics = len(alpha_beta)
        # TODO: do we really want to recompute statistics every time?
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
                raise ValueError(f"statistics alpha and beta do not have the same image, {v} is not a value of alpha, but of beta")
            self._statistics_fibers[v][1].append(b)

        # check compatibility
        if not all(len(fiber[0]) == len(fiber[1])
                   for fiber in self._statistics_fibers.values()):
            raise ValueError("statistics alpha and beta are not equidistributed")

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
            sage: def wex(p): return len(p.weak_excedences())
            sage: def fix(p): return len(p.fixed_points())
            sage: def des(p): return len(p.descents(final_descent=True)) if p else 0
            sage: def adj(p): return len([e for (e, f) in zip(p, p[1:]+[0]) if e == f+1])
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

        - ``header`` -- boolean (default: ``True``); whether to include a
          header with the standard Greek letters

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
            sage: def wex(p): return len(p.weak_excedences())
            sage: def fix(p): return len(p.fixed_points())
            sage: def des(p): return len(p.descents(final_descent=True)) if p else 0
            sage: def adj(p): return len([e for (e, f) in zip(p, p[1:]+[0]) if e == f+1])
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_statistics((wex, des), (fix, adj))
            sage: a, b = bij.statistics_table()
            sage: table(a, header_row=True, frame=True)
            ┌───────────┬────────┬────────┐
            │ a         │ α_1(a) │ α_2(a) │
            ╞═══════════╪════════╪════════╡
            │ []        │ 0      │ 0      │
            ├───────────┼────────┼────────┤
            │ [1]       │ 1      │ 1      │
            ├───────────┼────────┼────────┤
            │ [1, 2]    │ 2      │ 2      │
            ├───────────┼────────┼────────┤
            │ [2, 1]    │ 1      │ 0      │
            ├───────────┼────────┼────────┤
            │ [1, 2, 3] │ 3      │ 3      │
            ├───────────┼────────┼────────┤
            │ [1, 3, 2] │ 2      │ 1      │
            ├───────────┼────────┼────────┤
            │ [2, 1, 3] │ 2      │ 1      │
            ├───────────┼────────┼────────┤
            │ [2, 3, 1] │ 2      │ 0      │
            ├───────────┼────────┼────────┤
            │ [3, 1, 2] │ 1      │ 0      │
            ├───────────┼────────┼────────┤
            │ [3, 2, 1] │ 2      │ 1      │
            └───────────┴────────┴────────┘
            sage: table(b, header_row=True, frame=True)
            ┌───────────┬───┬────────┬────────┐
            │ b         │ τ │ β_1(b) │ β_2(b) │
            ╞═══════════╪═══╪════════╪════════╡
            │ []        │ 0 │ 0      │ 0      │
            ├───────────┼───┼────────┼────────┤
            │ [1]       │ 1 │ 1      │ 1      │
            ├───────────┼───┼────────┼────────┤
            │ [1, 2]    │ 2 │ 1      │ 0      │
            ├───────────┼───┼────────┼────────┤
            │ [2, 1]    │ 1 │ 2      │ 2      │
            ├───────────┼───┼────────┼────────┤
            │ [1, 2, 3] │ 3 │ 1      │ 0      │
            ├───────────┼───┼────────┼────────┤
            │ [1, 3, 2] │ 2 │ 2      │ 1      │
            ├───────────┼───┼────────┼────────┤
            │ [2, 1, 3] │ 2 │ 2      │ 1      │
            ├───────────┼───┼────────┼────────┤
            │ [2, 3, 1] │ 2 │ 2      │ 1      │
            ├───────────┼───┼────────┼────────┤
            │ [3, 1, 2] │ 2 │ 2      │ 0      │
            ├───────────┼───┼────────┼────────┤
            │ [3, 2, 1] │ 1 │ 3      │ 3      │
            └───────────┴───┴────────┴────────┘

        TESTS:

        If no statistics are given, the table should still be able to be generated::

            sage: A = B = [permutation for n in range(3) for permutation in Permutations(n)]
            sage: tau = Permutation.longest_increasing_subsequence_length
            sage: bij = Bijectionist(A, B, tau)
            sage: a, b = bij.statistics_table()
            sage: table(a, header_row=True, frame=True)
            ┌────────┐
            │ a      │
            ╞════════╡
            │ []     │
            ├────────┤
            │ [1]    │
            ├────────┤
            │ [1, 2] │
            ├────────┤
            │ [2, 1] │
            └────────┘
            sage: table(b, header_row=True, frame=True)
            ┌────────┬───┐
            │ b      │ τ │
            ╞════════╪═══╡
            │ []     │ 0 │
            ├────────┼───┤
            │ [1]    │ 1 │
            ├────────┼───┤
            │ [1, 2] │ 2 │
            ├────────┼───┤
            │ [2, 1] │ 1 │
            └────────┴───┘

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

    def set_value_restrictions(self, *value_restrictions):
        r"""
        Restrict the set of possible values `s(a)` for a given element
        `a`.

        .. WARNING::

             Any restriction imposed by a previous invocation of
             :meth:`set_value_restrictions` will be overwritten!

        INPUT:

        - ``value_restrictions`` -- one or more pairs `(a\in A, \tilde
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
            sage: for sol in sorted(bij.solutions_iterator(), key=lambda d: sorted(d.items())):
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

        TESTS::

            sage: A = B = [permutation for n in range(4) for permutation in Permutations(n)]
            sage: tau = Permutation.longest_increasing_subsequence_length
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_value_restrictions((Permutation([1, 2]), [4, 5]))
            sage: bij._compute_possible_block_values()
            Traceback (most recent call last):
            ...
            ValueError: no possible values found for singleton block [[1, 2]]

            sage: A = B = [permutation for n in range(4) for permutation in Permutations(n)]
            sage: tau = Permutation.longest_increasing_subsequence_length
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_constant_blocks([[permutation for permutation in Permutations(n)] for n in range(4)])
            sage: bij.set_value_restrictions((Permutation([1, 2]), [4, 5]))
            sage: bij._compute_possible_block_values()
            Traceback (most recent call last):
            ...
            ValueError: no possible values found for block [[1, 2], [2, 1]]

            sage: A = B = [permutation for n in range(4) for permutation in Permutations(n)]
            sage: tau = Permutation.longest_increasing_subsequence_length
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_value_restrictions(((1, 2), [4, 5, 6]))
            Traceback (most recent call last):
            ...
            AssertionError: element (1, 2) was not found in A
        """
        # it might be much cheaper to construct the sets as subsets
        # of _statistics_possible_values - however, we do not want to
        # insist that set_value_restrictions is called after
        # set_statistics
        self._bmilp = None
        set_Z = set(self._Z)
        self._restrictions_possible_values = {a: set_Z for a in self._A}
        for a, values in value_restrictions:
            assert a in self._A, f"element {a} was not found in A"
            self._restrictions_possible_values[a] = self._restrictions_possible_values[a].intersection(values)

    def _compute_possible_block_values(self):
        r"""
        Update the dictionary of possible values of each block.

        This has to be called whenever ``self._P`` was modified.

        It raises a :exc:`ValueError`, if the restrictions on a
        block are contradictory.

        TESTS::

            sage: A = B = [permutation for n in range(4) for permutation in Permutations(n)]
            sage: tau = Permutation.longest_increasing_subsequence_length
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_value_restrictions((Permutation([1, 2]), [4, 5]))
            sage: bij._compute_possible_block_values()
            Traceback (most recent call last):
            ...
            ValueError: no possible values found for singleton block [[1, 2]]
        """
        self._possible_block_values = {}  # P -> Power(Z)
        for p, block in self._P.root_to_elements_dict().items():
            sets = ([self._restrictions_possible_values[a] for a in block]
                    + [self._statistics_possible_values[a] for a in block])
            self._possible_block_values[p] = _non_copying_intersection(sets)
            if not self._possible_block_values[p]:
                if len(block) == 1:
                    raise ValueError(f"no possible values found for singleton block {block}")
                raise ValueError(f"no possible values found for block {block}")

    def set_distributions(self, *elements_distributions):
        r"""
        Specify the distribution of `s` for a subset of elements.

        .. WARNING::

             Any restriction imposed by a previous invocation of
             :meth:`set_distributions` will be overwritten!

        INPUT:

        - one or more pairs of `(\tilde A, \tilde Z)`, where `\tilde
          A\subseteq A` and `\tilde Z` is a list of values in `Z` of
          the same size as `\tilde A`

        This method specifies that `\{s(a) | a\in\tilde A\}` equals
        `\tilde Z` as a multiset for each of the pairs.

        When specifying several distributions, the subsets of `A` do
        not have to be disjoint.

        ALGORITHM:

        We add

        .. MATH::

            \sum_{a\in\tilde A} x_{p(a), z}t^z = \sum_{z\in\tilde Z} t^z,

        where `p(a)` is the block containing `a`, for each given
        distribution as a vector equation to the MILP.

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
             ([[2, 3, 1]], [2]),
             ([[1, 2], [2, 1]], [1, 2])]

        TESTS:

        Because of the current implementation of the output calculation, we do
        not improve our solution if we do not gain any unique solutions::

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
            sage: def alpha(p): return p(1) if len(p) > 0 else 0
            sage: def beta(p): return p(1) if len(p) > 0 else 0
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
            ValueError: element [1, 2, 3, 4] was not found in A
            sage: bij.set_distributions(([Permutation([1, 2, 3])], [-1]))
            Traceback (most recent call last):
            ...
            ValueError: value -1 was not found in tau(A)

        Note that the same error occurs when an element that is not the first element of the list is
        not in `A`.
        """
        self._bmilp = None
        for tA, tZ in elements_distributions:
            assert len(tA) == len(tZ), f"{tA} and {tZ} are not of the same size!"
            for a, z in zip(tA, tZ):
                if a not in self._A:
                    raise ValueError(f"element {a} was not found in A")
                if z not in self._Z:
                    raise ValueError(f"value {z} was not found in tau(A)")
        self._elements_distributions = tuple(elements_distributions)

    def set_intertwining_relations(self, *pi_rho):
        r"""
        Add restrictions of the form `s(\pi(a_1,\dots, a_k)) =
        \rho(s(a_1),\dots, s(a_k))`.

        .. WARNING::

             Any restriction imposed by a previous invocation of
             :meth:`set_intertwining_relations` will be overwritten!

        INPUT:

        - ``pi_rho`` -- one or more tuples `(k, \pi: A^k\to A, \rho:
          Z^k\to Z, \tilde A)` where `\tilde A` (optional) is a
          `k`-ary function that returns true if and only if a
          `k`-tuple of objects in `A` is in the domain of `\pi`

        ALGORITHM:

        The relation

        .. MATH::

            s(\pi(a_1,\dots, a_k)) = \rho(s(a_1),\dots, s(a_k))

        for each pair `(\pi, \rho)` implies immediately that
        `s(\pi(a_1,\dots, a_k))` only depends on the blocks of
        `a_1,\dots, a_k`.

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

        Note that `z` must be a possible value of `p` and each `z_i`
        must be a possible value of `p_i`.

        EXAMPLES:

        We can concatenate two permutations by increasing the values
        of the second permutation by the length of the first
        permutation::

            sage: def concat(p1, p2): return Permutation(p1 + [i + len(p1) for i in p2])

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
        self._bmilp = None
        Pi_Rho = namedtuple("Pi_Rho", "numargs pi rho domain")
        self._pi_rho = []

        for pi_rho_tuple in pi_rho:
            if len(pi_rho_tuple) == 3:
                k, pi, rho = pi_rho_tuple
                domain = None
            else:
                k, pi, rho, domain = pi_rho_tuple

            self._pi_rho.append(Pi_Rho(numargs=k, pi=pi, rho=rho, domain=domain))

    set_semi_conjugacy = set_intertwining_relations

    def set_quadratic_relation(self, *phi_psi):
        r"""
        Add restrictions of the form `s\circ\psi\circ s = \phi`.

        INPUT:

        - ``phi_psi`` -- (optional) a list of pairs `(\phi, \rho)` where `\phi:
          A\to Z` and `\psi: Z\to A`

        ALGORITHM:

        We add

        .. MATH::

            x_{p(a), z} = x_{p(\psi(z)), \phi(a)}

        for `a\in A` and `z\in Z` to the MILP, where `\phi:A\to Z`
        and `\psi:Z\to A`.  Note that, in particular, `\phi` must be
        constant on blocks.

        EXAMPLES::

            sage: A = B = DyckWords(3)
            sage: bij = Bijectionist(A, B)
            sage: bij.set_statistics((lambda D: D.number_of_touch_points(), lambda D: D.number_of_initial_rises()))
            sage: ascii_art(sorted(bij.minimal_subdistributions_iterator()))
            [ (             [   /\   ] )
            [ (             [  /  \  ] )  ( [    /\    /\    ]  [  /\      /\/\  ] )
            [ ( [ /\/\/\ ], [ /    \ ] ), ( [ /\/  \, /  \/\ ], [ /  \/\, /    \ ] ),
            <BLANKLINE>
             ( [           /\   ]                     ) ]
             ( [  /\/\    /  \  ]  [            /\  ] ) ]
             ( [ /    \, /    \ ], [ /\/\/\, /\/  \ ] ) ]
            sage: bij.set_quadratic_relation((lambda D: D, lambda D: D))
            sage: ascii_art(sorted(bij.minimal_subdistributions_iterator()))
            [ (             [   /\   ] )
            [ (             [  /  \  ] )  ( [    /\  ]  [  /\/\  ] )
            [ ( [ /\/\/\ ], [ /    \ ] ), ( [ /\/  \ ], [ /    \ ] ),
            <BLANKLINE>
            <BLANKLINE>
             ( [  /\    ]  [  /\    ] )  ( [  /\/\  ]  [    /\  ] )
             ( [ /  \/\ ], [ /  \/\ ] ), ( [ /    \ ], [ /\/  \ ] ),
            <BLANKLINE>
             ( [   /\   ]             ) ]
             ( [  /  \  ]             ) ]
             ( [ /    \ ], [ /\/\/\ ] ) ]
        """
        self._bmilp = None
        self._phi_psi = phi_psi

    def set_homomesic(self, Q):
        """
        Assert that the average of `s` on each block of `Q` is
        constant.

        INPUT:

        - ``Q`` -- set partition of ``A``

        EXAMPLES::

            sage: A = B = [1,2,3]
            sage: bij = Bijectionist(A, B, lambda b: b % 3)
            sage: bij.set_homomesic([[1,2], [3]])
            sage: list(bij.solutions_iterator())
            [{1: 2, 2: 0, 3: 1}, {1: 0, 2: 2, 3: 1}]
        """
        self._bmilp = None
        if Q is None:
            self._Q = None
        else:
            self._Q = SetPartition(Q)
            assert self._Q in SetPartitions(self._A), f"{Q} must be a set partition of A"

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
            sage: bij.constant_blocks(optimal=True)                             # indirect doctest
            {{[], [1], [1, 2], [2, 1]}}

        In this other example we look at permutations with length 2 and 3::

            sage: N = 4
            sage: A = B = [permutation for n in range(2, N) for permutation in Permutations(n)]
            sage: def tau(p): return p[0] if len(p) else 0
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

            sage: def concat(p1, p2): return Permutation(p1 + [i + len(p1) for i in p2])
            sage: def union(p1, p2): return Partition(sorted(list(p1) + list(p2), reverse=True))
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
            sage: def alpha1(p): return len(p.weak_excedences())
            sage: def alpha2(p): return len(p.fixed_points())
            sage: def beta1(p): return len(p.descents(final_descent=True)) if p else 0
            sage: def beta2(p): return len([e for (e, f) in zip(p, p[1:] + [0]) if e == f + 1])
            sage: tau = Permutation.longest_increasing_subsequence_length
            sage: def rotate_permutation(p):
            ....:    cycle = Permutation(tuple(range(1, len(p) + 1)))
            ....:    return Permutation([cycle.inverse()(p(cycle(i))) for i in range(1, len(p) + 1)])
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_statistics((alpha1, beta1), (alpha2, beta2))
            sage: from sage.combinat.cyclic_sieving_phenomenon import orbit_decomposition
            sage: bij.set_constant_blocks(orbit_decomposition(A, rotate_permutation))
            sage: P = bij.constant_blocks()
            sage: P = [sorted(p, key=lambda p: (len(p), p)) for p in P]
            sage: P = sorted(P, key=lambda p: (len(next(iter(p))), len(p)))
            sage: for p in P:
            ....:     print(p)
            [[1, 3, 2], [2, 1, 3], [3, 2, 1]]
            [[1, 4, 3, 2], [3, 2, 1, 4]]
            [[2, 1, 4, 3], [4, 3, 2, 1]]
            [[1, 2, 4, 3], [1, 3, 2, 4], [2, 1, 3, 4], [4, 2, 3, 1]]
            [[1, 3, 4, 2], [2, 3, 1, 4], [2, 4, 3, 1], [3, 2, 4, 1]]
            [[1, 4, 2, 3], [3, 1, 2, 4], [4, 1, 3, 2], [4, 2, 1, 3]]
            [[2, 4, 1, 3], [3, 1, 4, 2], [3, 4, 2, 1], [4, 3, 1, 2]]

            sage: P = bij.constant_blocks(optimal=True)
            sage: P = [sorted(p, key=lambda p: (len(p), p)) for p in P]
            sage: P = sorted(P, key=lambda p: (len(next(iter(p))), len(p)))
            sage: for p in P:
            ....:     print(p)
            [[1], [1, 2], [1, 2, 3], [1, 2, 3, 4]]
            [[1, 3, 2], [2, 1, 3], [3, 2, 1],
             [1, 2, 4, 3], [1, 3, 2, 4], [1, 3, 4, 2], [1, 4, 3, 2],
             [2, 1, 3, 4], [2, 1, 4, 3], [2, 3, 1, 4], [2, 3, 4, 1],
             [2, 4, 3, 1], [3, 2, 1, 4], [3, 2, 4, 1], [4, 2, 3, 1],
             [4, 3, 2, 1]]
            [[1, 4, 2, 3], [2, 4, 1, 3], [3, 1, 2, 4], [3, 1, 4, 2],
             [3, 4, 2, 1], [4, 1, 3, 2], [4, 2, 1, 3], [4, 3, 1, 2]]

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
        if self._bmilp is None:
            self._bmilp = _BijectionistMILP(self)

        solution = next(self._bmilp.solutions_iterator(True, []))
        # multiple_preimages[tZ] are the blocks p which have the same
        # value tZ[i] in the i-th known solution
        multiple_preimages = {(z,): tP
                              for z, tP in _invert_dict(solution).items()
                              if len(tP) > 1}

        # _P has to be copied to not mess with the solution process
        # since we do not want to regenerate the bmilp in each step,
        # so blocks have to stay consistent during the whole process
        tmp_P = copy(self._P)

        # check whether blocks p1 and p2 can have different values,
        # if so return such a solution
        def different_values(p1, p2):
            tmp_constraints = [self._bmilp._x[p1, z] + self._bmilp._x[p2, z] <= 1
                               for z in self._possible_block_values[p1]
                               if z in self._possible_block_values[p2]]
            return next(self._bmilp.solutions_iterator(True, tmp_constraints))

        # try to find a pair of blocks having the same value on all
        # known solutions, and a solution such that the values are
        # different on this solution
        def merge_until_split():
            for tZ in list(multiple_preimages):
                tP = multiple_preimages[tZ]
                for i2 in range(len(tP) - 1, -1, -1):
                    for i1 in range(i2):
                        try:
                            solution = different_values(tP[i1], tP[i2])
                        except StopIteration:
                            tmp_P.union(tP[i1], tP[i2])
                            if len(multiple_preimages[tZ]) == 2:
                                del multiple_preimages[tZ]
                            else:
                                tP.remove(tP[i2])
                            break  # skip all pairs (i, j) containing i2
                        return solution

        while True:
            solution = merge_until_split()
            if solution is None:
                self._P = tmp_P
                # recreate the MILP
                self._bmilp = _BijectionistMILP(self,
                                                self._bmilp._solution_cache)
                return

            updated_multiple_preimages = defaultdict(list)
            for tZ, tP in multiple_preimages.items():
                for p in tP:
                    updated_multiple_preimages[tZ + (solution[p],)].append(p)
            multiple_preimages = updated_multiple_preimages

    def possible_values(self, p=None, optimal=False):
        r"""
        Return for each block the values of `s` compatible with the
        imposed restrictions.

        INPUT:

        - ``p`` -- (optional) a block of `P`, or an element of a
          block of `P`, or a list of these

        - ``optimal`` -- boolean (default: ``False``); whether or not to
          compute the minimal possible set of statistic values

        .. NOTE::

            Computing the minimal possible set of statistic values
            may be computationally expensive.

        .. TODO::

            currently, calling this method with ``optimal=True`` does
            not update the internal dictionary, because this would
            interfere with the variables of the MILP.

        EXAMPLES::

            sage: A = B = ["a", "b", "c", "d"]
            sage: tau = {"a": 1, "b": 1, "c": 1, "d": 2}.get
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_constant_blocks([["a", "b"]])
            sage: bij.possible_values(A)
            {'a': {1, 2}, 'b': {1, 2}, 'c': {1, 2}, 'd': {1, 2}}
            sage: bij.possible_values(A, optimal=True)
            {'a': {1}, 'b': {1}, 'c': {1, 2}, 'd': {1, 2}}

        The internal dictionary is not updated::

            sage: bij.possible_values(A)
            {'a': {1, 2}, 'b': {1, 2}, 'c': {1, 2}, 'd': {1, 2}}

        TESTS::

            sage: A = B = ["a", "b", "c", "d"]
            sage: tau = {"a": 1, "b": 1, "c": 2, "d": 2}.get
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_constant_blocks([["a", "b"]])

        Test if all formats are really possible::

            sage: bij.possible_values(p='a')
            {'a': {1, 2}, 'b': {1, 2}}
            sage: bij.possible_values(p=["a", "b"])
            {'a': {1, 2}, 'b': {1, 2}}
            sage: bij.possible_values(p=[["a", "b"]])
            {'a': {1, 2}, 'b': {1, 2}}
            sage: bij.possible_values(p=[["a", "b"], ["c"]])
            {'a': {1, 2}, 'b': {1, 2}, 'c': {1, 2}}

        Test an unfeasible problem::

            sage: A = B = 'ab'
            sage: bij = Bijectionist(A, B, lambda x: B.index(x) % 2)
            sage: bij.set_constant_blocks([['a', 'b']])
            sage: bij.possible_values(p='a')
            {'a': {0, 1}, 'b': {0, 1}}
            sage: bij.possible_values(p='a', optimal=True)
            {'a': set(), 'b': set()}
        """
        # convert input to set of block representatives
        blocks = set()
        if p in self._A:
            blocks.add(self._P.find(p))
        elif isinstance(p, list):  # TODO: this looks very brittle
            for p1 in p:
                if p1 in self._A:
                    blocks.add(self._P.find(p1))
                elif isinstance(p1, list):
                    for p2 in p1:
                        blocks.add(self._P.find(p2))

        if optimal:
            if self._bmilp is None:
                self._bmilp = _BijectionistMILP(self)
            bmilp = self._bmilp
            solutions = defaultdict(set)
            try:
                solution = next(bmilp.solutions_iterator(True, []))
            except StopIteration:
                pass
            else:
                for p, z in solution.items():
                    solutions[p].add(z)
                for p in blocks:
                    tmp_constraints = [bmilp._x[p, z] == 0 for z in solutions[p]]
                    while True:
                        try:
                            solution = next(bmilp.solutions_iterator(True, tmp_constraints))
                        except StopIteration:
                            break
                        for p0, z in solution.items():
                            solutions[p0].add(z)
                        # veto new value and try again
                        tmp_constraints.append(bmilp._x[p, solution[p]] == 0)

            # create dictionary to return
            possible_values = {}
            for p in blocks:
                for a in self._P.root_to_elements_dict()[p]:
                    possible_values[a] = solutions[p]
        else:
            # create dictionary to return
            if self._possible_block_values is None:
                self._compute_possible_block_values()
            possible_values = {}
            for p in blocks:
                for a in self._P.root_to_elements_dict()[p]:
                    possible_values[a] = self._possible_block_values[p]

        return possible_values

    def minimal_subdistributions_iterator(self):
        r"""
        Return all minimal subsets `\tilde A` of `A`
        together with submultisets `\tilde Z` with `s(\tilde A) =
        \tilde Z` as multisets.

        EXAMPLES::

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
            sage: def tau(D): return D.number_of_touch_points()
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

        if self._bmilp is None:
            self._bmilp = _BijectionistMILP(self)
        s = next(self._bmilp.solutions_iterator(False, []))
        while True:
            for v in self._Z:
                minimal_subdistribution.add_constraint(sum(D[a] for a in self._A if s[a] == v) == V[v])
            try:
                minimal_subdistribution.solve()
            except MIPSolverException:
                return
            d = minimal_subdistribution.get_values(D, convert=bool, tolerance=0.1)  # a dict from A to {0, 1}
            new_s = self._find_counterexample(self._A, s, d, False)
            if new_s is None:
                values = self._sorter["Z"](s[a] for a in self._A if d[a])
                yield ([a for a in self._A if d[a]], values)

                # get all variables with value 1
                active_vars = [D[a] for a in self._A
                               if minimal_subdistribution.get_values(D[a], convert=bool, tolerance=0.1)]

                # add constraint that not all of these can be 1, thus vetoing
                # the current solution
                minimal_subdistribution.add_constraint(sum(active_vars) <= len(active_vars) - 1,
                                                       name='veto')
            else:
                s = new_s

    def _find_counterexample(self, P, s0, d, on_blocks):
        r"""
        Return a solution `s` such that ``d`` is not a subdistribution of
        `s0`.

        INPUT:

        - ``P`` -- the representatives of the blocks, or `A` if
          ``on_blocks`` is ``False``

        - ``s0`` -- a solution

        - ``d`` -- a subset of `A`, in the form of a dict from `A` to
          `\{0, 1\}`

        - ``on_blocks`` -- whether to return the counterexample on
          blocks or on elements

        EXAMPLES::

            sage: A = B = ["a", "b", "c", "d", "e"]
            sage: tau = {"a": 1, "b": 1, "c": 2, "d": 2, "e": 3}.get
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_constant_blocks([["a", "b"]])
            sage: bij.set_value_restrictions(("a", [1, 2]))
            sage: next(bij.solutions_iterator())
            {'a': 1, 'b': 1, 'c': 2, 'd': 3, 'e': 2}

            sage: s0 = {'a': 1, 'b': 1, 'c': 2, 'd': 3, 'e': 2}
            sage: d = {'a': 1, 'b': 0, 'c': 0, 'd': 0, 'e': 0}
            sage: bij._find_counterexample(bij._A, s0, d, False)
            {'a': 2, 'b': 2, 'c': 1, 'd': 3, 'e': 1}
        """
        bmilp = self._bmilp
        for z in self._Z:
            z_in_d_count = sum(d[p] for p in P if s0[p] == z)
            if not z_in_d_count:
                continue

            # try to find a solution which has a different
            # subdistribution on d than s0
            z_in_d = sum(d[p] * bmilp._x[self._P.find(p), z]
                         for p in P
                         if z in self._possible_block_values[self._P.find(p)])

            # it is sufficient to require that z occurs less often as
            # a value among {a | d[a] == 1} than it does in
            # z_in_d_count, because, if the distributions are
            # different, one such z must exist
            tmp_constraints = [z_in_d <= z_in_d_count - 1]
            try:
                solution = next(bmilp.solutions_iterator(on_blocks, tmp_constraints))
                return solution
            except StopIteration:
                pass

    def minimal_subdistributions_blocks_iterator(self):
        r"""
        Return all representatives of minimal subsets `\tilde P`
        of `P` together with submultisets `\tilde Z`
        with `s(\tilde P) = \tilde Z` as multisets.

        .. WARNING::

            If there are several solutions with the same support
            (i.e., the sets of block representatives are the same),
            only one of these will be found, even if the
            distributions are different, see the doctest below.  To
            find all solutions, use
            :meth:`minimal_subdistributions_iterator`, which is,
            however, computationally more expensive.

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
            sage: def tau(D): return D.number_of_touch_points()
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
            [(['b', 'b', 'c', 'd', 'e'], [1, 1, 2, 2, 3])]

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

            sage: sorted(bij.minimal_subdistributions_blocks_iterator())        # random
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
            minimal_subdistribution.add_constraint(X[p] * len(self._P.root_to_elements_dict()[p]) >= D[p] >= X[p])

        def add_counter_example_constraint(s):
            for v in self._Z:
                minimal_subdistribution.add_constraint(sum(D[p] for p in P
                                                           if s[p] == v) == V[v])

        if self._bmilp is None:
            self._bmilp = _BijectionistMILP(self)

        s = next(self._bmilp.solutions_iterator(True, []))
        add_counter_example_constraint(s)
        while True:
            try:
                minimal_subdistribution.solve()
            except MIPSolverException:
                return
            d = minimal_subdistribution.get_values(D, convert=ZZ, tolerance=0.1)  # a dict from P to multiplicities
            new_s = self._find_counterexample(P, s, d, True)
            if new_s is None:
                yield ([p for p in P for _ in range(ZZ(d[p]))],
                       self._sorter["Z"](s[p]
                                         for p in P
                                         for _ in range(ZZ(d[p]))))

                support = [X[p] for p in P if d[p]]
                # add constraint that the support is different
                minimal_subdistribution.add_constraint(sum(support) <= len(support) - 1,
                                                       name='veto')
            else:
                s = new_s
                add_counter_example_constraint(s)

    def _preprocess_intertwining_relations(self):
        r"""
        Make ``self._P`` be the finest set partition coarser
        than ``self._P`` such that composing elements preserves
        blocks.

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

        In particular, if `P` consists only if singletons, this
        method has no effect.

        .. TODO::

            it is not clear whether this method makes sense

        EXAMPLES::

            sage: A = B = 'abcd'
            sage: bij = Bijectionist(A, B, lambda x: B.index(x) % 2)
            sage: def pi(p1, p2): return 'abcdefgh'[A.index(p1) + A.index(p2)]
            sage: def rho(s1, s2): return (s1 + s2) % 2
            sage: bij.set_intertwining_relations((2, pi, rho))
            sage: bij._preprocess_intertwining_relations()
            sage: bij._P
            {{'a'}, {'b'}, {'c'}, {'d'}}

        However, adding that ``'a'`` and ``'c'`` are in the same block,
        we can infer that also ``'b'`` and ``'d'`` are in the same
        block::

            sage: bij.set_constant_blocks([['a', 'c']])
            sage: bij._P
            {{'a', 'c'}, {'b'}, {'d'}}
            sage: bij._preprocess_intertwining_relations()
            sage: bij._P
            {{'a', 'c'}, {'b', 'd'}}

        Let a group act on permutations::

            sage: A = B = Permutations(3)
            sage: bij = Bijectionist(A, B, lambda x: x[0])
            sage: bij.set_intertwining_relations((1, lambda pi: pi.reverse(), lambda z: z))
            sage: bij._preprocess_intertwining_relations()
            sage: bij._P
            {{[1, 2, 3]}, {[1, 3, 2]}, {[2, 1, 3]}, {[2, 3, 1]}, {[3, 1, 2]}, {[3, 2, 1]}}

        Thus, in this case we do not detect the constant blocks::

           sage: bij.constant_blocks(optimal=True)
           {{[1, 2, 3], [3, 2, 1]}, {[1, 3, 2], [2, 3, 1]}, {[2, 1, 3], [3, 1, 2]}}
        """
        A = self._A
        P = self._P
        images = defaultdict(set)  # A^k -> A, a_1,...,a_k +-> {pi(a_1,...,a_k) for all pi}
        for pi_rho in self._pi_rho:
            for a_tuple in itertools.product(*([A] * pi_rho.numargs)):
                if pi_rho.domain is not None and not pi_rho.domain(*a_tuple):
                    continue
                a = pi_rho.pi(*a_tuple)
                if a in A:
                    images[a_tuple].add(a)

        # merge blocks
        something_changed = True
        while something_changed:
            something_changed = False
            # collect (preimage, image) pairs by (representatives) of
            # the blocks of the elements of the preimage
            updated_images = defaultdict(set)  # (p_1,...,p_k) to {a_1,....}
            for a_tuple, image_set in images.items():
                representatives = tuple(P.find(a) for a in a_tuple)
                updated_images[representatives].update(image_set)

            # merge blocks
            for a_tuple, image_set in updated_images.items():
                image = image_set.pop()
                while image_set:
                    P.union(image, image_set.pop())
                    something_changed = True
                # we keep a representative
                image_set.add(image)

            images = updated_images

    def solutions_iterator(self):
        r"""
        An iterator over all solutions of the problem.

        OUTPUT: an iterator over all possible mappings `s: A\to Z`

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
          `\tilde A` and a distribution of values given by integers
          `d_z` representing the multiplicity of each `z \in Z`, and
          `r_p = |p \cap\tilde A|` indicating the relative size of
          block `p` in the set of elements of the distribution,

        .. MATH::

            \sum_p r_p x_{p, z} = d_z.

        EXAMPLES::

            sage: A = B = 'abc'
            sage: bij = Bijectionist(A, B, lambda x: B.index(x) % 2, solver='GLPK')
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

            sage: bij = Bijectionist(A, B, tau, solver='GLPK')
            sage: bij.set_statistics((len, len))
            sage: bij.set_constant_blocks(P)
            sage: for solution in bij.solutions_iterator():
            ....:     print(solution)
            {[]: 0, [1]: 0, [1, 2]: 1, [2, 1]: 0, [1, 2, 3]: 0, [1, 3, 2]: 1, [2, 1, 3]: 1, [3, 2, 1]: 1, [2, 3, 1]: 2, [3, 1, 2]: 2}
            {[]: 0, [1]: 0, [1, 2]: 0, [2, 1]: 1, [1, 2, 3]: 0, [1, 3, 2]: 1, [2, 1, 3]: 1, [3, 2, 1]: 1, [2, 3, 1]: 2, [3, 1, 2]: 2}

        Changing or re-setting problem parameters clears the internal
        cache.  Setting the verbosity prints the MILP which is solved.::

            sage: set_verbose(2)
            sage: bij.set_constant_blocks(P)
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
                statistics: 0 <= <= 0
                statistics: 0 <= <= 0
                statistics: 1 <= x_1 <= 1
                statistics: 0 <= <= 0
                statistics: 0 <= <= 0
                statistics: 1 <= x_2 + x_4 <= 1
                statistics: 1 <= x_3 + x_5 <= 1
                statistics: 0 <= <= 0
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
                statistics: 0 <= <= 0
                statistics: 0 <= <= 0
                statistics: 1 <= x_1 <= 1
                statistics: 0 <= <= 0
                statistics: 0 <= <= 0
                statistics: 1 <= x_2 + x_4 <= 1
                statistics: 1 <= x_3 + x_5 <= 1
                statistics: 0 <= <= 0
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
                statistics: 0 <= <= 0
                statistics: 0 <= <= 0
                statistics: 1 <= x_1 <= 1
                statistics: 0 <= <= 0
                statistics: 0 <= <= 0
                statistics: 1 <= x_2 + x_4 <= 1
                statistics: 1 <= x_3 + x_5 <= 1
                statistics: 0 <= <= 0
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

        Testing interactions between multiple instances using Fedor Petrov's example from https://mathoverflow.net/q/424187::

            sage: A = B = ["a"+str(i) for i in range(1, 9)] + ["b"+str(i) for i in range(3, 9)] + ["d"]
            sage: tau = {b: 0 if i < 10 else 1 for i, b in enumerate(B)}.get
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_constant_blocks([["a"+str(i), "b"+str(i)] for i in range(1, 9) if "b"+str(i) in A])
            sage: d = [0]*8+[1]*4
            sage: bij.set_distributions((A[:8] + A[8+2:-1], d), (A[:8] + A[8:-3], d))
            sage: iterator1 = bij.solutions_iterator()
            sage: iterator2 = bij.solutions_iterator()

        Generate a solution in iterator1, iterator2 should generate the same solution and vice versa::

            sage: s1_1 = next(iterator1)
            sage: s2_1 = next(iterator2)
            sage: s1_1 == s2_1
            True
            sage: s2_2 = next(iterator2)
            sage: s1_2 = next(iterator1)
            sage: s1_2 == s2_2
            True

        Re-setting the distribution resets the cache, so a new
        iterator will generate the first solutions again, but the old
        iterator continues::

            sage: bij.set_distributions((A[:8] + A[8+2:-1], d), (A[:8] + A[8:-3], d))
            sage: iterator3 = bij.solutions_iterator()

            sage: s3_1 = next(iterator3)
            sage: s1_1 == s3_1
            True

            sage: s1_3 = next(iterator1)
            sage: len(set([tuple(sorted(s.items())) for s in [s1_1, s1_2, s1_3]]))
            3
        """
        if self._bmilp is None:
            self._bmilp = _BijectionistMILP(self)
        yield from self._bmilp.solutions_iterator(False, [])


class _BijectionistMILP:
    r"""
    Wrapper class for the MixedIntegerLinearProgram (MILP).

    This class is used to manage the MILP, add constraints, solve the
    problem and check for uniqueness of solution values.
    """
    def __init__(self, bijectionist: Bijectionist, solutions=None):
        r"""
        Initialize the mixed integer linear program.

        INPUT:

        - ``bijectionist`` -- an instance of :class:`Bijectionist`

        - ``solutions`` -- (optional) a list of solutions of the
          problem, each provided as a dictionary mapping `(a, z)` to
          a boolean, such that at least one element from each block
          of `P` appears as `a`.

        .. TODO::

            it might be cleaner not to pass the full bijectionist
            instance, but only those attributes we actually use

        TESTS::

            sage: A = B = ["a", "b", "c", "d"]
            sage: bij = Bijectionist(A, B)
            sage: from sage.combinat.bijectionist import _BijectionistMILP
            sage: _BijectionistMILP(bij)
            <sage.combinat.bijectionist._BijectionistMILP object at ...>
        """
        # the attributes of the bijectionist class we actually use:
        # _possible_block_values
        # _elements_distributions
        # _W, _Z, _A, _B, _P, _alpha, _beta, _tau, _pi_rho, _phi_psi
        self._bijectionist = bijectionist
        # the variables of the MILP are indexed by pairs (p, z), for
        # p in _P and z an element of _posible_block_values[p].
        # Thus, _P and _posible_block_values have to be fixed before
        # creating the MILP.
        bijectionist._preprocess_intertwining_relations()
        bijectionist._compute_possible_block_values()

        self.milp = MixedIntegerLinearProgram(solver=bijectionist._solver)
        self.milp.set_objective(None)
        indices = [(p, z)
                   for p, tZ in bijectionist._possible_block_values.items()
                   for z in tZ]
        self._x = self.milp.new_variable(binary=True, indices=indices)

        tZ = bijectionist._possible_block_values
        P = bijectionist._P
        for p in _disjoint_set_roots(P):
            self.milp.add_constraint(sum(self._x[p, z] for z in tZ[p]) == 1,
                                     name=f"block {p}"[:50])
        self.add_alpha_beta_constraints()
        self.add_distribution_constraints()
        self.add_quadratic_relation_constraints()
        self.add_homomesic_constraints()
        self.add_intertwining_relation_constraints()
        if get_verbose() >= 2:
            self.show()

        self._solution_cache = []
        if solutions is not None:
            for solution in solutions:
                self._add_solution({(P.find(a), z): value
                                    for (a, z), value in solution.items()})

    def show(self, variables=True):
        r"""
        Print the constraints and variables of the MILP together
        with some explanations.

        EXAMPLES::

            sage: A = B = ["a", "b", "c"]
            sage: bij = Bijectionist(A, B, lambda x: A.index(x) % 2, solver='GLPK')
            sage: bij.set_constant_blocks([["a", "b"]])
            sage: next(bij.solutions_iterator())
            {'a': 0, 'b': 0, 'c': 1}
            sage: bij._bmilp.show()
            Constraints are:
                block a: 1 <= x_0 + x_1 <= 1
                block c: 1 <= x_2 + x_3 <= 1
                statistics: 2 <= 2 x_0 + x_2 <= 2
                statistics: 1 <= 2 x_1 + x_3 <= 1
                veto: x_0 + x_3 <= 1
            Variables are:
                x_0: s(a) = s(b) = 0
                x_1: s(a) = s(b) = 1
                x_2: s(c) = 0
                x_3: s(c) = 1
        """
        print("Constraints are:")
        b = self.milp.get_backend()
        varid_name = {}
        for i in range(b.ncols()):
            s = b.col_name(i)
            default_name = str(self.milp.linear_functions_parent()({i: 1}))
            if s and s != default_name:
                varid_name[i] = s
            else:
                varid_name[i] = default_name
        for i, (lb, (indices, values), ub) in enumerate(self.milp.constraints()):
            if b.row_name(i):
                print("    " + b.row_name(i) + ":", end=" ")
            if lb is not None:
                print(str(ZZ(lb)) + " <=", end=" ")
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
            print("<= " + str(ZZ(ub)) if ub is not None else "")

        if variables:
            print("Variables are:")
            P = self._bijectionist._P.root_to_elements_dict()
            for (p, z), v in self._x.items():
                print(f"    {v}: " + "".join([f"s({a}) = "
                                              for a in P[p]]) + f"{z}")

    def _prepare_solution(self, on_blocks, solution):
        r"""
        Return the solution as a dictionary from `A` (or `P`) to
        `Z`.

        INPUT:

        - ``on_blocks`` -- whether to return the solution on blocks
          or on all elements

        TESTS::

            sage: A = B = ["a", "b", "c"]
            sage: bij = Bijectionist(A, B, lambda x: 0)
            sage: bij.set_constant_blocks([["a", "b"]])
            sage: next(bij.solutions_iterator())
            {'a': 0, 'b': 0, 'c': 0}
            sage: bmilp = bij._bmilp
            sage: bmilp._prepare_solution(True, bmilp._solution_cache[0])
            {'a': 0, 'c': 0}
        """
        P = self._bijectionist._P
        tZ = self._bijectionist._possible_block_values
        mapping = {}  # A -> Z or P -> Z, a +-> s(a)
        for p, block in P.root_to_elements_dict().items():
            for z in tZ[p]:
                if solution[p, z] == 1:
                    if on_blocks:
                        mapping[p] = z
                    else:
                        for a in block:
                            mapping[a] = z
                    break
        return mapping

    def solutions_iterator(self, on_blocks, additional_constraints):
        r"""
        Iterate over all solutions satisfying the additional constraints.

        INPUT:

        - ``additional_constraints`` -- list of constraints for the
          underlying MILP

        - ``on_blocks`` -- whether to return the solution on blocks or
          on all elements

        TESTS::

            sage: A = B = 'abc'
            sage: bij = Bijectionist(A, B, lambda x: B.index(x) % 2, solver='GLPK')
            sage: from sage.combinat.bijectionist import _BijectionistMILP
            sage: bmilp = _BijectionistMILP(bij)
            sage: it = bmilp.solutions_iterator(False, [])
            sage: it2 = bmilp.solutions_iterator(False, [bmilp._x[('c', 1)] == 1])
            sage: next(it)
            {'a': 0, 'b': 1, 'c': 0}
            sage: next(it2)
            {'a': 0, 'b': 0, 'c': 1}
            sage: next(it)
            {'a': 0, 'b': 0, 'c': 1}
            sage: next(it)
            {'a': 1, 'b': 0, 'c': 0}
        """
        i = 0  # the first unconsidered element of _solution_cache
        while True:
            # skip solutions which do not satisfy additional_constraints
            while i < len(self._solution_cache):
                solution = self._solution_cache[i]
                i += 1
                if all(self._is_solution(constraint, solution)
                       for constraint in additional_constraints):
                    yield self._prepare_solution(on_blocks, solution)
                    break
            else:
                new_indices = []
                for constraint in additional_constraints:
                    new_indices.extend(self.milp.add_constraint(constraint,
                                                                return_indices=True))
                try:
                    self.milp.solve()
                    # moving this out of the try...finally block breaks SCIP
                    solution = self.milp.get_values(self._x,
                                                    convert=bool, tolerance=0.1)
                except MIPSolverException:
                    return
                finally:
                    b = self.milp.get_backend()
                    if hasattr(b, "_get_model"):
                        m = b._get_model()
                        if m.getStatus() != 'unknown':
                            m.freeTransform()
                    self.milp.remove_constraints(new_indices)

                self._add_solution(solution)
                i += 1
                assert i == len(self._solution_cache)
                yield self._prepare_solution(on_blocks, solution)
                if get_verbose() >= 2:
                    print("after vetoing")
                    self.show(variables=False)

    def _add_solution(self, solution):
        r"""
        Add the ``solution`` to the cache and an appropriate
        veto constraint to the MILP.

        INPUT:

        - ``solution`` -- dictionary from the indices of the MILP to
          a boolean

        EXAMPLES::

            sage: A = B = ["a", "b"]
            sage: bij = Bijectionist(A, B)
            sage: from sage.combinat.bijectionist import _BijectionistMILP
            sage: bmilp = _BijectionistMILP(bij)
            sage: bmilp._add_solution({(a, b): a == b for a in A for b in B})
            sage: bmilp.show()                                                  # random
            Constraints are:
                block a: 1 <= x_0 + x_1 <= 1
                block b: 1 <= x_2 + x_3 <= 1
                statistics: 1 <= x_1 + x_3 <= 1
                statistics: 1 <= x_0 + x_2 <= 1
                veto: x_1 + x_2 <= 1
            Variables are:
                x_0: s(a) = b
                x_1: s(a) = a
                x_2: s(b) = b
                x_3: s(b) = a
        """
        active_vars = [self._x[p, z]
                       for p in _disjoint_set_roots(self._bijectionist._P)
                       for z in self._bijectionist._possible_block_values[p]
                       if solution[(p, z)]]
        self.milp.add_constraint(sum(active_vars) <= len(active_vars) - 1,
                                 name='veto')
        self._solution_cache.append(solution)

    def _is_solution(self, constraint, values):
        r"""
        Evaluate the given function at the given values.

        INPUT:

        - ``constraint`` -- a
          :class:`sage.numerical.linear_functions.LinearConstraint`

        - ``values`` -- a candidate for a solution of the MILP as a
          dictionary from pairs `(a, z)\in A\times Z` to `0` or `1`,
          specifying whether `a` is mapped to `z`.

        EXAMPLES::

            sage: A = B = ["a", "b"]
            sage: bij = Bijectionist(A, B)
            sage: from sage.combinat.bijectionist import _BijectionistMILP
            sage: bmilp = _BijectionistMILP(bij)
            sage: f = bmilp._x["a", "a"] + bmilp._x["b", "a"] >= bmilp._x["b", "b"] + 1
            sage: v = {('a', 'a'): 1, ('a', 'b'): 0, ('b', 'a'): 1, ('b', 'b'): 1}
            sage: bmilp._is_solution(f, v)
            True
            sage: v = {('a', 'a'): 0, ('a', 'b'): 0, ('b', 'a'): 1, ('b', 'b'): 1}
            sage: bmilp._is_solution(f, v)
            False
        """
        index_block_value_dict = {}
        for (p, z), v in self._x.items():
            variable_index = next(iter(v.dict()))
            index_block_value_dict[variable_index] = (p, z)

        def evaluate(f):
            return sum(coeff if index == -1 else
                       coeff * values[index_block_value_dict[index]]
                       for index, coeff in f.dict().items())

        for lhs, rhs in constraint.equations():
            if evaluate(lhs - rhs):
                return False
        for lhs, rhs in constraint.inequalities():
            if evaluate(lhs - rhs) > 0:
                return False
        return True

    def add_alpha_beta_constraints(self):
        r"""
        Add constraints enforcing that `(alpha, s)` is equidistributed
        with `(beta, tau)` and `S` is the intertwining bijection.

        We do this by adding

        .. MATH::

            \sum_{a\in A, z\in Z} x_{p(a), z} s^z t^{\alpha(a)}
            = \sum_{b\in B} s^{\tau(b)} t(\beta(b))

        as a matrix equation.

        EXAMPLES::

            sage: A = B = [permutation for n in range(3) for permutation in Permutations(n)]
            sage: bij = Bijectionist(A, B, len)
            sage: bij.set_statistics((len, len))
            sage: from sage.combinat.bijectionist import _BijectionistMILP
            sage: bmilp = _BijectionistMILP(bij)                                # indirect doctest
            sage: next(bmilp.solutions_iterator(False, []))
            {[]: 0, [1]: 1, [1, 2]: 2, [2, 1]: 2}
        """
        W = self._bijectionist._W
        Z = self._bijectionist._Z
        zero = self.milp.linear_functions_parent().zero()
        AZ_matrix = [[zero] * len(W) for _ in range(len(Z))]
        B_matrix = [[zero] * len(W) for _ in range(len(Z))]

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

        for w in range(len(W)):
            for z in range(len(Z)):
                self.milp.add_constraint(AZ_matrix[z][w] == B_matrix[z][w],
                                         name='statistics')

    def add_distribution_constraints(self):
        r"""
        Add constraints so the distributions given by
        :meth:`set_distributions` are fulfilled.

        To accomplish this we add

        .. MATH::

            \sum_{a\in\tilde A} x_{p(a), z}t^z = \sum_{z\in\tilde Z} t^z,

        where `p(a)` is the block containing `a`, for each given
        distribution as a vector equation.

        EXAMPLES::

            sage: A = B = Permutations(3)
            sage: tau = Permutation.longest_increasing_subsequence_length
            sage: bij = Bijectionist(A, B, tau)
            sage: bij.set_distributions(([Permutation([1, 2, 3]), Permutation([1, 3, 2])], [1, 3]))
            sage: from sage.combinat.bijectionist import _BijectionistMILP
            sage: bmilp = _BijectionistMILP(bij)                                # indirect doctest
            sage: next(bmilp.solutions_iterator(False, []))
            {[1, 2, 3]: 3,
             [1, 3, 2]: 1,
             [2, 1, 3]: 2,
             [2, 3, 1]: 2,
             [3, 1, 2]: 2,
             [3, 2, 1]: 2}
        """
        Z = self._bijectionist._Z
        Z_dict = {z: i for i, z in enumerate(Z)}
        zero = self.milp.linear_functions_parent().zero()
        for tA, tZ in self._bijectionist._elements_distributions:
            tA_sum = [zero] * len(Z_dict)
            tZ_sum = [zero] * len(Z_dict)
            for a in tA:
                p = self._bijectionist._P.find(a)
                for z in self._bijectionist._possible_block_values[p]:
                    tA_sum[Z_dict[z]] += self._x[p, z]
            for z in tZ:
                tZ_sum[Z_dict[z]] += 1

            for a, z in zip(tA_sum, tZ_sum):
                self.milp.add_constraint(a == z, name=f"d: {a} == {z}")

    def add_intertwining_relation_constraints(self):
        r"""
        Add constraints corresponding to the given intertwining
        relations.

        INPUT:

        - origins, a list of triples `((\pi/\rho, p,
          (p_1,\dots,p_k))`, where `p` is the block of
          `\rho(s(a_1),\dots, s(a_k))`, for any `a_i\in p_i`.

        This adds the constraints imposed by
        :meth:`set_intertwining_relations`,

        .. MATH::

            s(\pi(a_1,\dots, a_k)) = \rho(s(a_1),\dots, s(a_k))

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

        Note that `z` must be a possible value of `p` and each `z_i`
        must be a possible value of `p_i`.

        EXAMPLES::

            sage: A = B = 'abcd'
            sage: bij = Bijectionist(A, B, lambda x: B.index(x) % 2)
            sage: def pi(p1, p2): return 'abcdefgh'[A.index(p1) + A.index(p2)]
            sage: def rho(s1, s2): return (s1 + s2) % 2
            sage: bij.set_intertwining_relations((2, pi, rho))
            sage: from sage.combinat.bijectionist import _BijectionistMILP
            sage: bmilp = _BijectionistMILP(bij)                                # indirect doctest
            sage: next(bmilp.solutions_iterator(False, []))
            {'a': 0, 'b': 1, 'c': 0, 'd': 1}
        """
        A = self._bijectionist._A
        tZ = self._bijectionist._possible_block_values
        P = self._bijectionist._P
        for composition_index, pi_rho in enumerate(self._bijectionist._pi_rho):
            pi_blocks = set()
            for a_tuple in itertools.product(A, repeat=pi_rho.numargs):
                if pi_rho.domain is not None and not pi_rho.domain(*a_tuple):
                    continue
                a = pi_rho.pi(*a_tuple)
                if a in A:
                    p_tuple = tuple(P.find(a) for a in a_tuple)
                    p = P.find(a)
                    if (p_tuple, p) not in pi_blocks:
                        pi_blocks.add((p_tuple, p))
                        for z_tuple in itertools.product(*[tZ[p] for p in p_tuple]):
                            rhs = (1 - pi_rho.numargs
                                   + sum(self._x[p_i, z_i]
                                         for p_i, z_i in zip(p_tuple, z_tuple)))
                            z = pi_rho.rho(*z_tuple)
                            if z in tZ[p]:
                                c = self._x[p, z] - rhs
                                if c.is_zero():
                                    continue
                                self.milp.add_constraint(c >= 0,
                                                         name=f"pi/rho({composition_index})")
                            else:
                                self.milp.add_constraint(rhs <= 0,
                                                         name=f"pi/rho({composition_index})")

    def add_quadratic_relation_constraints(self):
        r"""
        Add constraints enforcing that `s\circ\phi\circ s =
        \psi`.

        We do this by adding

        .. MATH::

            x_{p(a), z} = x_{p(\psi(z)), \phi(a)}

        for `a\in A` and `z\in Z`, where `\phi:A\to Z` and `\psi:Z\to
        A`.  Note that, in particular, `\phi` must be constant on
        blocks.

        EXAMPLES::

            sage: A = B = DyckWords(3)
            sage: bij = Bijectionist(A, B)
            sage: bij.set_statistics((lambda D: D.number_of_touch_points(), lambda D: D.number_of_initial_rises()))
            sage: ascii_art(sorted(bij.minimal_subdistributions_iterator()))
            [ (             [   /\   ] )
            [ (             [  /  \  ] )  ( [    /\    /\    ]  [  /\      /\/\  ] )
            [ ( [ /\/\/\ ], [ /    \ ] ), ( [ /\/  \, /  \/\ ], [ /  \/\, /    \ ] ),
            <BLANKLINE>
             ( [           /\   ]                     ) ]
             ( [  /\/\    /  \  ]  [            /\  ] ) ]
             ( [ /    \, /    \ ], [ /\/\/\, /\/  \ ] ) ]
            sage: bij.set_quadratic_relation((lambda D: D, lambda D: D))   # indirect doctest
            sage: ascii_art(sorted(bij.minimal_subdistributions_iterator()))
            [ (             [   /\   ] )
            [ (             [  /  \  ] )  ( [    /\  ]  [  /\/\  ] )
            [ ( [ /\/\/\ ], [ /    \ ] ), ( [ /\/  \ ], [ /    \ ] ),
            <BLANKLINE>
            <BLANKLINE>
             ( [  /\    ]  [  /\    ] )  ( [  /\/\  ]  [    /\  ] )
             ( [ /  \/\ ], [ /  \/\ ] ), ( [ /    \ ], [ /\/  \ ] ),
            <BLANKLINE>
             ( [   /\   ]             ) ]
             ( [  /  \  ]             ) ]
             ( [ /    \ ], [ /\/\/\ ] ) ]
        """
        P = self._bijectionist._P
        for phi, psi in self._bijectionist._phi_psi:
            for p, block in P.root_to_elements_dict().items():
                z0 = phi(p)
                assert all(phi(a) == z0 for a in block), "phi must be constant on the block %s" % block
                for z in self._bijectionist._possible_block_values[p]:
                    p0 = P.find(psi(z))
                    if z0 in self._bijectionist._possible_block_values[p0]:
                        c = self._x[p, z] - self._x[p0, z0]
                        if c.is_zero():
                            continue
                        self.milp.add_constraint(c == 0, name=f"i: s({p})={z}<->s(psi({z})=phi({p})")
                    else:
                        self.milp.add_constraint(self._x[p, z] == 0, name=f"i: s({p})!={z}")

    def add_homomesic_constraints(self):
        r"""
        Add constraints enforcing that `s` has constant average
        on the blocks of `Q`.

        We do this by adding

        .. MATH::

            \frac{1}{|q|}\sum_{a\in q} \sum_z z x_{p(a), z} =
            \frac{1}{|q_0|}\sum_{a\in q_0} \sum_z z x_{p(a), z},

        for `q\in Q`, where `q_0` is some fixed block of `Q`.

        EXAMPLES::

            sage: A = B = [1,2,3]
            sage: bij = Bijectionist(A, B, lambda b: b % 3)
            sage: bij.set_homomesic([[1,2], [3]])                               # indirect doctest
            sage: list(bij.solutions_iterator())
            [{1: 2, 2: 0, 3: 1}, {1: 0, 2: 2, 3: 1}]
        """
        Q = self._bijectionist._Q
        if Q is None:
            return
        P = self._bijectionist._P
        tZ = self._bijectionist._possible_block_values

        def sum_q(q):
            return sum(sum(z * self._x[P.find(a), z] for z in tZ[P.find(a)])
                       for a in q)
        q0 = Q[0]
        v0 = sum_q(q0)
        for q in Q[1:]:
            self.milp.add_constraint(len(q0) * sum_q(q) == len(q) * v0,
                                     name=f"h: ({q})~({q0})")


def _invert_dict(d):
    """
    Return the dictionary whose keys are the values of the input and
    whose values are the lists of preimages.

    INPUT:

    - ``d`` -- dictionary

    EXAMPLES::

        sage: from sage.combinat.bijectionist import _invert_dict
        sage: _invert_dict({1: "a", 2: "a", 3:"b"})
        {'a': [1, 2], 'b': [3]}

        sage: _invert_dict({})
        {}
    """
    preimages = {}
    for k, v in d.items():
        preimages[v] = preimages.get(v, []) + [k]
    return preimages


def _disjoint_set_roots(d):
    """
    Return the representatives of the blocks of the disjoint set.

    INPUT:

    - ``d`` -- a :class:`sage.sets.disjoint_set.DisjointSet_of_hashables`

    EXAMPLES::

        sage: from sage.combinat.bijectionist import _disjoint_set_roots
        sage: d = DisjointSet('abcde')
        sage: d.union("a", "b")
        sage: d.union("a", "c")
        sage: d.union("e", "d")
        sage: _disjoint_set_roots(d)
        dict_keys(['a', 'e'])
    """
    return d.root_to_elements_dict().keys()


def _non_copying_intersection(sets):
    """
    Return the intersection of the sets.

    If the intersection is equal to one of the sets, return this
    set.

    EXAMPLES::

        sage: from sage.combinat.bijectionist import _non_copying_intersection
        sage: A = set(range(7000)); B = set(range(8000));
        sage: _non_copying_intersection([A, B]) is A
        True

        sage: A = set([1,2]); B = set([2,3])
        sage: _non_copying_intersection([A, B])
        {2}
    """
    sets = sorted(sets, key=len)
    result = set.intersection(*sets)
    n = len(result)
    for s in sets:
        N = len(s)
        if n < N:
            return result
        if s == result:
            return s


"""
TESTS::

    sage: As = Bs = [[],
    ....:            [(1,i,j) for i in [-1,0,1] for j in [-1,1]],
    ....:            [(2,i,j) for i in [-1,0,1] for j in [-1,1]],
    ....:            [(3,i,j) for i in [-2,-1,0,1,2] for j in [-1,1]]]

Note that adding ``[(2,-2,-1), (2,2,-1), (2,-2,1), (2,2,1)]`` makes
it take (seemingly) forever::

    sage: def c1(a, b): return (a[0]+b[0], a[1]*b[1], a[2]*b[2])
    sage: def c2(a): return (a[0], -a[1], a[2])

    sage: bij = Bijectionist(sum(As, []), sum(Bs, []))
    sage: bij.set_statistics((lambda x: x[0], lambda x: x[0]))
    sage: bij.set_intertwining_relations((2, c1, c1), (1, c2, c2))
    sage: l = list(bij.solutions_iterator()); len(l)                            # long time -- (2.7 seconds with SCIP on AMD Ryzen 5 PRO 3500U w/ Radeon Vega Mobile Gfx)
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
    sage: l2 = [s for s in it if respects_c1(s) and respects_c2(s)]             # long time -- (17 seconds on AMD Ryzen 5 PRO 3500U w/ Radeon Vega Mobile Gfx)
    sage: sorted(l1, key=lambda s: tuple(s.items())) == l2                      # long time
    True

Our benchmark example::

    sage: from sage.combinat.cyclic_sieving_phenomenon import orbit_decomposition
    sage: def alpha1(p): return len(p.weak_excedences())
    sage: def alpha2(p): return len(p.fixed_points())
    sage: def beta1(p): return len(p.descents(final_descent=True)) if p else 0
    sage: def beta2(p): return len([e for (e, f) in zip(p, p[1:]+[0]) if e == f+1])
    sage: gamma = Permutation.longest_increasing_subsequence_length
    sage: def rotate_permutation(p):
    ....:    cycle = Permutation(tuple(range(1, len(p)+1)))
    ....:    return Permutation([cycle.inverse()(p(cycle(i))) for i in range(1, len(p)+1)])

    sage: N = 5
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
    ([[2, 3, 1, 5, 4], [2, 4, 5, 3, 1], [2, 5, 4, 1, 3], [3, 4, 1, 5, 2]], [2, 3, 3, 3])
    ([[3, 1, 2, 5, 4], [4, 1, 2, 5, 3], [3, 5, 2, 1, 4], [4, 1, 5, 2, 3]], [3, 3, 4, 4])
    ([[2, 1, 3, 5, 4], [2, 4, 1, 3, 5], [2, 5, 3, 1, 4], [3, 4, 1, 2, 5], [3, 1, 5, 4, 2], [2, 5, 1, 4, 3], [2, 1, 5, 4, 3]], [2, 2, 3, 3, 3, 3, 3])

    sage: l = list(bij.solutions_iterator()); len(l)                            # not tested -- (17 seconds with SCIP on AMD Ryzen 5 PRO 3500U w/ Radeon Vega Mobile Gfx)
    504

    sage: for a, d in bij.minimal_subdistributions_iterator():                  # not tested
    ....:     print(sorted(a), sorted(d))
"""
