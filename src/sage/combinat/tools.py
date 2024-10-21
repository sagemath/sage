# sage_setup: distribution = sagemath-categories
r"""
Transitive ideal closure tool
"""
# ****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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
# ****************************************************************************


def transitive_ideal(f, x):
    r"""
    Return a list of all elements reachable from `x` in the abstract
    reduction system whose reduction relation is given by the function
    `f`.

    In more elementary terms:

    If `S` is a set, and `f` is a function sending every element of `S`
    to a list of elements of `S`, then we can define a digraph on the
    vertex set `S` by drawing an edge from `s` to `t` for every
    `s \in S` and every `t \in f(s)`.

    If `x \in S`, then an element `y \in S` is said to be reachable
    from `x` if there is a path `x \to y` in this graph.

    Given `f` and `x`, this method computes the list of all elements of
    `S` reachable from `x`.

    Note that if there are infinitely many such elements, then this
    method will never halt.

    For more powerful versions, see :mod:`sage.combinat.backtrack`

    EXAMPLES::

        sage: f = lambda x: [x-1] if x > 0 else []
        sage: sage.combinat.tools.transitive_ideal(f, 10)
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    """
    known = [x]
    todo = [x]

    while todo:
        y = todo.pop(0)
        for z in f(y):
            if z in known:
                continue
            known.append(z)
            todo.append(z)

    known.sort()
    return known
