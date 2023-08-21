r"""
Pseudolines

This module gathers everything that has to do with pseudolines, and for a start
a :class:`PseudolineArrangement` class that can be used to describe an
arrangement of pseudolines in several different ways, and to translate one
description into another, as well as to display *Wiring diagrams* via the
:meth:`show <sage.geometry.pseudolines.PseudolineArrangement.show>` method.

In the following, we try to stick to the terminology given in [Fe1997]_, which
can be checked in case of doubt. And please fix this module's documentation
afterwards :-)

**Definition**

A *pseudoline* can not be defined by itself, though it can be thought of as a
`x`-monotone curve in the plane. A *set* of pseudolines, however, represents a
set of such curves that pairwise intersect exactly once (and hence mimic the
behaviour of straight lines in general position). We also assume that those
pseudolines are in general position, that is that no three of them cross at the
same point.

The present class is made to deal with a combinatorial encoding of a pseudolines
arrangement, that is the ordering in which a pseudoline `l_i` of an arrangement
`l_0, ..., l_{n-1}` crosses the `n-1` other lines.

.. WARNING::

    It is assumed through all the methods that the given lines are numbered
    according to their `y`-coordinate on the vertical line `x=-\infty`.
    For instance, it is not possible that the first transposition be ``(0,2)``
    (or equivalently that the first line `l_0` crosses is `l_2` and conversely),
    because one of them would have to cross `l_1` first.

Encodings
----------

**Permutations**

An arrangement of pseudolines can be described by a sequence of `n` lists of
length `n-1`, where the `i` list is a permutation of `\{0, ..., n-1\} \backslash
i` representing the ordering in which the `i` th pseudoline meets the other
ones.

::

    sage: from sage.geometry.pseudolines import PseudolineArrangement
    sage: permutations = [[3, 2, 1], [3, 2, 0], [3, 1, 0], [2, 1, 0]]
    sage: p = PseudolineArrangement(permutations)
    sage: p
    Arrangement of pseudolines of size 4
    sage: p.show()                                                                      # needs sage.plot

**Sequence of transpositions**

An arrangement of pseudolines can also be described as a sequence of `\binom n
2` transpositions (permutations of two elements). In this sequence, the
transposition `(2,3)` appears before `(8, 2)` iif `l_2` crosses `l_3` before it
crosses `l_8`. This encoding is easy to obtain by reading the wiring diagram
from left to right (see the :meth:`show
<sage.geometry.pseudolines.PseudolineArrangement.show>` method).

::

    sage: from sage.geometry.pseudolines import PseudolineArrangement
    sage: transpositions = [(3, 2), (3, 1), (0, 3), (2, 1), (0, 2), (0, 1)]
    sage: p = PseudolineArrangement(transpositions)
    sage: p
    Arrangement of pseudolines of size 4
    sage: p.show()                                                                      # needs sage.plot


Note that this ordering is not necessarily unique.

**Felsner's Matrix**

Felser gave an encoding of an arrangement of pseudolines that takes `n^2` bits
instead of the `n^2\log(n)` bits required by the two previous encodings.

Instead of storing the permutation ``[3, 2, 1]`` to remember that line `l_0`
crosses `l_3` then `l_2` then `l_1`, it is sufficient to remember the positions
for which each line `l_i` meets a line `l_j` with `j < i`. As `l_0` -- the first
of the lines -- can only meet pseudolines with higher index, we can store ``[0,
0, 0]`` instead of ``[3, 2, 1]`` stored previously. For `l_1`'s permutation
``[3, 2, 0]`` we only need to remember that `l_1` first crosses 2 pseudolines of
higher index, and then a pseudoline with smaller index, which yields the bit
vector ``[0, 0, 1]``. Hence we can transform the list of permutations above into
a list of `n` bit vectors of length `n-1`, that is

.. MATH::

    \begin{array}{ccc}
      3 & 2 & 1\\
      3 & 2 & 0\\
      3 & 1 & 0\\
      2 & 1 & 0\\
    \end{array}
    \Rightarrow
    \begin{array}{ccc}
      0 & 0 & 0\\
      0 & 0 & 1\\
      0 & 1 & 1\\
      1 & 1 & 1\\
    \end{array}

In order to go back from Felsner's matrix to an encoding by a sequence of
transpositions, it is sufficient to look for occurrences of
`\begin{array}{c}0\\1\end{array}` in the first column of the matrix, as it
corresponds in the wiring diagram to a line going up while the line immediately
above it goes down -- those two lines cross. Each time such a pattern is found
it yields a new transposition, and the matrix can be updated so that this
pattern disappears. A more detailed description of this algorithm is given in
[Fe1997]_.

::

    sage: from sage.geometry.pseudolines import PseudolineArrangement
    sage: felsner_matrix = [[0, 0, 0], [0, 0, 1], [0, 1, 1], [1, 1, 1]]
    sage: p = PseudolineArrangement(felsner_matrix)
    sage: p
    Arrangement of pseudolines of size 4

Example
-------

Let us define in the plane several lines `l_i` of equation `y = a x+b` by
picking a coefficient `a` and `b` for each of them. We make sure that no two of
them are parallel by making sure all of the `a` chosen are different, and we
avoid a common crossing of three lines by adding a random noise to `b`::

    sage: n = 20
    sage: l = sorted(zip(Subsets(20*n, n).random_element(),
    ....:                [randint(0, 20*n) + random() for i in range(n)]))
    sage: print(l[:5])                          # not tested                            # needs sage.combinat
    [(96, 278.0130613051349), (74, 332.92512282478714), (13, 155.65820951249867),
     (209, 34.753946221755307), (147, 193.51376457741441)]

We can now compute for each `i` the order in which line `i` meets the other lines::

    sage: permutations = [[0..i-1] + [i+1..n-1] for i in range(n)]
    sage: def a(x): return l[x][0]
    sage: def b(x): return l[x][1]
    sage: for i, perm in enumerate(permutations):
    ....:     perm.sort(key=lambda j: (b(j)-b(i))/(a(i)-a(j)))

And finally build the line arrangement::

    sage: from sage.geometry.pseudolines import PseudolineArrangement
    sage: p = PseudolineArrangement(permutations)
    sage: print(p)
    Arrangement of pseudolines of size 20
    sage: p.show(figsize=[20,8])                                                        # needs sage.combinat sage.plot

Author
^^^^^^
Nathann Cohen

Methods
-------
"""
##############################################################################
#       Copyright (C) 2011 Nathann Cohen <nathann.cohen@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  https://www.gnu.org/licenses/
##############################################################################

from copy import deepcopy


class PseudolineArrangement:

    def __init__(self, seq, encoding="auto"):
        r"""
        Creates an arrangement of pseudolines.

        INPUT:

        - ``seq`` (a sequence describing the line arrangement). It can be:

          - A list of `n` permutations of size `n-1`.
          - A list of `\binom n 2` transpositions
          - A Felsner matrix, given as a sequence of `n` binary vectors of
            length `n-1`.

        - ``encoding`` (information on how the data should be interpreted), and
          can assume any value among 'transpositions', 'permutations', 'Felsner'
          or 'auto'. In the latter case, the type will be guessed (default
          behaviour).

        .. NOTE::

           * The pseudolines are assumed to be integers `0,\dots,n-1`.

           * For more information on the different encodings, see the
             :mod:`pseudolines module <sage.geometry.pseudolines>`'s
             documentation.

        TESTS:

        From permutations::

            sage: from sage.geometry.pseudolines import PseudolineArrangement
            sage: permutations = [[3, 2, 1], [3, 2, 0], [3, 1, 0], [2, 1, 0]]
            sage: PseudolineArrangement(permutations)
            Arrangement of pseudolines of size 4

        From transpositions ::

            sage: from sage.geometry.pseudolines import PseudolineArrangement
            sage: transpositions = [(3, 2), (3, 1), (0, 3), (2, 1), (0, 2), (0, 1)]
            sage: PseudolineArrangement(transpositions)
            Arrangement of pseudolines of size 4

        From a Felsner matrix::

            sage: from sage.geometry.pseudolines import PseudolineArrangement
            sage: permutations = [[3, 2, 1], [3, 2, 0], [3, 1, 0], [2, 1, 0]]
            sage: p = PseudolineArrangement(permutations)
            sage: matrix = p.felsner_matrix()
            sage: PseudolineArrangement(matrix) == p
            True

        Wrong input::

            sage: PseudolineArrangement([[5, 2, 1], [3, 2, 0], [3, 1, 0], [2, 1, 0]])
            Traceback (most recent call last):
            ...
            ValueError: Are the lines really numbered from 0 to n-1?
            sage: PseudolineArrangement([(3, 2), (3, 1), (0, 3), (2, 1), (0, 2)])
            Traceback (most recent call last):
            ...
            ValueError: A line is numbered 3 but the number of transpositions ...
        """

        # Sequence of transpositions
        if (encoding == "transpositions" or
            (encoding == "auto" and len(seq[0]) == 2 and len(seq) > 3)):

            self._n = max(map(max, seq)) + 1
            if (self._n * (self._n-1))/2 != len(seq):
                raise ValueError(
                    "A line is numbered "+str(self._n-1)+" but the number"+
                    " of transpositions is different from binomial("+
                    str(self._n-1)+",2). Are the lines numbered from 0 to n-1?"+
                    " Are they really non-parallel? Please check the documentation.")

            self._permutations = [[] for i in range(self._n)]

            for i,j in seq:
                self._permutations[i].append(j)
                self._permutations[j].append(i)

        # Sequence of permutations
        elif (encoding == "permutations" or
            (encoding == "auto" and (len(seq[0]) == len(seq)-1) and max(seq[0]) > 1)):

            self._n = len(seq)
            self._permutations = [list(_) for _ in seq]

            if max(map(max, seq)) != self._n -1 :
                raise ValueError("Are the lines really numbered from 0 to n-1?")

        # Felsner encoding
        elif (encoding == "Felsner" or
            (encoding == "auto" and len(seq[0]) == len(seq) -1)):

            seq = deepcopy(seq)
            self._n = len(seq)
            ordering = list(range(self._n))

            self._permutations = [[] for i in range(self._n)]

            crossings = (self._n * (self._n-1))/2

            i = 0
            while crossings > 0:
                if (seq[i] != [] and
                    (seq[i][0] == 0 and
                     seq[i+1][0] == 1)):

                    crossings -= 1

                    self._permutations[ordering[i]].append(ordering[i+1])
                    self._permutations[ordering[i+1]].append(ordering[i])

                    ordering[i], ordering[i+1] = ordering[i+1], ordering[i]
                    seq[i], seq[i+1] = seq[i+1], seq[i]

                    seq[i].pop(0)
                    seq[i+1].pop(0)

                    if i > 0 and seq[i-1] is not []:
                        i -= 1
                    else:
                        i += 1
                else:
                    i += 1
        else:

            if encoding != "auto":
                raise ValueError("The value of encoding must be one of 'transpositions', 'permutations', 'Felsner' or 'auto'.")

            raise ValueError("The encoding you used could not be guessed. Your input string is probably badly formatted, or you have at most 3 lines and we cannot distinguish the encoding. Please specify the encoding you used.")

    def transpositions(self):
        r"""
        Return the arrangement as `\binom n 2` transpositions.

        See the :mod:`pseudolines module <sage.geometry.pseudolines>`'s
        documentation for more information on this encoding.

        EXAMPLES::

            sage: from sage.geometry.pseudolines import PseudolineArrangement
            sage: permutations = [[3, 2, 1], [3, 2, 0], [3, 1, 0], [2, 1, 0]]
            sage: p1 = PseudolineArrangement(permutations)
            sage: transpositions = [(3, 2), (3, 1), (0, 3), (2, 1), (0, 2), (0, 1)]
            sage: p2 = PseudolineArrangement(transpositions)
            sage: p1 == p2
            True
            sage: p1.transpositions()
            [(3, 2), (3, 1), (0, 3), (2, 1), (0, 2), (0, 1)]
            sage: p2.transpositions()
            [(3, 2), (3, 1), (0, 3), (2, 1), (0, 2), (0, 1)]
        """
        t = []
        perm = deepcopy(self._permutations)

        crossings = (self._n * (self._n-1))/2

        while crossings > 0:

            i = 0

            while perm[i] == []:
                i += 1

            k = 0
            while i != perm[perm[i][0]][0]:
                i = perm[i][0]
                k+= 1

                if k > self._n:
                    raise ValueError(
                        "It looks like the data does not correspond to a"+
                        "pseudoline arrangement. We have found k>2 lines"+
                        "such that the ith line meets the (i+1)th before"+
                        " the (i-1)th (this creates a cyclic dependency)"+
                        " which is totally impossible.")

            t.append((i, perm[i][0]))
            perm[perm[i][0]].pop(0)
            perm[i].pop(0)

            crossings -= 1

        if max(map(len, perm)) != 0:
            raise ValueError("There has been an error while computing the transpositions.")

        return t

    def permutations(self):
        r"""
        Return the arrangements as `n` permutations of size `n-1`.

        See the :mod:`pseudolines module <sage.geometry.pseudolines>`'s
        documentation for more information on this encoding.

        EXAMPLES::

            sage: from sage.geometry.pseudolines import PseudolineArrangement
            sage: permutations = [[3, 2, 1], [3, 2, 0], [3, 1, 0], [2, 1, 0]]
            sage: p = PseudolineArrangement(permutations)
            sage: p.permutations()
            [[3, 2, 1], [3, 2, 0], [3, 1, 0], [2, 1, 0]]
        """
        return deepcopy(self._permutations)

    def felsner_matrix(self):
        r"""
        Return a Felsner matrix describing the arrangement.

        See the :mod:`pseudolines module <sage.geometry.pseudolines>`'s
        documentation for more information on this encoding.

        EXAMPLES::

            sage: from sage.geometry.pseudolines import PseudolineArrangement
            sage: permutations = [[3, 2, 1], [3, 2, 0], [3, 1, 0], [2, 1, 0]]
            sage: p = PseudolineArrangement(permutations)
            sage: p.felsner_matrix()
            [[0, 0, 0], [0, 0, 1], [0, 1, 1], [1, 1, 1]]
        """

        m = [[] for i in range(self._n)]

        for i,j in self.transpositions():
            if i < j:
                i, j = j, i

            m[j].append(0)
            m[i].append(1)

        return m

    def show(self, **args):
        r"""
        Displays the pseudoline arrangement as a wiring diagram.

        INPUT:

        - ``**args`` -- any arguments to be forwarded to the ``show`` method. In
          particular, to tune the dimensions, use the ``figsize`` argument
          (example below).

        EXAMPLES::

            sage: from sage.geometry.pseudolines import PseudolineArrangement
            sage: permutations = [[3, 2, 1], [3, 2, 0], [3, 1, 0], [2, 1, 0]]
            sage: p = PseudolineArrangement(permutations)
            sage: p.show(figsize=[7,5])                                                 # needs sage.plot

        TESTS::

            sage: from sage.geometry.pseudolines import PseudolineArrangement
            sage: permutations = [[3, 2, 1], [3, 2, 0], [3, 0, 1], [2, 0, 1]]
            sage: p = PseudolineArrangement(permutations)
            sage: p.show()                                                              # needs sage.plot
            Traceback (most recent call last):
            ...
            ValueError: There has been a problem while plotting the figure...
        """
        x = 1
        from sage.plot.line import line
        from sage.plot.text import text

        lines = [[(0,self._n-1-i)] for i in range(self._n)]

        for i,j in self.transpositions():
            iy = lines[i][-1][1]
            jy = lines[j][-1][1]

            lines[i].append((x, iy))
            lines[j].append((x, jy))

            if abs(iy-jy) != 1:
                raise ValueError(
                    "There has been a problem while plotting the figure. It "+
                    "seems that the lines are not correctly ordered. Please "+
                    "check the pseudolines modules documentation, there is a "
                    +"warning about that. ")

            lines[i].append((x+2,jy))
            lines[j].append((x+2,iy))

            x += 2

        L = line([(1,1)])

        for i, l in enumerate(lines):
            l.append((x+2, l[-1][1]))
            L += line(l)

            L += text(str(i), (0, l[0][1]+.3), horizontal_alignment="right")
            L += text(str(i), (x+2, l[-1][1]+.3), horizontal_alignment="left")

        return L.show(axes=False, **args)

    def __repr__(self):
        r"""
        A short txt description of the pseudoline arrangement.

        EXAMPLES::

            sage: from sage.geometry.pseudolines import PseudolineArrangement
            sage: permutations = [[3, 2, 1], [3, 2, 0], [3, 1, 0], [2, 1, 0]]
            sage: p = PseudolineArrangement(permutations)
            sage: p
            Arrangement of pseudolines of size 4
        """
        return "Arrangement of pseudolines of size " + str(self._n)

    def __eq__(self, other):
        r"""
        Test of equality.

        TESTS::

            sage: from sage.geometry.pseudolines import PseudolineArrangement
            sage: permutations = [[3, 2, 1], [3, 2, 0], [3, 1, 0], [2, 1, 0]]
            sage: p1 = PseudolineArrangement(permutations)
            sage: transpositions = [(3, 2), (3, 1), (0, 3), (2, 1), (0, 2), (0, 1)]
            sage: p2 = PseudolineArrangement(transpositions)
            sage: p1 == p2
            True
        """
        return (self._n == other._n) and (self._permutations == other._permutations)

    def __ne__(self, other):
        """
        Test for non-equality.

        TESTS::

            sage: from sage.geometry.pseudolines import PseudolineArrangement
            sage: permutations = [[3, 2, 1], [3, 2, 0], [3, 1, 0], [2, 1, 0]]
            sage: p1 = PseudolineArrangement(permutations)
            sage: transpositions = [(3, 2), (3, 1), (0, 3), (2, 1), (0, 2), (0, 1)]
            sage: p2 = PseudolineArrangement(transpositions)
            sage: p1 != p2
            False
        """
        return not (self == other)
