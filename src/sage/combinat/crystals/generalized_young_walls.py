# sage.doctest: needs sage.combinat sage.modules
r"""
Crystals of Generalized Young Walls

AUTHORS:

- Lucas David-Roesler: Initial version

- Ben Salisbury: Initial version

- Travis Scrimshaw: Initial version

Generalized Young walls are certain generalizations of Young tableaux
introduced in [KS2010]_ and designed to be a realization of the crystals
`\mathcal{B}(\infty)` and `\mathcal{B}(\lambda)` in type `A_n^{(1)}`.

REFERENCES:

- [KLRS2016]_
- [KS2010]_
"""

# *****************************************************************************
#  Copyright (C) 2013
#
#  Lucas David-Roesler (roesler at lvc dot edu)
#  Ben Salisbury (bsalisbury at ccny dot cuny dot edu)
#  Travis Scrimshaw (tscrim at ucdavis dot edu)
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# *****************************************************************************

import re
from copy import deepcopy
from sage.combinat.root_system.cartan_type import CartanType
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.combinat.combinat import CombinatorialElement
from sage.categories.regular_crystals import RegularCrystals
from sage.categories.highest_weight_crystals import HighestWeightCrystals
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.combinat.root_system.root_system import RootSystem


class GeneralizedYoungWall(CombinatorialElement):
    r"""
    A generalized Young wall.

    For more information, see
    :class:`~sage.combinat.crystals.generalized_young_walls.InfinityCrystalOfGeneralizedYoungWalls`.

    EXAMPLES::

        sage: Y = crystals.infinity.GeneralizedYoungWalls(4)
        sage: mg = Y.module_generators[0]; mg.pp()
        0
        sage: mg.f_string([1,2,0,1]).pp()
        1|2|
        0|1|
           |
    """

    def __init__(self, parent, data):
        r"""
        EXAMPLES::

            sage: Y = crystals.infinity.GeneralizedYoungWalls(2)
            sage: mg = Y.module_generators[0]
            sage: TestSuite(mg).run()
        """
        i = len(data) - 1
        while i >= 0 and not data[i]:
            data.pop()
            i -= 1
        self.rows = len(data)
        if not data:
            self.cols = 0
        else:
            self.cols = max(len(r) for r in data)
        self.data = data
        CombinatorialElement.__init__(self, parent, data)

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: y = crystals.infinity.GeneralizedYoungWalls(3)([[0],[1,0,3,2],[2,1],[3,2,1,0,3,2],[0],[],[2]])
            sage: y
            [[0], [1, 0, 3, 2], [2, 1], [3, 2, 1, 0, 3, 2], [0], [], [2]]
        """
        return repr(self.data)

    def _repr_diagram(self):
        r"""
        Return a string representation of the diagram of ``self``.

        EXAMPLES::

            sage: y = crystals.infinity.GeneralizedYoungWalls(2)([[0,2,1],[1,0,2,1,0],[],[0],[1,0,2],[],[],[1]])
            sage: print(y._repr_diagram())
                    1|
                     |
                     |
                2|0|1|
                    0|
                     |
            0|1|2|0|1|
                1|2|0|
        """
        if not self.data:
            return '0'
        ret = ""
        for row in reversed(self.data):
            wall = ''
            for elem in reversed(row):
                wall += str(elem)
                wall += '|'
            if row == []:
                wall += '|'
            ret += wall.rjust(2 * self.cols + 1) + "\n"
        return ret

    def _ascii_art_(self):
        r"""
        Return an ascii art representation of ``self``.

        EXAMPLES::

            sage: y = crystals.infinity.GeneralizedYoungWalls(2)([[0,2,1],[1,0,2,1,0],[],[0],[1,0,2],[],[],[1]])
            sage: ascii_art(y)
                    1|
                     |
                     |
                2|0|1|
                    0|
                     |
            0|1|2|0|1|
                1|2|0|
        """
        from sage.typeset.ascii_art import AsciiArt
        return AsciiArt(self._repr_diagram().splitlines())

    def _unicode_art_(self):
        """
        Return a unicode art representation of ``self``.

        TESTS::

            sage: y = crystals.infinity.GeneralizedYoungWalls(2)([[0,2,1],[1,0,2,1,0],[],[0],[1,0,2],[],[],[1]])
            sage: unicode_art(y)
                            ┌───┐
                            │ 1 │
                            └───┘
                                │
                                ┤
                                │
                    ┌───┬───┬───┐
                    │ 1 │ 0 │ 2 │
                    └───┴───┼───┤
                            │ 0 │
                            └───┘
                                │
            ┌───┬───┬───┬───┬───┐
            │ 1 │ 0 │ 2 │ 1 │ 0 │
            └───┴───┼───┼───┼───┤
                    │ 0 │ 2 │ 1 │
                    └───┴───┴───┘
        """
        from sage.typeset.unicode_art import UnicodeArt
        if not self.data:
            return UnicodeArt(["0"])

        from sage.combinat.output import ascii_art_table
        import unicodedata
        v = unicodedata.lookup('BOX DRAWINGS LIGHT VERTICAL')
        vl = unicodedata.lookup('BOX DRAWINGS LIGHT VERTICAL AND LEFT')
        table = [[None] * (self.cols - len(row)) + row for row in reversed(self)]
        ret = []
        for i, row in enumerate(ascii_art_table(table, use_unicode=True).splitlines()):
            if row[-1] == " ":
                if i % 2 == 0:
                    ret.append(row[:-1] + vl)
                else:
                    ret.append(row[:-1] + v)
            else:
                ret.append(row)
        return UnicodeArt(ret)

    def __eq__(self, other):
        r"""
        EXAMPLES::

            sage: GYW = crystals.infinity.GeneralizedYoungWalls(2)
            sage: y = GYW([[],[1,0],[2,1]])
            sage: x = GYW([[],[1,0],[2,1]])
            sage: z = GYW([[],[1],[2]])
            sage: x == y
            True
            sage: x == z
            False
        """
        if isinstance(other, GeneralizedYoungWall):
            return self.data == other.data
        return self.data == other

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: GYW = crystals.infinity.GeneralizedYoungWalls(2)
            sage: h = hash(GYW)
        """
        return hash(tuple(tuple(u) for u in self.data))

    def raw_signature(self, i):
        r"""
        Return the sequence from `\{+,-\}` obtained from all `i`-admissible
        slots and removable `i`-boxes  without canceling any `(+,-)`-pairs.
        The result also notes the row and column of the sign.

        EXAMPLES::

            sage: x = crystals.infinity.GeneralizedYoungWalls(3)([[],[1,0,3,2],[2,1],[3,2,1,0,3,2],[],[],[2]])
            sage: x.raw_signature(2)
            [['-', 3, 6], ['-', 1, 4], ['-', 6, 1]]
        """
        sig = []
        rank = self.parent().cartan_type().rank()  # n+1
        for row in range(self.rows):
            if self.data[row] == [] and i == (row % rank):
                sig.append(['+', row, 0])
            elif self.data[row] == []:
                continue
            elif self.data[row][-1] == ((i + 1) % rank):
                sig.append(['+', row, len(self.data[row]) + 1])
            elif self.data[row][-1] == i:
                sig.append(['-', row, len(self.data[row])])
        return sorted(sig, key=self._sig_sort)

    def _sig_sort(self, a):
        r"""
        Internal command used to appropriately sort the output
        from :meth:`raw_signature()`.

        INPUT:

        - ``a`` -- list of the form ``['s',j,k]`` where `s` is a string, `j` is
          an integer and `k` is an integer

        EXAMPLES::

            sage: hw = crystals.infinity.GeneralizedYoungWalls(5)([])
            sage: hw._sig_sort(['+',1,0])
            (0, 1)
        """
        return (-a[2], a[1])

    def generate_signature(self, i):
        r"""
        The `i`-signature of ``self`` (with whitespace where cancellation
        occurs) together with the unreduced sequence from `\{+,-\}`.  The
        result also records to the row and column position of the sign.

        EXAMPLES::

            sage: y = crystals.infinity.GeneralizedYoungWalls(2)([[0],[1,0],[2,1,0,2],[],[1]])
            sage: y.generate_signature(1)
            ([['+', 2, 5], ['-', 4, 1]], '  ')
        """
        sig = []
        rank = self.parent().cartan_type().classical().rank()
        for row in range(self.rows):
            if self.data[row] == [] and i == (row % (rank + 1)):
                sig.append(['+', row, 0])
            elif self.data[row] == []:
                continue
            elif self.data[row][-1] == ((i + 1) % (rank + 1)):
                sig.append(['+', row, len(self.data[row]) + 1])
            elif self.data[row][-1] == i:
                sig.append(['-', row, len(self.data[row])])
        sig = sorted(sig, key=self._sig_sort)
        strsig = ''.join(x[0] for x in sig)
        reducedsig = strsig
        while re.search(r"\+\s*-", reducedsig):
            reducedsig = re.sub(r"\+\s*-", lambda match: ''.ljust(len(match.group(0))), reducedsig)
        return (sig, reducedsig)

    def signature(self, i):
        r"""
        Return the `i`-signature of ``self``.

        The signature is obtained by reading ``self`` in columns bottom to top starting from the left.
        Then add a `-` at every `i`-box which may be removed from ``self`` and still obtain a legal
        generalized Young wall, and add a `+` at each site for which an `i`-box may be added and still
        obtain a valid generalized Young wall.  Then successively cancel any `(+,-)`-pair to obtain a
        sequence of the form `- \cdots -+ \cdots +`.  This resulting sequence is the output.

        EXAMPLES::

            sage: y = crystals.infinity.GeneralizedYoungWalls(2)([[0],[1,0],[2,1,0,2],[],[1]])
            sage: y.signature(1)
            ''

            sage: x = crystals.infinity.GeneralizedYoungWalls(3)([[],[1,0,3,2],[2,1],[3,2,1,0,3,2],[],[],[2]])
            sage: x.signature(2)
            '---'
        """
        return self.generate_signature(i)[1].strip()

    def pp(self):
        r"""
        Pretty print ``self``.

        EXAMPLES::

            sage: y = crystals.infinity.GeneralizedYoungWalls(2)([[0,2,1],[1,0,2,1,0],[],[0],[1,0,2],[],[],[1]])
            sage: y.pp()
                    1|
                     |
                     |
                2|0|1|
                    0|
                     |
            0|1|2|0|1|
                1|2|0|
        """
        print(self._repr_diagram())

    def content(self):
        r"""
        Return total number of blocks in ``self``.

        EXAMPLES::

            sage: y = crystals.infinity.GeneralizedYoungWalls(2)([[0],[1,0],[2,1,0,2],[],[1]])
            sage: y.content()
            8

            sage: x = crystals.infinity.GeneralizedYoungWalls(3)([[],[1,0,3,2],[2,1],[3,2,1,0,3,2],[],[],[2]])
            sage: x.content()
            13
        """
        return sum(len(r) for r in self.data)

    def number_of_parts(self):
        r"""
        Return the value of `\mathscr{N}` on ``self``.

        In [KLRS2016]_, the statistic `\mathscr{N}` was defined on elements in
        `\mathcal{Y}(\infty)` which counts how many parts are in the
        corresponding Kostant partition.  Specifically, the computation of
        `\mathscr{N}(Y)` is done using the following algorithm:

        - If `Y` has no rows whose right-most box is colored `n` and such that
          the length of this row is a multiple of `n+1`, then `\mathscr{N}(Y)`
          is the total number of distinct rows in `Y`, not counting multiplicity.

        - Otherwise, search `Y` for the longest row such that the right-most box
          is colored `n` and such that the total number of boxes in the row is
          `k(n+1)` for some `k\ge 1`.  Replace this row by `n+1` distinct rows
          of length `k`, reordering all rows, if necessary, so that the result
          is a proper wall.  (Note that the resulting wall may no longer be
          reduced.) Repeat the search and replace process for all other rows of
          the above form for each `k' < k`.  Then `\mathscr{N}(Y)` is the number
          of distinct rows, not counting multiplicity, in the wall resulting
          from this process.

        EXAMPLES::

            sage: Y = crystals.infinity.GeneralizedYoungWalls(3)
            sage: y = Y([[0],[],[],[],[0],[],[],[],[0]])
            sage: y.number_of_parts()
            1

            sage: Y = crystals.infinity.GeneralizedYoungWalls(3)
            sage: y = Y([[0,3,2],[1,0],[],[],[0,3],[1,0],[],[],[0]])
            sage: y.number_of_parts()
            4

            sage: Y = crystals.infinity.GeneralizedYoungWalls(2)
            sage: y = Y([[0,2,1],[1,0],[2,1,0,2,1,0,2,1,0],[],[2,1,0,2,1,0]])
            sage: y.number_of_parts()
            8
        """
        n = self.parent().cartan_type().rank() - 1
        new = self.data[:]
        i = 0
        while i < len(new):
            r = new[i]
            if not r or r in new[i + 1:]:
                new.pop(i)
            elif r[0] == n and not len(r) % (n + 1):
                for j in range(n + 1):
                    temp = [k % (n + 1)
                            for k in range(j + len(r) // (n + 1) - 1, j - 1, -1)]
                    if temp not in new:
                        new.insert(i + 1, temp)
                new.pop(i)
            else:
                i += 1
        return len(new)

    def sum_of_weighted_row_lengths(self):
        r"""
        Return the value of `\mathscr{M}` on ``self``.

        Let `\mathcal{Y}_0 \subset \mathcal{Y}(\infty)` be the set of
        generalized Young walls which have no rows whose right-most box is
        colored `n`.  For `Y \in \mathcal{Y}_0`,

        .. MATH::

            \mathscr{M}(Y) = \sum_{i=1}^n (i+1)M_i(Y),

        where `M_i(Y)` is the number of nonempty rows in `Y` whose right-most
        box is colored `i-1`.

        EXAMPLES::

            sage: Y = crystals.infinity.GeneralizedYoungWalls(2)
            sage: y = Y([[0,2,1,0,2],[1,0,2],[],[0,2],[1,0],[],[0],[1,0]])
            sage: y.sum_of_weighted_row_lengths()
            15
        """
        n = self.parent().cartan_type().rank() - 1
        m = lambda i: len([1 for r in self.data if r and r[0] == (i - 1) % (n + 1)])
        for r in self.data:
            if r and r[0] == n:
                raise ValueError('Statistic only valid for generalized Young walls in Y_0')
        return sum((i + 1) * m(i) for i in range(1, n + 1))

    def e(self, i):
        r"""
        Return the application of the Kashiwara raising operator
        `e_i` on ``self``.

        This will remove the `i`-colored box corresponding to the
        rightmost `+` in ``self.signature(i)``.

        EXAMPLES::

            sage: x = crystals.infinity.GeneralizedYoungWalls(3)([[],[1,0,3,2],[2,1],[3,2,1,0,3,2],[],[],[2]])
            sage: x.e(2)
            [[], [1, 0, 3, 2], [2, 1], [3, 2, 1, 0, 3, 2]]
            sage: _.e(2)
            [[], [1, 0, 3], [2, 1], [3, 2, 1, 0, 3, 2]]
            sage: _.e(2)
            [[], [1, 0, 3], [2, 1], [3, 2, 1, 0, 3]]
            sage: _.e(2)
        """
        signature = self.generate_signature(i)
        raw_signature = signature[0]
        lastminus = signature[1].rfind('-')
        newdata = []
        if lastminus > -1:
            deletionrow = raw_signature[lastminus][1]
            for r in range(self.rows):
                if r == deletionrow:
                    newdata.append(list(self.data[r][:-1]))
                else:
                    newdata.append(list(self.data[r]))
            return self.__class__(self.parent(), newdata)
        else:
            return None

    def f(self, i):
        r"""
        Return the application of the Kashiwara lowering operator
        `f_i` on ``self``.

        This will add an `i`-colored colored box to the site corresponding
        to the leftmost plus in ``self.signature(i)``.

        EXAMPLES::

            sage: hw = crystals.infinity.GeneralizedYoungWalls(2)([])
            sage: hw.f(1)
            [[], [1]]
            sage: _.f(2)
            [[], [1], [2]]
            sage: _.f(0)
            [[], [1, 0], [2]]
            sage: _.f(0)
            [[0], [1, 0], [2]]
        """
        signature = self.generate_signature(i)
        raw_signature = signature[0]
        firstplus = signature[1].find('+')
        newdata = deepcopy(self.data)
        if firstplus > -1:
            additionrow = raw_signature[firstplus][1]
            newdata[additionrow].append(i)
        else:
            while len(newdata) % self.cartan_type().rank() != i:
                newdata.append([])
            newdata.append([i])
        return self.__class__(self.parent(), newdata)

    def latex_large(self):
        r"""
        Generate LaTeX code for ``self`` but the output is larger.
        Requires TikZ.

        EXAMPLES::

            sage: x = crystals.infinity.GeneralizedYoungWalls(3)([[],[1,0,3,2],[2,1],[3,2,1,0,3,2],[],[],[2]])
            sage: x.latex_large()
            '\\begin{tikzpicture}[baseline=5,scale=.45] \n \\foreach \\x [count=\\s from 0] in \n{{},{1,0,3,2},{2,1},{3,2,1,0,3,2},{},{},{2}} \n{\\foreach \\y [count=\\t from 0] in \\x {  \\node[font=\\scriptsize] at (-\\t,\\s) {$\\y$}; \n \\draw (-\\t+.5,\\s+.5) to (-\\t-.5,\\s+.5); \n \\draw (-\\t+.5,\\s-.5) to (-\\t-.5,\\s-.5); \n \\draw (-\\t-.5,\\s-.5) to (-\\t-.5,\\s+.5);  } \n \\draw[-,thick] (.5,\\s+1) to (.5,-.5) to (-\\t-1,-.5); } \n \\end{tikzpicture} \n'
        """
        s = ""
        if not self.data:
            s += "\\emptyset"
        else:
            s += "\\begin{tikzpicture}[baseline=5,scale=.45] \n \\foreach \\x [count=\\s from 0] in \n"
            s += "{" + ','.join("{" + ','.join(str(i) for i in r) + "}"
                                for r in self.data) + "} \n"
            s += "{\\foreach \\y [count=\\t from 0] in \\x {  \\node[font=\\scriptsize] at (-\\t,\\s) {$\\y$}; \n \\draw (-\\t+.5,\\s+.5) to (-\\t-.5,\\s+.5); \n \\draw (-\\t+.5,\\s-.5) to (-\\t-.5,\\s-.5); \n \\draw (-\\t-.5,\\s-.5) to (-\\t-.5,\\s+.5);  } \n \\draw[-,thick] (.5,\\s+1) to (.5,-.5) to (-\\t-1,-.5); } \n \\end{tikzpicture} \n"
        return s

    def _latex_(self):
        r"""
        Generate LaTeX code for ``self``.  Requires TikZ.

        EXAMPLES::

            sage: x = crystals.infinity.GeneralizedYoungWalls(3)([[],[1,0,3,2],[2,1],[3,2,1,0,3,2],[],[],[2]])
            sage: x._latex_()
            '\\begin{tikzpicture}[baseline=5,scale=.25] \\foreach \\x [count=\\s from 0] in \n{{},{1,0,3,2},{2,1},{3,2,1,0,3,2},{},{},{2}} \n{\\foreach \\y [count=\\t from 0] in \\x {  \\node[font=\\tiny] at (-\\t,\\s) {$\\y$}; \n \\draw (-\\t+.5,\\s+.5) to (-\\t-.5,\\s+.5); \n \\draw (-\\t+.5,\\s-.5) to (-\\t-.5,\\s-.5); \n \\draw (-\\t-.5,\\s-.5) to (-\\t-.5,\\s+.5);  } \n \\draw[-] (.5,\\s+1) to (.5,-.5) to (-\\t-1,-.5); } \n \\end{tikzpicture} \n'
        """
        s = ""
        if not self.data:
            s += "\\emptyset"
        else:
            s += "\\begin{tikzpicture}[baseline=5,scale=.25] \\foreach \\x [count=\\s from 0] in \n"
            s += "{" + ','.join("{" + ','.join(str(i) for i in r) + "}"
                                for r in self.data) + "} \n"
            s += "{\\foreach \\y [count=\\t from 0] in \\x {  \\node[font=\\tiny] at (-\\t,\\s) {$\\y$}; \n \\draw (-\\t+.5,\\s+.5) to (-\\t-.5,\\s+.5); \n \\draw (-\\t+.5,\\s-.5) to (-\\t-.5,\\s-.5); \n \\draw (-\\t-.5,\\s-.5) to (-\\t-.5,\\s+.5);  } \n \\draw[-] (.5,\\s+1) to (.5,-.5) to (-\\t-1,-.5); } \n \\end{tikzpicture} \n"
        return s

    def weight(self, root_lattice=False):
        r"""
        Return the weight of ``self``.

        INPUT:

        - ``root_lattice`` -- boolean determining whether weight should appear
          in root lattice or not in extended affine weight lattice

        EXAMPLES::

            sage: x = crystals.infinity.GeneralizedYoungWalls(3)([[],[1,0,3,2],[2,1],[3,2,1,0,3,2],[],[],[2]])
            sage: x.weight()
            2*Lambda[0] + Lambda[1] - 4*Lambda[2] + Lambda[3] - 2*delta
            sage: x.weight(root_lattice=True)
            -2*alpha[0] - 3*alpha[1] - 5*alpha[2] - 3*alpha[3]
        """
        L = self.cartan_type().root_system().root_lattice()
        alpha = L.simple_roots()
        W = sum(-1 * alpha[i] for r in self.data for i in r)
        if not root_lattice:
            E = self.cartan_type().root_system().weight_lattice(extended=True)
            return E(W)
        return L(W)

    def epsilon(self, i):
        r"""
        Return the number of `i`-colored arrows in the `i`-string above
        ``self`` in the crystal graph.

        EXAMPLES::

            sage: y = crystals.infinity.GeneralizedYoungWalls(3)([[],[1,0,3,2],[2,1],[3,2,1,0,3,2],[],[],[2]])
            sage: y.epsilon(1)
            0
            sage: y.epsilon(2)
            3
            sage: y.epsilon(0)
            0
        """
        if i not in self.index_set():
            raise ValueError("i must be in the index set")
        eps = 0
        while True:
            self = self.e(i)
            if self is None:
                break
            eps += 1
        return eps

    def Epsilon(self):
        r"""
        Return `\sum_{i=0}^n \varepsilon_i(Y) \Lambda_i` where `Y` is ``self``.

        EXAMPLES::

            sage: y = crystals.infinity.GeneralizedYoungWalls(3)([[0],[1,0,3,2],[2,1],[3,2,1,0,3,2],[0],[],[2]])
            sage: y.Epsilon()
            Lambda[0] + 3*Lambda[2]
        """
        La = self.cartan_type().root_system().weight_lattice().fundamental_weights()
        return sum(self.epsilon(i) * La[i] for i in self.index_set())

    def phi(self, i):
        r"""
        Return the value `\varepsilon_i(Y) + \langle h_i,
        \mathrm{wt}(Y)\rangle`, where `h_i` is the `i`-th simple
        coroot and `Y` is ``self``.

        EXAMPLES::

            sage: y = crystals.infinity.GeneralizedYoungWalls(3)([[0],[1,0,3,2],[2,1],[3,2,1,0,3,2],[0],[],[2]])
            sage: y.phi(1)
            3
            sage: y.phi(2)
            -1
        """
        h = self.parent().weight_lattice_realization().simple_coroots()
        return self.epsilon(i) + self.weight(root_lattice=False).scalar(h[i])

    def Phi(self):
        r"""
        Return `\sum_{i=0}^n \varphi_i(Y) \Lambda_i` where `Y` is ``self``.

        EXAMPLES::

            sage: y = crystals.infinity.GeneralizedYoungWalls(3)([[0],[1,0,3,2],[2,1],[3,2,1,0,3,2],[0],[],[2]])
            sage: y.Phi()
            -Lambda[0] + 3*Lambda[1] - Lambda[2] + 3*Lambda[3]

            sage: x = crystals.infinity.GeneralizedYoungWalls(3)([[],[1,0,3,2],[2,1],[3,2,1,0,3,2],[],[],[2]])
            sage: x.Phi()
            2*Lambda[0] + Lambda[1] - Lambda[2] + Lambda[3]
        """
        La = self.cartan_type().root_system().weight_lattice(extended=True).fundamental_weights()
        return sum(self.phi(i) * La[i] for i in self.index_set())

    def column(self, k):
        r"""
        Return the list of boxes from the ``k``-th column of ``self``.

        EXAMPLES::

            sage: y = crystals.infinity.GeneralizedYoungWalls(3)([[0],[1,0,3,2],[2,1],[3,2,1,0,3,2],[0],[],[2]])
            sage: y.column(2)
            [None, 0, 1, 2, None, None, None]

            sage: hw = crystals.infinity.GeneralizedYoungWalls(5)([])
            sage: hw.column(1)
            []
        """
        return [row[k - 1] if k - 1 < len(row) else None
                for row in self.data]

    def a(self, i, k):
        r"""
        Return the number `a_i(k)` of `i`-colored boxes in the ``k``-th
        column of ``self``.

        EXAMPLES::

            sage: y = crystals.infinity.GeneralizedYoungWalls(3)([[0],[1,0,3,2],[2,1],[3,2,1,0,3,2],[0],[],[2]])
            sage: y.a(1,2)
            1
            sage: y.a(0,2)
            1
            sage: y.a(3,2)
            0
        """
        A = [1 for c in range(len(self.column(k)))
             if self.column(k)[c] == i]
        return len(A)

    def in_highest_weight_crystal(self, La):
        r"""
        Return a boolean indicating if the generalized Young wall element
        is in the highest weight crystal cut out by the given highest weight
        ``La``.

        By Theorem 4.1 of [KS2010]_, a generalized Young wall `Y` represents a
        vertex in the highest weight crystal `Y(\lambda)`, with
        `\lambda = \Lambda_{i_1} + \Lambda_{i_2} + \cdots + \Lambda_{i_\ell}`
        a dominant integral weight of level `\ell > 0`, if it satisfies the
        following condition. For each positive integer `k`, if there exists
        `j \in I` such that `a_j(k) - a_{j-1}(k) > 0`, then for some
        `p = 1, \ldots, \ell`,

        .. MATH::

            j + k \equiv i_p + 1 \bmod n+1 \text{ and } a_j(k) - a_{j-1}(k)
            \le \lambda(h_{i_p}),

        where `\{h_0, h_1, \ldots, h_n\}` is the set of simple coroots attached
        to `A_n^{(1)}`.

        EXAMPLES::

            sage: La = RootSystem(['A',2,1]).weight_lattice(extended=True).fundamental_weights()[1]
            sage: GYW = crystals.infinity.GeneralizedYoungWalls(2)
            sage: y = GYW([[],[1,0],[2,1]])
            sage: y.in_highest_weight_crystal(La)
            True
            sage: x = GYW([[],[1],[2],[],[],[2],[],[],[2]])
            sage: x.in_highest_weight_crystal(La)
            False
        """
        if La not in self.parent().weight_lattice_realization():
            raise TypeError("Must be an element in the weight lattice realization")
        ac = self.parent().weight_lattice_realization().simple_coroots()
        n = self.cartan_type().classical().rank()
        index_set = self.index_set()
        for k in range(1, self.cols + 1):
            for j in index_set:
                if self.a(j, k) - self.a((j - 1) % (n + 1), k) <= 0:
                    continue
                else:
                    p_not_found = True
                    for p in index_set:
                        if (j + k - p - 1) % (n + 1) == 0 and self.a(j, k) - self.a((j - 1) % (n + 1), k) <= La.scalar(ac[p]):
                            p_not_found = False
                            continue
                        else:
                            continue
                    if p_not_found:
                        return False
        return True


class InfinityCrystalOfGeneralizedYoungWalls(UniqueRepresentation, Parent):
    r"""
    The crystal `\mathcal{Y}(\infty)` of generalized Young walls of
    type `A_n^{(1)}` as defined in [KS2010]_.

    A generalized Young wall is a collection of boxes stacked on a fixed board,
    such that color of the box at the site located in the `j`-th row from the
    bottom and the `i`-th column from the right is `j-1 \bmod n+1`.  There are
    several growth conditions on elements in `Y \in \mathcal{Y}(\infty)`:

    - Walls grow in rows from right to left.  That is, for every box `y\in Y`
      that is not in the rightmost column, there must be a box immediately to
      the right of `y`.

    - For all `p>q` such that `p-q \equiv 0 \bmod n+1`, the `p`-th row has
      most as many boxes as the `q`-th row.

    - There does not exist a column in the wall such that if one `i`-colored
      box, for every `i = 0,1,\ldots,n`, is removed from that column, then the
      result satisfies the above conditions.

    There is a crystal structure on `\mathcal{Y}(\infty)` defined as follows.
    Define maps

    .. MATH::

        e_i,\ f_i \colon \mathcal{Y}(\infty)
        \longrightarrow \mathcal{Y}(\infty) \sqcup \{0\}, \qquad
        \varepsilon_i,\ \varphi_i \colon \mathcal{Y}(\infty)
        \longrightarrow \ZZ, \qquad
        \mathrm{wt}\colon \mathcal{Y}(\infty) \longrightarrow
        \bigoplus_{i=0}^n \ZZ \Lambda_i \oplus \ZZ \delta,

    by

    .. MATH::

        \mathrm{wt}(Y) = -\sum_{i=0}^n m_i(Y) \alpha_i,

    where `m_i(Y)` is the number of `i`-boxes in `Y`, `\varepsilon_i(Y)`
    is the number of `-` in the `i`-signature of `Y`, and

    .. MATH::

        \varphi_i(Y)  = \varepsilon_i(Y) + \langle h_i, \mathrm{wt}(Y) \rangle.

    See :meth:`GeneralizedYoungWall.e()`, :meth:`GeneralizedYoungWall.f()`,
    and :meth:`GeneralizedYoungWall.signature()` for more about
    `e_i`, `f_i`, and `i`-signatures.


    INPUT:

    - ``n`` -- type `A_n^{(1)}`

    EXAMPLES::

        sage: Yinf = crystals.infinity.GeneralizedYoungWalls(3)
        sage: y = Yinf([[0],[1,0,3,2],[],[3,2,1],[0],[1,0]])
        sage: y.pp()
            0|1|
              0|
          1|2|3|
               |
        2|3|0|1|
              0|
        sage: y.weight(root_lattice=True)
        -4*alpha[0] - 3*alpha[1] - 2*alpha[2] - 2*alpha[3]
        sage: y.f(0)
        [[0], [1, 0, 3, 2], [], [3, 2, 1], [0], [1, 0], [], [], [0]]
        sage: y.e(0).pp()
            0|1|
               |
          1|2|3|
               |
        2|3|0|1|
              0|

    To display the crystal down to depth 3::

        sage: S = Yinf.subcrystal(max_depth=3)
        sage: G = Yinf.digraph(subset=S) # long time
        sage: view(G) # not tested
    """

    @staticmethod
    def __classcall_private__(cls, n, category=None):
        r"""
        Normalize input to ensure a unique representation.

        INPUT:

        - ``n`` -- type `A_n^{(1)}`

        EXAMPLES::

            sage: Yinf = crystals.infinity.GeneralizedYoungWalls(3)
            sage: Yinf2 = crystals.infinity.GeneralizedYoungWalls(int(3))
            sage: Yinf is Yinf2
            True
        """
        return super().__classcall__(cls, n, category)

    def __init__(self, n, category):
        r"""
        EXAMPLES::

            sage: Yinf = crystals.infinity.GeneralizedYoungWalls(3)
            sage: TestSuite(Yinf).run()
        """
        self._cartan_type = CartanType(['A', n, 1])
        if category is None:
            category = (HighestWeightCrystals(), InfiniteEnumeratedSets())
        Parent.__init__(self, category=category)
        self.module_generators = (self.element_class(self, []),)

    Element = GeneralizedYoungWall

    def _element_constructor_(self, data):
        r"""
        Construct an element of ``self`` from ``data``.

        INPUT:

        - ``data`` -- a multilist

        EXAMPLES::

            sage: GYW = crystals.infinity.GeneralizedYoungWalls(2)
            sage: y = GYW([[],[1,0],[2,1]]) # indirect doctest
            sage: y
            [[], [1, 0], [2, 1]]
        """
        return self.element_class(self, data)

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: Y = crystals.infinity.GeneralizedYoungWalls(4)
            sage: Y
            Crystal of generalized Young walls of type ['A', 4, 1]
        """
        return "Crystal of generalized Young walls of type {}".format(self._cartan_type)


########################
#  Highest weight GYW  #
########################

class CrystalOfGeneralizedYoungWallsElement(GeneralizedYoungWall):
    r"""
    Element of the highest weight crystal of generalized Young walls.
    """

    def e(self, i):
        r"""
        Compute the action of `e_i` restricted to the highest weight crystal.

        EXAMPLES::

            sage: La = RootSystem(['A',2,1]).weight_lattice(extended=True).fundamental_weights()[1]
            sage: hwy = crystals.GeneralizedYoungWalls(2,La)([[],[1,0],[2,1]])
            sage: hwy.e(1)
            [[], [1, 0], [2]]
            sage: hwy.e(2)
            sage: hwy.e(3)
        """
        ret = GeneralizedYoungWall.e(self, i)
        if ret is None:
            return None
        if ret.in_highest_weight_crystal(self.parent().hw):
            return self.__class__(self.parent(), ret.data)
        return None

    def f(self, i):
        r"""
        Compute the action of `f_i` restricted to the highest weight crystal.

        EXAMPLES::

            sage: La = RootSystem(['A',2,1]).weight_lattice(extended=True).fundamental_weights()[1]
            sage: GYW = crystals.infinity.GeneralizedYoungWalls(2)
            sage: y = GYW([[],[1,0],[2,1]])
            sage: y.f(1)
            [[], [1, 0], [2, 1], [], [1]]
            sage: hwy = crystals.GeneralizedYoungWalls(2,La)([[],[1,0],[2,1]])
            sage: hwy.f(1)
        """
        ret = GeneralizedYoungWall.f(self, i)
        if ret.in_highest_weight_crystal(self.parent().hw):
            return self.__class__(self.parent(), ret.data)
        return None

    def weight(self):
        r"""
        Return the weight of ``self`` in the highest weight crystal as an
        element of the weight lattice `\bigoplus_{i=0}^n \ZZ \Lambda_i`.

        EXAMPLES::

            sage: La = RootSystem(['A',2,1]).weight_lattice(extended=True).fundamental_weights()[1]
            sage: hwy = crystals.GeneralizedYoungWalls(2,La)([[],[1,0],[2,1]])
            sage: hwy.weight()
            Lambda[0] - Lambda[1] + Lambda[2] - delta
        """
        return self.parent().weight_lattice_realization()(self.parent().hw + GeneralizedYoungWall.weight(self))

    def phi(self, i):
        r"""
        Return the value `\varepsilon_i(Y) + \langle h_i,
        \mathrm{wt}(Y)\rangle`, where `h_i` is the `i`-th simple
        coroot and `Y` is ``self``.

        EXAMPLES::

            sage: La = RootSystem(['A',3,1]).weight_lattice(extended=True).fundamental_weights()
            sage: y = crystals.GeneralizedYoungWalls(3,La[0])([])
            sage: y.phi(1)
            0
            sage: y.phi(2)
            0
        """
        h = self.parent().weight_lattice_realization().simple_coroots()
        return self.epsilon(i) + self.weight().scalar(h[i])


class CrystalOfGeneralizedYoungWalls(InfinityCrystalOfGeneralizedYoungWalls):
    r"""
    The crystal `\mathcal{Y}(\lambda)` of generalized Young walls of the given
    type with highest weight `\lambda`.

    These were characterized in Theorem 4.1 of [KS2010]_.
    See :meth:`GeneralizedYoungWall.in_highest_weight_crystal()`.

    INPUT:

    - ``n`` -- type `A_n^{(1)}`

    - ``weight`` -- dominant integral weight

    EXAMPLES::

        sage: La = RootSystem(['A',3,1]).weight_lattice(extended=True).fundamental_weights()[1]
        sage: YLa = crystals.GeneralizedYoungWalls(3,La)
        sage: y = YLa([[0],[1,0,3,2,1],[2,1,0],[3]])
        sage: y.pp()
                3|
            0|1|2|
        1|2|3|0|1|
                0|
        sage: y.weight()
        -Lambda[0] + Lambda[2] + Lambda[3] - 3*delta
        sage: y.in_highest_weight_crystal(La)
        True
        sage: y.f(1)
        [[0], [1, 0, 3, 2, 1], [2, 1, 0], [3], [], [1]]
        sage: y.f(1).f(1)
        sage: yy = crystals.infinity.GeneralizedYoungWalls(3)([[0], [1, 0, 3, 2, 1], [2, 1, 0], [3], [], [1]])
        sage: yy.f(1)
        [[0], [1, 0, 3, 2, 1], [2, 1, 0], [3], [], [1], [], [], [], [1]]
        sage: yyy = yy.f(1)
        sage: yyy.in_highest_weight_crystal(La)
        False

        sage: LS = crystals.LSPaths(['A',3,1],[1,0,0,0])
        sage: C = LS.subcrystal(max_depth=4)
        sage: G = LS.digraph(subset=C)
        sage: P = RootSystem(['A',3,1]).weight_lattice(extended=True)
        sage: La = P.fundamental_weights()
        sage: YW = crystals.GeneralizedYoungWalls(3,La[0])
        sage: CW = YW.subcrystal(max_depth=4)
        sage: GW = YW.digraph(subset=CW)
        sage: GW.is_isomorphic(G,edge_labels=True)
        True

    To display the crystal down to a specified depth::

        sage: S = YLa.subcrystal(max_depth=4)
        sage: G = YLa.digraph(subset=S)
        sage: view(G) # not tested
    """
    @staticmethod
    def __classcall_private__(cls, n, La):
        r"""
        EXAMPLES::

            sage: La = RootSystem(['A',2,1]).weight_lattice(extended=True).fundamental_weights()[2]
            sage: Al = RootSystem(['A',2,1]).weight_lattice(extended=True).monomial(2)
            sage: Y = crystals.GeneralizedYoungWalls(2,La)
            sage: Y1 = crystals.GeneralizedYoungWalls(int(2),Al)
            sage: Y is Y1
            True
        """
        La = RootSystem(['A', n, 1]).weight_lattice(extended=True)(La)
        return super().__classcall__(cls, n, La)

    def __init__(self, n, La):
        r"""
        EXAMPLES::

            sage: La = RootSystem(['A',2,1]).weight_lattice(extended=True).fundamental_weights()[1]
            sage: YLa = crystals.GeneralizedYoungWalls(2,La)

        We skip the two tests because they take a very long time::

            sage: TestSuite(YLa).run(skip=["_test_enumerated_set_contains","_test_stembridge_local_axioms"]) # long time
        """
        InfinityCrystalOfGeneralizedYoungWalls.__init__(self, n,
                category=(RegularCrystals(), HighestWeightCrystals(), InfiniteEnumeratedSets()))
        self.hw = La

    Element = CrystalOfGeneralizedYoungWallsElement

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: La = RootSystem(['A',5,1]).weight_lattice(extended=True).fundamental_weights()[2]
            sage: Y = crystals.GeneralizedYoungWalls(5,La)
            sage: Y
            Highest weight crystal of generalized Young walls of Cartan type ['A', 5, 1] and highest weight Lambda[2]
        """
        return "Highest weight crystal of generalized Young walls of Cartan type {1!s} and highest weight {0!s}".format(self.hw, self._cartan_type)

    def __iter__(self):
        r"""
        EXAMPLES::

            sage: y = crystals.infinity.GeneralizedYoungWalls(3)([[0],[1,0,3,2],[2,1],[3,2,1,0,3,2],[0],[],[2]])
            sage: x = y.__iter__()
            sage: next(x)
            [0]
        """
        for c in super().__iter__():
            if c.in_highest_weight_crystal(self.hw):
                yield c
