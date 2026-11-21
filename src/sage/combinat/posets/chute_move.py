r"""
Chute move lattices

Chute move posets were defined by Rubey in [Rub2011]_ as a
generalization of the chute and ladder posets of Bergeron and Billey
[BB1993]_.  They were shown to be lattices independently by Sara
Billey, Connor McCausland and Clare Minnerath in [BMM2025]_ and by
Ilani Axelrod-Freed, Colin Defant, Hanna Mularczyk, Son Nguyen and
Katherine Tung in [ADMNT2025]_.

TESTS:

We check Corollary 5.8 of [BMM2025]_::

    sage: def l0(w):
    ....:    return posets.ChuteMoveLattice(w).maximal_chain_length()

    sage: def R(x, y, w):
    ....:     return sum(1 for i in range(1, w.inverse()(y) + 1) if w(i) <= x)

    sage: def l1(w):
    ....:     return 1 + sum(R(x, y, w) for x, y in w.inverse().inversions())

    sage: def l2(w):
    ....:     return 1 + len(w.pattern_positions([1,3,2]))

    sage: all(l0(w) == l1(w) == l2(w) for n in range(7) for w in Permutations(n))
    True
"""
from sage.categories.finite_lattice_posets import FiniteLatticePosets
from sage.combinat.permutation import Permutation
from sage.combinat.posets.lattices import LatticePoset
from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet
from sage.structure.sage_object import SageObject


def intervals_to_polyomino(intervals) -> list[tuple[int, int]]:
    r"""
    Create a polyomino from a list of intervals.

    The intervals correspond to the rows of the polyomino.

    INPUT:

    - list of pairs of integers `(i, j)` with `i \leq j`

    OUTPUT:

    cells as a list of pairs of integers `(x, y)`

    EXAMPLES::

        sage: from sage.combinat.posets.chute_move import intervals_to_polyomino
        sage: intervals_to_polyomino([[0,4],[1,3],[2,2]])
        ((1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (2, 1), (2, 2), (2, 3), (3, 2))

    TESTS::

        sage: intervals_to_polyomino([[3,2]])
        Traceback (most recent call last):
        ...
        AssertionError: interval=[3, 2] is not a proper interval
    """
    M = set()
    for x, interval in enumerate(intervals, 1):
        low, high = interval
        assert low <= high, f"interval={interval} is not a proper interval"
        M.update((x, y) for y in range(low, high + 1))
    return tuple(sorted(M))


class PolyominoFilling(SageObject):
    r"""
    Class for polyominoes with some filled cells.

    A filling is represented by a list of matrix coordinates.
    """
    def __init__(self, P, B, check=True) -> None:
        r"""
        Initialize a filling of a polyomino.

        INPUT:

        - ``P`` -- an iterable of cells (i,j)

        - ``B`` -- a subset of ``P``

        - ``check`` -- boolean (default: ``True``) whether to
          transform ``P`` and ``Q`` into tuples of sorted coordinates

        EXAMPLES::

            sage: from sage.combinat.posets.chute_move import PolyominoFilling
            sage: PolyominoFilling([(1,1),(1,2),(2,1)], [(1,2)])
            shape=((1, 1), (1, 2), (2, 1)), filling=((1, 2),)
        """
        if check:
            self._P = tuple(sorted(P))
            self._B = tuple(sorted(B))
        else:
            self._P = P
            self._B = B

    def _repr_(self):
        r"""
        Return a text representation of a filling.

        EXAMPLES::

            sage: from sage.combinat.posets.chute_move import PolyominoFilling
            sage: PolyominoFilling([(1,1),(1,2),(2,1)], [(1,2)])
            shape=((1, 1), (1, 2), (2, 1)), filling=((1, 2),)
        """
        return f"shape={self._P}, filling={self._B}"

    def _array(self):
        r"""
        Return the filling as an array.

        EXAMPLES::

            sage: from sage.combinat.posets.chute_move import PolyominoFilling
            sage: P = PolyominoFilling([(1,1),(1,2),(2,1)], [(1,2)])
            sage: P._array()
            [['', 'o'], ['', None]]
        """
        xs = [i for i, _ in self._P]
        ys = [j for _, j in self._P]
        min_i, max_i = min(xs), max(xs)
        min_j, max_j = min(ys), max(ys)

        array = []
        for i in range(min_i, max_i + 1):
            row = []
            for j in range(min_j, max_j + 1):
                if (i, j) not in self._P:
                    row.append(None)
                elif (i, j) in self._B:
                    row.append("o")
                else:
                    row.append("")
            array.append(row)

        return array

    def _ascii_art_(self):
        r"""
        Return a pretty representation.

        EXAMPLES::

            sage: from sage.combinat.posets.chute_move import PolyominoFilling
            sage: P = PolyominoFilling([(1,1),(1,2),(2,1)], [(1,2)])
            sage: ascii_art(P)
            +---+---+
            |   | o |
            +---+---+
            |   |
            +---+
        """
        from sage.combinat.output import ascii_art_table
        from sage.typeset.ascii_art import AsciiArt
        return AsciiArt(ascii_art_table(self._array(), use_unicode=False).splitlines())

    def _unicode_art_(self):
        r"""
        Return a pretty representation.

        EXAMPLES::

            sage: from sage.combinat.posets.chute_move import PolyominoFilling
            sage: P = PolyominoFilling([(1,1),(1,2),(2,1)], [(1,2)])
            sage: unicode_art(P)
            ┌───┬───┐
            │   │ o │
            ├───┼───┘
            │   │
            └───┘
        """
        from sage.combinat.output import ascii_art_table
        from sage.typeset.unicode_art import UnicodeArt
        return UnicodeArt(ascii_art_table(self._array(), use_unicode=True).splitlines())

    def _latex_(self):
        r"""
        Return a LaTeX representation.

        EXAMPLES::

            sage: from sage.combinat.posets.chute_move import PolyominoFilling
            sage: P = PolyominoFilling([(1,1),(1,2),(2,1)], [(1,2)])
            sage: latex(P)
            \begin{tikzpicture}[x=0.5cm,y=0.5cm]
            ...
            \end{tikzpicture}
        """
        from sage.misc.latex import latex
        latex.add_package_to_preamble_if_available("tikz")

        xs = [i for i, _ in self._P]
        ys = [j for _, j in self._P]
        min_i = min(xs)
        max_i = max(xs)
        min_j = min(ys)
        height = max_i - min_i + 1

        cell_size = 0.5
        bullet = r"\bullet"

        lines = []
        lines.append(r"\begin{tikzpicture}[x=%scm,y=%scm]" % (cell_size, cell_size))

        for (i, j) in self._P:
            x = j - min_j
            y = height - (i - min_i) - 1

            lines.append(f"\\draw ({x},{y}) rectangle ({x+1},{y+1});")

            if (i, j) in self._B:
                lines.append(f"\\node at ({x+0.5},{y+0.5}) {{$ {bullet} $}};")

        lines.append(r"\end{tikzpicture}")
        return "\n".join(lines)

    def __hash__(self):
        r"""
        Return the hash of ``self``.

        EXAMPLES::

            sage: from sage.combinat.posets.chute_move import PolyominoFilling
            sage: P = PolyominoFilling([(1,1),(1,2),(2,1)], [(1,2)])
            sage: hash(P)  # random
            5139392050573932802
        """
        return hash((self._P, self._B))

    def __eq__(self, other):
        r"""
        Check whether ``self`` is equal to ``other``.

        EXAMPLES::

            sage: from sage.combinat.posets.chute_move import PolyominoFilling
            sage: P = PolyominoFilling([(1,1),(1,2),(2,1)], [(1,2)])
            sage: Q = PolyominoFilling([(1,1),(1,2),(2,1)], [(2,1)])
            sage: P == P
            True
            sage: P == Q
            False
        """
        return (isinstance(other, PolyominoFilling)
                and self._P == other._P
                and self._B == other._B)

    def __ne__(self, other):
        r"""
        Check whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: from sage.combinat.posets.chute_move import PolyominoFilling
            sage: P = PolyominoFilling([(1,1),(1,2),(2,1)], [(1,2)])
            sage: Q = PolyominoFilling([(1,1),(1,2),(2,1)], [(2,1)])
            sage: P != P
            False
            sage: P != Q
            True
        """
        return not (self == other)


def ChuteMoveLattice(M, n=None):
    r"""
    Return the chute move lattice.

    INPUT:

    - ``M`` -- a permutation, or an L-convex polyomino
    - ``n`` -- a positive integer, if ``M`` is a polyomino, or ``None``

    EXAMPLES::

        sage: from sage.combinat.posets.chute_move import intervals_to_polyomino
        sage: M = intervals_to_polyomino([(1,4),(1,3),(1,2),(1,1)])
        sage: L = posets.ChuteMoveLattice(M, 1); L
        Finite lattice containing 5 elements

        sage: from sage.combinat.tamari_lattices import GeneralizedTamariLattice
        sage: b = 4; m = 2; T = GeneralizedTamariLattice(m*b+1,b,m=2)
        sage: C = posets.ChuteMoveLattice(Permutation([1,8,6,4,2,3,5,7]))
        sage: C.is_isomorphic(T)
        True
    """
    from sage.groups.perm_gps.permgroup_element import SymmetricGroupElement
    if isinstance(M, SymmetricGroupElement):
        M = Permutation(M)

    if isinstance(M, Permutation):
        n = len(M)
        Minv = M.inverse()

        def above_left(j):
            return sum(M(i) > j for i in range(1, Minv(j)))

        top = tuple(sorted([(j, c)
                            for j in range(1, n + 1)
                            for c in range(1 + above_left(j), n + 2 - j)]))

        M = intervals_to_polyomino([(1, n - i) for i in range(n)])

    else:
        boundary_cells = [(x, y) for x, y in M
                          if (x-1, y) not in M
                          or (x, y+1) not in M
                          or (x-1, y+1) not in M]

        top = tuple(sorted(set((cx, cy)
                               for x, y in boundary_cells
                               for cx in range(x, x + n)
                               for cy in range(y - n + 1, y + 1)
                               if (cx, cy) in M)))

    def chutable(i, l):
        r"""
        Return the other coordinate, if ``l[i]`` is a chutable
        coordinate and the result is in the polyomino, otherwise
        return ``None``.

        The rectangle is given (in Cartesian coordinates) by

        # a,d - c,d
        #  |     |
        # a,b - c,b

        with `a < c` and `d < b`.

        `(a, b)` must not be the first or last element of ``l``, and
        the predecessor `(a, d)` of `(a, b)` must have the same `x`
        coordinate.
        """
        a, b = l[i]

        i -= 1
        a1, d = l[i]
        if a1 != a:
            return

        i += 2
        c, b1 = l[i]
        while i + 1 < len(l) and (c == a or b1 < d or b1 > b):
            i += 1
            c, b1 = l[i]

        if i >= len(l) or b1 != b or (c, d) not in M:
            return

        return c, d

    def children(l):
        r"""
        Return the children of a given filling.
        """
        return [tuple(sorted(l[:i] + (c,) + l[i+1:]))
                for i in range(1, len(l)-1)
                if (c := chutable(i, l)) is not None]

    S = RecursivelyEnumeratedSet([top], children,
                                 structure=None, enumeration="naive")

    F = PolyominoFilling
    d = {F(M, f, check=False): [F(M, g, check=False) for g in children(f)]
         for f in S}

    return LatticePoset(d, category=FiniteLatticePosets().CongruenceUniform())
