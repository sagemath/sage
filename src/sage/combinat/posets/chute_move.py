r"""
Chute move lattices

Chute move posets were defined by Rubey in [Rub2011]_.  They were
shown to be lattices independently by Sara C. Billey, Connor
McCausland, Clare Minnerath in [Bil2025]_ and by Ilani Axelrod-Freed,
Colin Defant, Hanna Mularczyk, Son Nguyen, Katherine Tung [Axe2025]_.
"""
from sage.categories.finite_lattice_posets import FiniteLatticePosets
from sage.combinat.permutation import Permutation
from sage.combinat.posets.lattices import LatticePoset
from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet
from sage.structure.sage_object import SageObject


def intervals_to_polyomino(intervals):
    M = set()
    for x, interval in enumerate(intervals, 1):
        low, high = interval
        assert low <= high, f"interval={interval} is not a proper interval"
        M.update((x, y) for y in range(low, high + 1))
    return sorted(M)


class PolyominoFilling(SageObject):
    def __init__(self, P, F):
        self._P = tuple(sorted(P))
        self._F = tuple(sorted(F))

    def __hash__(self):
        return hash((self._P, self._F))

    def __eq__(self, other):
        return self._P == other._P and self._F == other._F

    def _repr_(self):
        return f"shape={self._P}, filling={self._F}"

    def _latex_(self):
        draw_grid_lines = True,
        stroke_width = "" # very thin"
        bullet_tex = r"\bullet"
        cell_size = 0.5

        from sage.misc.latex import latex
        latex.add_package_to_preamble_if_available("tikz")

        # bounding box
        xs = [x for x,_ in self._P]
        ys = [y for _,y in self._P]
        minx, maxx = min(xs), max(xs)
        miny, maxy = min(ys), max(ys)

        # shift so lower-left is (0,0)
        shift_x = -minx
        shift_y = -miny
        width = maxx - minx + 1
        height = maxy - miny + 1

        # Build TikZ code
        tikz_lines = []
        tikz_lines.append(r"\begin{tikzpicture}[x=%scm,y=%scm]" % (cell_size, cell_size))

        # Optionally draw a faint background rectangle for the whole bounding box
        tikz_lines.append(r"  % bounding box (for reference, not stroked)")
        tikz_lines.append(r"  \path[use as bounding box] (0,0) rectangle (%d,%d);" % (width, height))
        tikz_lines.append("")

        # Draw each cell as a unit square at shifted coordinates
        if draw_grid_lines:
            tikz_lines.append(r"  % draw cells")
            for (x, y) in self._P:
                sx, sy = x + shift_x, y + shift_y
                tikz_lines.append(f"  \\draw[{stroke_width}] ({sx},{sy}) rectangle ({sx+1},{sy+1});")
        else:
            tikz_lines.append(r"  % no grid lines; draw filled white squares to define shape")
            for (x, y) in self._P:
                sx, sy = x + shift_x, y + shift_y
                tikz_lines.append(f"  \\fill[white] ({sx},{sy}) rectangle ({sx+1},{sy+1});")

        tikz_lines.append("")

        # Place bullets (or arbitrary TeX) at the center of cells F
        tikz_lines.append(r"  % bullets in selected cells")
        for (x, y) in self._F:
            if (x, y) not in self._P:
                # skip or optionally raise error: here we skip invalid bullets
                continue
            sx, sy = x + shift_x, y + shift_y
            # place a node at the center; use baseline to center math mode bullet nicely
            tikz_lines.append(f"  \\node at ({sx+0.5},{sy+0.5}) {{${bullet_tex}$}};")

        tikz_lines.append(r"\end{tikzpicture}")
        return "\n".join(tikz_lines)


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

        top = tuple([(j, c)
                     for j in range(1, n + 1)
                     for c in range(1 + above_left(j), n + 2 - j)])
        M = intervals_to_polyomino([(1, n - i) for i in range(n)])
    else:
        boundary_cells = [(x, y) for x, y in M
                          if (x-1, y) not in M or (x, y+1) not in M or (x-1, y+1) not in M]
        top = tuple(sorted(set((cx, cy)
                               for x, y in boundary_cells
                               for cx in range(x, x + n)
                               for cy in range(y - n + 1, y + 1)
                               if (cx, cy) in M)))

    def chutable(i, l):
        # Return the other coordinate, if ``l[i]`` is a chutable
        # coordinate and the result is in ``M``, otherwise ``False``.

        # the rectangle is given (in Cartesian coordinates) by
        # a,d - c,d
        #  |     |
        # a,b - c,b
        # with a < c and d < b
        a, b = l[i]
        # (a, b) must not be the first or last element of l
        # the predecessor (a, d) of (a, b) must have same x coordinate
        i -= 1
        a1, d = l[i]
        if a1 != a:
            return

        # the first element within the rectangle must be (c, b)
        i += 2
        c, b1 = l[i]
        while i + 1 < len(l) and (c == a or b1 < d or b1 > b):
            i += 1
            c, b1 = l[i]

        if i >= len(l) or b1 != b or (c, d) not in M:
            return

        return c, d

    def children(l):
        return [tuple(sorted(l[:i] + (c,) + l[i+1:]))
                for i in range(1, len(l)-1)
                if (c := chutable(i, l)) is not None]

    S = RecursivelyEnumeratedSet([top], children,
                                 structure=None, enumeration="naive")
    d = {PolyominoFilling(M, f): [PolyominoFilling(M, g) for g in children(f)]
         for f in S}
    cat = FiniteLatticePosets().CongruenceUniform()
    return LatticePoset(d, category=cat)
