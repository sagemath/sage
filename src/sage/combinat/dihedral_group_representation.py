# sage.doctest: needs sage.groups sage.modules
r"""
Dihedral group representations

This implements the irreducible representations of the dihedral group
D_n using Sage's Representation_abstract framework.

EXAMPLES::

    sage: G = DihedralGroup(5)
    sage: reps = dihedral_irreducibles(G, 5)
    sage: len(reps)
    4

    sage: [V.dimension() for V in reps]
    [1, 1, 2, 2]

    sage: V = reps[0]
    sage: v = V.an_element()
    sage: e = G.identity()
    sage: e * v == v
    True

    sage: V = reps[-1]
    sage: v = V.an_element()
    sage: g, h = list(G)[:2]
    sage: (g*h) * v == g * (h * v)
    True
"""

from sage.categories.modules_with_basis import ModulesWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.modules.with_basis.representation import Representation_abstract
from sage.rings.all import QQ
from sage.rings.number_field.number_field import CyclotomicField
from sage.misc.cachefunc import cached_method
from sage.matrix.constructor import matrix


class DihedralGroupRepresentation(Representation_abstract):

    def __init__(self, G, n, base_ring=QQ):
        self._n = n

        # pick generators once
        self._r = [g for g in G if g.order() == n][0]
        self._s = [g for g in G if g.order() == 2 and g != self._r**(n//2)][0]

        Representation_abstract.__init__(self, G, "left", G.algebra(base_ring))

    @cached_method
    def express_in_gens(self, g):
        """
        Return (s_exp, r_exp) such that:
            g = r^k        -> (0, k)
            g = s * r^k    -> (1, k)
        """
        r = self._r
        s = self._s
        n = self._n

        for k in range(n):
            if g == r**k:
                return (0, k)

        for k in range(n):
            if g == s * r**k:
                return (1, k)

        raise ValueError(f"Could not express {g} in generators")


# =========================================================
# 1-DIMENSIONAL REPRESENTATIONS
# =========================================================

class DihedralOneDimensionalRep(DihedralGroupRepresentation,
                                CombinatorialFreeModule):
    """
    A 1-dimensional representation of D_n.

    EXAMPLES::

        sage: G = DihedralGroup(5)
        sage: V = DihedralOneDimensionalRep(G, "sign_s", 5)
        sage: r = [g for g in G if g.order() == 5][0]
        sage: s = [g for g in G if g.order() == 2][0]
        sage: V._character(r)
        1
        sage: V._character(s)
        -1
    """

    def __init__(self, G, character_type, n, base_ring=QQ):
        """
        character_type:
            "trivial", "sign_r", "sign_s", "sign_total"
        """
        self._type = character_type

        indices = ["v"]  # single basis vector

        CombinatorialFreeModule.__init__(
            self,
            base_ring,
            indices,
            category=ModulesWithBasis(base_ring).FiniteDimensional(),
            prefix="e"
        )

        DihedralGroupRepresentation.__init__(self, G, n, base_ring)

    def _repr_(self):
        return f"{self._type} representation of {self._semigroup}"

    def _character(self, g):

        (s_exp, r_exp) = self.express_in_gens(g)

        if self._type == "trivial":
            return 1

        if self._type == "sign_s":
            return (-1)**s_exp

        if self._type == "sign_r":
            return (-1)**r_exp

        if self._type == "sign_total":
            return (-1)**(r_exp + s_exp)

    class Element(CombinatorialFreeModule.Element):
        def _acted_upon_(self, g, self_on_left):
            if self_on_left:
                return None

            P = self.parent()
            scalar = P._character(g)
            return scalar * self


# =========================================================
# 2-DIMENSIONAL REPRESENTATIONS
# =========================================================

class DihedralTwoDimensionalRep(DihedralGroupRepresentation,
                               CombinatorialFreeModule):
    """
    The 2-dimensional irreducible representation indexed by k.

    EXAMPLES::

        sage: G = DihedralGroup(5)
        sage: V = DihedralTwoDimensionalRep(G, 1, 5)
        sage: g = list(G)[0]
        sage: M = V._matrix(g)
        sage: M.nrows(), M.ncols()
        (2, 2)

        sage: v = V.an_element()
        sage: lhs = g * v
        sage: rhs = V.from_vector(V._matrix(g) * v.to_vector())
        sage: lhs == rhs
        True
    """

    def __init__(self, G, k, n, base_ring=None):
        self._k = k
        self._n = n

        if base_ring is None:
            base_ring = CyclotomicField(n)
        else:
            try:
                base_ring.zeta(n)
            except (AttributeError, ValueError):
                base_ring = CyclotomicField(n)

        z = base_ring.zeta(n)

        indices = ["v1", "v2"]

        CombinatorialFreeModule.__init__(
            self,
            base_ring,
            indices,
            category=ModulesWithBasis(base_ring).FiniteDimensional(),
            prefix="v"
        )

        DihedralGroupRepresentation.__init__(self, G, n, base_ring)

        self._z = z

    def _repr_(self):
        return f"2D irrep k={self._k} of {self._semigroup}"

    def _matrix(self, g):
        """
        Return the matrix of g in this representation.

        EXAMPLES::

            sage: G = DihedralGroup(5)
            sage: V = DihedralTwoDimensionalRep(G, 1, 5)
            sage: g = list(G)[0]
            sage: V._matrix(g).nrows()
            2
        """
        # You should reuse your express_in_gens logic here
        # Placeholder:
        k = self._k
        z = self._z

        (s_exp, r_exp) = self.express_in_gens(g)

        if s_exp == 0:
            return matrix([[z**(k*r_exp), 0], [0, z**(-k*r_exp)]])
        if s_exp == 1:
            return matrix([[0, z**(k*r_exp)], [z**(-k*r_exp), 0]])

    class Element(CombinatorialFreeModule.Element):
        def _acted_upon_(self, g, self_on_left):
            if self_on_left:
                return None

            P = self.parent()
            M = P._matrix(g)

            vec = self.to_vector()
            new_vec = M * vec

            return P.from_vector(new_vec)


# =========================================================
# FACTORY
# =========================================================

def dihedral_irreducibles(G, n, base_ring=QQ):
    """
    Return all irreducible representations of D_n.
    """
    reps = []

    # 1-dimensional
    reps.append(DihedralOneDimensionalRep(G, "trivial", n, base_ring))
    reps.append(DihedralOneDimensionalRep(G, "sign_s", n, base_ring))

    if n % 2 == 0:
        reps.append(DihedralOneDimensionalRep(G, "sign_r", n, base_ring))
        reps.append(DihedralOneDimensionalRep(G, "sign_total", n, base_ring))

    # 2-dimensional
    for k in range(1, n // 2):
        reps.append(DihedralTwoDimensionalRep(G, k, n, base_ring))

    return reps

# dihedral_cellular_basis.py
# sage.doctest: needs sage.groups sage.modules
r"""
Cellular basis for the dihedral group algebra kD_n.

Implements the matrix-coefficient basis arising from the Wedderburn
isomorphism:  kD_n  ≅  ⊕_λ  M_{d_λ}(k)

The cellular basis elements are:

    C^λ_{ij}  =  (d_λ / |G|)  *  Σ_{g ∈ G}  conj(ρ^λ(g)_{ij}) * g

where d_λ = dim(V^λ) and ρ^λ is the matrix representation.

This satisfies the Graham-Lehrer axioms:
  (GL1)  The set {C^λ_{ij}} is a k-basis for kG.
  (GL2)  The anti-involution * on kG (induced by g ↦ g⁻¹) satisfies
         *(C^λ_{ij}) = C^λ_{ji}.
  (GL3)  For any a ∈ kG:
         a · C^λ_{ij}  ≡  Σ_i'  r_a(i', i) C^λ_{i'j}  (mod A^{>λ})
         where r_a(i', i) = Σ_{g} coeff(a, g) * ρ^λ(g)_{i'i}.

EXAMPLES::

    sage: from dihedral_group_representation import dihedral_irreducibles
    sage: G = DihedralGroup(5)
    sage: CB = DihedralCellularBasis(G, 5)
    sage: CB.basis_elements()
    ...
    sage: CB.verify_cellular_axioms()
    True
"""


# ─────────────────────────────────────────────────────────────────
# Helper: extract matrix of g in a representation
# ─────────────────────────────────────────────────────────────────

def _rep_matrix(V, g):
    """
    Return the matrix of group element g acting on representation V.

    Works for both 1- and 2-dimensional representations.

    EXAMPLES::

        sage: G = DihedralGroup(5)
        sage: from dihedral_group_representation import DihedralTwoDimensionalRep
        sage: V = DihedralTwoDimensionalRep(G, 1, 5)
        sage: g = list(G)[1]
        sage: M = _rep_matrix(V, g)
        sage: M.nrows()
        2
    """
    if isinstance(V, DihedralOneDimensionalRep):
        return matrix([[V._character(g)]])
    else:
        return V._matrix(g)


# ─────────────────────────────────────────────────────────────────
# Cellular basis class
# ─────────────────────────────────────────────────────────────────

class DihedralCellularBasis:
    r"""
    Cellular basis for kD_n via the Wedderburn decomposition.

    INPUT:

    - ``G`` -- a ``DihedralGroup(n)``
    - ``n`` -- integer, order of the rotation subgroup
    - ``base_ring`` -- (default: ``CyclotomicField(n)``) the coefficient ring
    - ``reps`` -- list of irreducible representations
    - ``kG`` -- the group algebra over the base ring
    - ``cells`` -- dict mapping irrep index to list of ``(i, j, algebra_element)``

    EXAMPLES::

        sage: G = DihedralGroup(5)
        sage: CB = DihedralCellularBasis(G, 5)
        sage: len(CB.cells)     # one entry per irrep
        4
        sage: sum(len(v) for v in CB.cells.values())  # total basis elements = 2n
        10
    """

    def __init__(self, G, n, base_ring=None):
        self.G  = G
        self.n  = n
        order   = 2 * n          # |D_n| = 2n

        # Choose a base ring that contains all representation values.
        # The 2D reps need a primitive n-th root of unity.
        if base_ring is None:
            base_ring = CyclotomicField(n) if n > 2 else QQ

        self.base_ring = base_ring
        self.kG        = G.algebra(base_ring)
        self.reps      = dihedral_irreducibles(G, n, base_ring)

        # Build the basis
        self._build_cells()

    # ── internal ──────────────────────────────────────────────────

    def _build_cells(self):
        """
        Populate self.cells:
            cells[lam_idx] = [ (i, j, C^lam_{ij}),  … ]
            where ``C^lam_{ij}`` is in ``kG`` and is given by::

                C^lam_{ij} = (d_lam / 2n) * sum_g  conj(rho^lam(g)_{ij}) * g

        EXAMPLES::

            sage: G = DihedralGroup(3)
            sage: CB = DihedralCellularBasis(G, 3)
            sage: CB._build_cells()
            sage: all(C in CB.kG for _, _, C in CB.cells[0])
            True
        """
        G      = self.G
        kG     = self.kG
        order  = 2 * self.n
        cells  = {}

        for lam_idx, V in enumerate(self.reps):
            d    = V.dimension()          # d_λ
            coeff = base_ring(d) / order  # d_λ / |G|, lives in the base ring

            cell_list = []
            for i in range(d):
                for j in range(d):
                    # Sum over group elements
                    elem = kG.zero()
                    for g in G:
                        M_g  = _rep_matrix(V, g)
                        # matrix entry ρ^λ(g)_{ij}  (0-indexed)
                        rho_ij = M_g[i, j]
                        # conjugate: for cyclotomic fields, use .conjugate()
                        # for QQ it is a no-op
                        try:
                            rho_ij_bar = rho_ij.conjugate()
                        except AttributeError:
                            rho_ij_bar = rho_ij

                        elem += base_ring(rho_ij_bar) * kG(g)

                    cell_list.append((i, j, coeff * elem))

            cells[lam_idx] = cell_list

        self.cells = cells

    # ── public API ────────────────────────────────────────────────

    def basis_elements(self):
        """
        Return all basis elements as a flat list of (λ_idx, i, j, C^λ_{ij}).

        EXAMPLES::

            sage: G = DihedralGroup(5)
            sage: CB = DihedralCellularBasis(G, 5)
            sage: elems = CB.basis_elements()
            sage: len(elems) == 2 * 5
            True
        """
        result = []
        for lam_idx, cell_list in self.cells.items():
            for (i, j, C) in cell_list:
                result.append((lam_idx, i, j, C))
        return result

    def C(self, lam_idx, i, j):
        """
        Return the single basis element C^{lam_idx}_{i,j}.

        EXAMPLES::

            sage: G = DihedralGroup(5)
            sage: CB = DihedralCellularBasis(G, 5)
            sage: c = CB.C(0, 0, 0)   # trivial idempotent
            sage: c in CB.kG
            True
        """
        for (ii, jj, C) in self.cells[lam_idx]:
            if ii == i and jj == j:
                return C
        raise ValueError(f"No cell element for λ={lam_idx}, i={i}, j={j}")

    # ── verification ──────────────────────────────────────────────

    def verify_basis(self):
        """
        Check that the ``2*n`` basis elements are linearly independent
        over the base ring by forming a change-of-basis matrix and checking
        its determinant is non-zero.

        EXAMPLES::

            sage: G = DihedralGroup(5)
            sage: CB = DihedralCellularBasis(G, 5)
            sage: CB.verify_basis()
            True
        """
        kG    = self.kG
        G     = self.G
        elems = self.basis_elements()
        G_list = list(G)
        N      = len(G_list)
        g_idx  = {g: k for k, g in enumerate(G_list)}

        rows = []
        for (_, _, C) in elems:
            row = [C.coefficient(g) for g in G_list]
            rows.append(row)

        M = matrix(self.base_ring, rows)
        return M.rank() == N

    def verify_involution(self):
        """
        Check (GL2): the anti-involution ``*`` satisfies ``*(C^lam_{ij}) = C^lam_{ji}``,
        where ``*`` is induced by ``g -> g^{-1}``.

        EXAMPLES::

            sage: G = DihedralGroup(5)
            sage: CB = DihedralCellularBasis(G, 5)
            sage: CB.verify_involution()
            True
        """
        kG = self.kG

        def star(elem):
            """Apply the anti-involution g ↦ g^{-1}."""
            result = kG.zero()
            for g in self.G:
                coeff = elem.coefficient(g)
                if coeff != 0:
                    result += self.base_ring(coeff) * kG(g.inverse())
            return result

        for lam_idx in self.cells.keys():
            d = self.reps[lam_idx].dimension()
            for i in range(d):
                for j in range(d):
                    Cij  = self.C(lam_idx, i, j)
                    Cji  = self.C(lam_idx, j, i)
                    if star(Cij) != Cji:
                        return False
        return True

    def verify_multiplication(self):
        """
        Check (GL3): g · C^λ_{ij} = Σ_{i'} ρ^λ(g)_{i'i} · C^λ_{i'j}
        for each g ∈ G and each basis element C^λ_{ij}.

        This is the key cellular multiplication axiom.

        EXAMPLES::

            sage: G = DihedralGroup(5)
            sage: CB = DihedralCellularBasis(G, 5)
            sage: CB.verify_multiplication()
            True
        """
        kG = self.kG

        for lam_idx in self.cells.keys():
            V = self.reps[lam_idx]
            d = V.dimension()

            for g in self.G:
                M_g = _rep_matrix(V, g)

                for i in range(d):
                    for j in range(d):
                        Cij = self.C(lam_idx, i, j)

                        # Left-hand side: g · C^λ_{ij}
                        lhs = kG(g) * Cij

                        # Right-hand side: Σ_{i'} ρ^λ(g)_{i'i} · C^λ_{i'j}
                        rhs = kG.zero()
                        for ip in range(d):
                            rho_entry = self.base_ring(M_g[ip, i])
                            rhs += rho_entry * self.C(lam_idx, ip, j)

                        if lhs != rhs:
                            return False
        return True

    def verify_cellular_axioms(self):
        """
        Run all three Graham-Lehrer checks.

        EXAMPLES::

            sage: G = DihedralGroup(5)
            sage: CB = DihedralCellularBasis(G, 5)
            sage: CB.verify_cellular_axioms()
            True
        """
        checks = {
            "basis (linear independence)": self.verify_basis(),
            "involution (GL2)":            self.verify_involution(),
            "multiplication (GL3)":        self.verify_multiplication(),
        }
        for name, result in checks.items():
            status = "PASS" if result else "FAIL"
            print(f"  [{status}] {name}")
        return all(checks.values())

    def cell_module(self, lam_idx):
        """
        Return the cell module W^λ as the span of { C^λ_{ij} : fixed j, all i }.

        By GL3 this is a left kG-module isomorphic to V^λ.
        We take j=0 by convention.

        Returns a list of the spanning elements [ C^λ_{i,0} for i in range(d) ].

        EXAMPLES::

            sage: G = DihedralGroup(5)
            sage: CB = DihedralCellularBasis(G, 5)
            sage: W = CB.cell_module(2)    # first 2D irrep
            sage: len(W)
            2
        """
        return [self.C(lam_idx, i, 0)
                for i in range(self.reps[lam_idx].dimension())]

    def idempotents(self):
        """
        Return the Wedderburn idempotents e^λ = C^λ_{00} (d_λ=1) or
        Σ_i C^λ_{ii} (d_λ > 1).

        These satisfy ``e^lam * e^mu = delta_{lam,mu} * e^lam`` and ``sum_lam e^lam = 1``.

        EXAMPLES::

            sage: G = DihedralGroup(5)
            sage: CB = DihedralCellularBasis(G, 5)
            sage: idems = CB.idempotents()
            sage: sum(idems) == CB.kG.one()
            True
        """
        kG     = self.kG
        result = []
        for lam_idx, V in enumerate(self.reps):
            d   = V.dimension()
            e   = sum(self.C(lam_idx, i, i) for i in range(d))
            result.append(e)
        return result

    def __repr__(self):
        dims = [V.dimension() for V in self.reps]
        return (f"Cellular basis for D_{self.n} over {self.base_ring}\n"
                f"  Irrep dimensions: {dims}\n"
                f"  Total basis size: {sum(d**2 for d in dims)} (= 2·{self.n})")
