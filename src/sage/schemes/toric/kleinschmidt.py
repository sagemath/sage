r"""
Kleinschmidt's classification of smooth complete toric varieties of Picard rank 2.

This module provides a small factory that constructs the toric variety X(d; a)
from the data:
    - d >= 1: ambient lattice dimension
    - a = (a_0, ..., a_{r-1}) with 1 <= r <= d-1  (so s := d - r + 1 >= 2)

Rays (in N = Z^d) are ordered as in the Macaulay2 implementation:
    v_0  = -(e_0 + ... + e_{r-1})
    v_1, ..., v_r            = e_0, ..., e_{r-1}
    v_{r+1}, ..., v_{r+s-1}  = e_r, ..., e_{d-1}   (there are s-1 of these)
    v_{r+s}                  = sum_{i=0}^{r-1} a_i e_i  -  sum_{k=r}^{d-1} e_k

Maximal cones are all complements of a pair (i, j) with
    i in {0,1,...,r} and j in {r+1, ..., r+s}  (indices are in the ray list above).
Equivalently, each maximal cone omits exactly one index from {0,...,r} and
one from {r+1,...,r+s}. The fan is complete and (under the classification
assumptions) smooth of Picard rank 2.

We also provide the 2 x (d+2) "degree matrix"
   deg = [ [  0,    -a_0, ..., -a_{r-1},   1, ..., 1 ],
           [  1,     1,   ...,    1,       0, ..., 0 ] ],
whose columns are aligned with the ordered rays above:
  [v_0 | v_1..v_r | v_{r+1}..v_{r+s-1} | v_{r+s} ].

Sage computes Cl(X) from the fan; its chosen basis may differ by GL(2, Z).
The provided degree matrix is attached to the returned variety as an attribute
for reference:  ``X._kleinschmidt_degree_matrix``

TODO: add references
"""

from sage.structure.sage_object import SageObject
from sage.geometry.fan import Fan
from sage.geometry.toric_lattice import ToricLattice
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.matrix.constructor import matrix
from sage.schemes.toric.variety import (DEFAULT_PREFIX,
                                        ToricVariety,
                                        normalize_names)
from sage.schemes.toric.batyrev import SmoothFanoToricVariety
from sage.schemes.toric.fano_variety import FanoToricVariety

class KleinschidmtFactory(SageObject):
    r"""
    Factory for Kleinschmidt rank-2 toric varieties.

    User entry point:
        get_variety(d: int, a: list[int], coordinate_names=None, base_ring=QQ, check=True) -> SmoothFanoToricVariety

    Examples
    --------
    Hirzebruch surface F_m (d=2, a=[m])::

        sage: X = KTF.get_variety(2, [3])           # F_3
        sage: X.is_complete()
        True
        sage: X.rational_class_group().ngens()     # Picard rank = 2
        2
        sage: X.fan().nrays()
        4
        sage: X.is_smooth()
        True
        sage: X._kleinschmidt_degree_matrix         # 2 x 4 matrix aligned with rays
        [ 0 -3  1  1]
        [ 1  1  0  0]

    A 3-fold example (d=3, a=[2])::

        sage: X = KTF.get_variety(3, [2])
        sage: (X.dimension(), X.rational_class_group().ngens(), X.fan().nrays())
        (3, 2, 5)
    """

    def _validate_input(self, d, a):
        if d < 1:
            raise ValueError("d must be a positive integer (d >= 1).")
        if not isinstance(a, (list, tuple)):
            raise TypeError("a must be a list or tuple of integers.")
        r = len(a)
        if r < 1:
            raise ValueError("the list a must be nonempty (r >= 1).")
        if r >= d:
            raise ValueError("list a is too long: require len(a) <= d-1.")
        # s = d - r + 1 >= 2 automatically follows from r <= d-1
        return r, d - r + 1

    @staticmethod
    def _degree_matrix(d, a):
        """
        Return the 2 x (d+2) degree matrix aligned with the ordered rays:
        [v0 | v1..vr | v_{r+1}..v_{r+s-1} | v_{r+s}],
        with r = len(a) and s = d - r + 1.

        Columns are:
        (0,1), (-a_i,1) for i=0..r-1, and s copies of (1,0).
        """
        r = len(a)
        s = int(d) - r + 1
        if s <= 0:
            raise ValueError("Require len(a) <= d-1 so that s = d - r + 1 >= 2.")

        row1 = [ZZ(0)] + [-ZZ(ai) for ai in a] + [ZZ(1)] * s
        row2 = [ZZ(1)] + [ZZ(1)] * r          + [ZZ(0)] * s

        # sanity: total columns must be d+2
        ncols = 1 + r + s
        assert ncols == int(d) + 2, f"degree matrix should have {int(d)+2} columns, got {ncols}"

        return matrix(ZZ, [row1, row2])


    def get_variety(self, d, a, coordinate_names=None, base_ring=QQ, check=True):
        r"""
        Construct the Kleinschmidt rank-2 toric variety X(d; a).

        INPUT:

        - ``d`` -- integer (`d \geq 1`); the ambient lattice dimension.

        - ``a`` -- list (or tuple) of integers `a = (a_0,\dots,a_{r-1})` with
        `1 \le r \le d-1`; determines the family parameter in Kleinschmidt's
        classification. We set `s := d - r + 1`.

        - ``coordinate_names`` -- either a string prefix (e.g. `'x'` yields
        `x0, x1, ...`), a list of names of length ``d+2`` (the number of rays),
        or ``None`` (default: automatically generated names).
        See :func:`~sage.schemes.toric.variety.normalize_names`.

        - ``base_ring`` -- a field (default: :math:`\QQ`); the base ring for the
        toric variety. Must lie in :class:`~sage.categories.fields.Fields`.

        - ``check`` -- boolean (default: ``True``); whether to validate the fan data.

        OUTPUT:

        A :class:`toric variety
        <sage.schemes.toric.variety.ToricVariety_field>`.

        The attribute ``X._kleinschmidt_degree_matrix`` stores the degree matrix aligned with the ordered rays used here

        EXAMPLES::

            sage: from sage.schemes.toric.kleinschmidt import KleinschidmtFactory
            sage: KTF = KleinschidmtFactory()

        A Hirzebruch surface (Fano case ``F_1``)::

            sage: X = KTF.get_variety(2, [1])
            sage: X.dimension(), X.rational_class_group().ngens(), X.fan().nrays()
            (2, 2, 4)
            sage: X._kleinschmidt_degree_matrix.nrows(), X._kleinschmidt_degree_matrix.ncols()
            (2, 4)

        A 3-fold example::

            sage: Y = KTF.get_variety(3, [2])
            sage: Y.dimension(), Y.rational_class_group().ngens(), Y.fan().nrays()
            (3, 2, 5)

        """

        r, s = self._validate_input(d, a)

        # Lattice and standard basis
        N = ToricLattice(d)
        e = N.basis()  # e[0], ..., e[d-1]

        # Rays (ordered exactly as in the original paper)
        rays = []
        v0 = -(sum(e[i] for i in range(r)))                  # v_0
        rays.append(v0)
        rays.extend(e[i] for i in range(r))                  # v_1..v_r
        rays.extend(e[r + j] for j in range(s - 1))          # v_{r+1}..v_{r+s-1}
        last = sum(a[i] * e[i] for i in range(r)) - sum(e[r + j] for j in range(s - 1))
        rays.append(last)                                     # v_{r+s}
        assert len(rays) == d + 2

        # Cones: for each i in 0..r and j in r+1..r+s, take the complement of {i, j}
        L = list(range(0, r + s + 1))  # 0..(r+s) inclusive = d+1 inclusive
        cones = []
        for i in range(0, r + 1):
            for j in range(r + 1, r + s + 1):
                cone = tuple(k for k in L if k != i and k != j)
                cones.append(cone)

        fan = Fan(cones, rays=rays, check=check)

        # Coordinate names
        num_rays = d + 2
        if coordinate_names is not None:
            coordinate_names = normalize_names(coordinate_names, num_rays, DEFAULT_PREFIX)
        # AL: So far I cannot find a way to check if a fan is Fano cheaply. So we always return a ToricVariety.
        # instead of a SmoothFanoToricVariety.
        # if fan.is_smooth() and is_fano_fan(fan):
        #     X = SmoothFanoToricVariety(fan=fan, coordinate_names=coordinate_names, base_ring=base_ring, check=check)
        # else:
        X = ToricVariety(fan=fan, coordinate_names=coordinate_names, base_ring=base_ring)

        # attach the 2 x (d+2) degree matrix for reference TODO: maybe delete?
        X._kleinschmidt_degree_matrix = self._degree_matrix(d, a)
        return X

KTF = KleinschidmtFactory()

# def is_fano_fan(fan):
#     """
#     Helper function to return True if the complete fan defines a Fano toric variety.
    
#     TESTS::
#         sage: from sage.schemes.toric.kleinschmidt import is_fano_fan
#         sage: N = ToricLattice(2)
#         sage: e1, e2 = N.basis()
#         sage: rays = [e1, e2, -e1-e2]
#         sage: cones = [(0,1), (1,2), (0,2)]
#         sage: fan = Fan(cones, rays=rays)
#         sage: is_fano_fan(fan)
#         True
#         sage: rays = [e2, e1, -e2, -e1 + 3*e2] # F_3 is not Fano
#         sage: fan = Fan(cones, rays=rays)
#         sage: is_fano_fan(fan)
#         False

#     """
#     if not fan.is_complete():
#         return False

#     N = fan.lattice()
#     # Ensure primitive generators
#     rays = [N(v) for v in fan.rays()]

#     # Build inequalities: <m, v> + 1 >= 0
#     ieqs = [[QQ(1)] + [QQ(c) for c in v] for v in rays]
#     P = Polyhedron(ieqs=ieqs)

#     for v in P.vertices():
#         vec = v.vector()
#         if any(c.denominator() != 1 for c in vec):
#             return False

#     # Fano iff bounded: no rays and no lines in recession cone
#     return P.n_rays() == 0 and P.n_lines() == 0
