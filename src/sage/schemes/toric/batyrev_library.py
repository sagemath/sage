r"""
Library for Batyrev's classification of smooth toric varieties.
TODO: Batyrev's construction for complete toric varieties of Picard rank 3. [Batyrev1991]
TODO: Batyrev's construction for 18 toric Fano 4-folds [Batyrev1999]
TODO: Batyrev's construction for 123 toric Fano 4-folds [Batyrev1999]
TODO: Maybe add the correspoding well-known name for each variety.
"""

from sage.structure.sage_object import SageObject
from sage.geometry.cone import Cone
from sage.geometry.fan import Fan, FaceFan
from sage.geometry.lattice_polytope import LatticePolytope
from sage.geometry.toric_lattice import ToricLattice
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.arith.misc import GCD as gcd
from sage.schemes.toric.variety import (DEFAULT_PREFIX,
                                        ToricVariety,
                                        normalize_names)
from sage.schemes.toric.batyrev import SmoothFanoToricVariety

from .poly_db_3d import polytopes_3d
from .poly_db_4d import polytopes_4d

class BatyrevToricFactory(SageObject):
    r"""
    The methods of this class construct Batyrev' classification of smooth and complete toric varieties.`
    """

    _check = True # AL: what does this do?

    def _make_SmoothFanoToricVariety(self, name, coordinate_names, base_ring):
        r"""
        Construct a toric variety and cache the result.

        INPUT:

        - ``name`` -- string; one of the pre-defined names in the
          ``toric_varieties_rays_cones`` data structure

        - ``coordinate_names`` -- string describing the names of the
          homogeneous coordinates of the toric variety

        - ``base_ring`` -- a ring (default: `\QQ`); the base ring for
          the toric variety

        OUTPUT: a :class:`toric variety
        <sage.schemes.toric.variety.ToricVariety_field>`.

        EXAMPLES::

            sage: toric_varieties.A1()           # indirect doctest
            1-d affine toric variety
        """
        if name in polytopes_3d:
            Delta = LatticePolytope(polytopes_3d[name])
        elif name in polytopes_4d:
            Delta = LatticePolytope(polytopes_4d[name])
        else:
            raise ValueError("unknown family name '%s'" % name)
        
        if coordinate_names is None:
            dict_key = (name, base_ring)
        else:
            coordinate_names = normalize_names(coordinate_names, len(Delta),
                                               DEFAULT_PREFIX)
            dict_key = (name, base_ring) + tuple(coordinate_names)
        if dict_key not in self.__dict__:
            # TODO: to define the points, ray2point, charts etc explicitly
            self.__dict__[dict_key] = \
                SmoothFanoToricVariety(Delta=Delta,
                                    coordinate_names=coordinate_names,
                                    base_ring=base_ring,
                                    check=self._check)
        return self.__dict__[dict_key]
    
    def Fano3fold000(self, coordinate_names=None, base_ring=QQ):
        r"""
        Return the smooth toric Fano 3-fold of type `\mathrm{F}.3\mathrm{D}.0000`.

        INPUT:

        - ``coordinate_names`` -- string describing the names of the
          homogeneous coordinates of the toric variety

        - ``base_ring`` -- a ring (default: `\QQ`); the base ring for
          the toric variety

        OUTPUT: a :class:`toric variety
        <sage.schemes.toric.variety.ToricVariety_field>`.

        EXAMPLES::

            sage: Fano3fold000 = smooth_fano_toric_varieties.Fano3fold000()
            sage: Fano3fold000.is_smooth()
            True

        """
        return self._make_SmoothFanoToricVariety('F.3D.0000', coordinate_names, base_ring)

    def Fano4fold000(self, coordinate_names=None, base_ring=QQ):
        r"""
        Return the smooth toric Fano 4-fold of type `\mathrm{F}.4\mathrm{D}.0000`.

        INPUT:

        - ``coordinate_names`` -- string describing the names of the
          homogeneous coordinates of the toric variety

        - ``base_ring`` -- a ring (default: `\QQ`); the base ring for
          the toric variety

        OUTPUT: a :class:`toric variety
        <sage.schemes.toric.variety.ToricVariety_field>`.

        EXAMPLES::

            sage: Fano4fold000 = smooth_fano_toric_varieties.Fano4fold000()
            sage: Fano4fold000.is_smooth()
            True
        """
        return self._make_SmoothFanoToricVariety('F.4D.0000', coordinate_names, base_ring)
    

smooth_fano_toric_varieties = BatyrevToricFactory()


### AL: change it to the other construction of Batyrev for smooth toric varieties of Picard rank 3 because they use different logic
# def _make_batyrev_matrix(self, P, B, C):
#     p0, p1, p2, p3, p4 = P
#     Cprime = [0] + list(C)
#     R1 = ([1]*p0) + ([1]*p1) + ([-c for c in Cprime]) + ([-(b+1) for b in B]) + ([0]*p4)
#     R2 = ([0]*p0) + ([1]*p1) + ([1]*p2) + ([0]*p3) + ([-1]*p4)
#     R3 = ([0]*p0) + ([0]*p1) + ([1]*p2) + ([1]*p3) + ([0]*p4)
#     R4 = ([0]*p0) + ([-1]*p1) + ([0]*p2) + ([1]*p3) + ([1]*p4)
#     R5 = ([1]*p0) + ([0]*p1) + ([-c for c in Cprime]) + ([-b for b in B]) + ([1]*p4)
#     return matrix(ZZ, [R1, R2, R3, R4, R5])

# def add_family(self, name, P, B, C):
#     """
#     Register a family by its Batyrev parameters P, B, C.
#     P = [p0, p1, p2, p3, p4]
#     B is a list of length p3
#     C is a list of length p2 - 1
#     """
#     A = self._make_batyrev_matrix(P, B, C)
#     K = A.right_kernel()
#     rays = [v for v in K.basis()]
#     p0, p1, p2, p3, p4 = P
#     X0 = list(range(p0))
#     X1 = list(range(p0, p0+p1))
#     X2 = list(range(p0+p1, p0+p1+p2))
#     X3 = list(range(p0+p1+p2, p0+p1+p2+p3))
#     X4 = list(range(p0+p1+p2+p3, p0+p1+p2+p3+p4))
#     X5 = X0 + X4
#     primitive_cols = [X0, X1, X2, X3, X5]
#     d = len(rays[0])
#     cones = []
#     from itertools import combinations
#     for combo in combinations(range(len(rays)), d):
#         if any(set(pc) <= set(combo) for pc in primitive_cols):
#             continue
#         cones.append(list(combo))
#     self.toric_varieties_rays_cones[name] = (rays, cones)

# def make_variety(self, name, P, B, C, coordinate_names=None, base_ring=QQ):
#     """
#     Convenience method: add the family and return its ToricVariety.
#     """
#     self.add_family(name, P, B, C)
#     return self._make_ToricVariety(name, coordinate_names, base_ring)

# # instantiate and register example families
# btf = BatyrevToricFactory()
# # rank-3 Fano 3-folds examples
# btf.add_family('Batyrev_3_26', [2,1,1,1,1], [-2], [])
# btf.add_family('Batyrev_3_29', [1,1,1,1,1], [2], [])
