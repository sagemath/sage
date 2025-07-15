r"""
Library for Batyrev's classification of smooth toric varieties.
TODO: Batyrev's construction for complete toric varieties of Picard rank 3. [Batyrev1991]
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

    # whether to check the input data when constructing a variety, here we by default trust the database of smooth Fano polytopes
    _check = False 

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
            self.__dict__[dict_key] = \
                SmoothFanoToricVariety(Delta=Delta,
                                    coordinate_names=coordinate_names,
                                    base_ring=base_ring,
                                    check=self._check)
        return self.__dict__[dict_key]
       
    def get_smooth_fano(self, dim: int, index: int, coordinate_names=None, base_ring=QQ):
      r"""
      Return the smooth toric Fano ``dim``-fold with the given index.

      The dataset includes all smooth toric Fano 3-folds (124 in total) and 4-folds (124 known examples).
      These varieties are indexed in the Batyrev database using the label format ``F.{dim}D.{index:04d}``.

      INPUT:

      - ``dim`` -- dimension of the Fano toric variety (either 3 or 4)

      - ``index`` -- integer index for the Fano variety (e.g., 0 for F.3D.0000)

      - ``coordinate_names`` -- optional names for homogeneous coordinates

      - ``base_ring`` -- base ring (default: `\QQ`)

      OUTPUT:

      - A smooth Fano toric variety of the specified dimension and index.

      EXAMPLES::

          sage: X = BTF.get_smooth_fano(3, 14)
          sage: X
          3-d smooth Fano toric variety covered by ...
          sage: X.is_smooth()
          True

          sage: Y = BTF.get_smooth_fano(4, 23)
          sage: Y.dimension()
          4
      """
      if dim not in (3, 4):
          raise ValueError("Only dimensions 3 and 4 are supported.")
      if dim == 3 and not (0 <= index < 18):
          raise ValueError("Index for 3D must be in the range [0, 18).")
      if dim == 4 and not (0 <= index < 124):
          raise ValueError("Index for 4D must be in the range [0, 124).")  
      
      label = f"F.{dim}D.{index:04d}"
      return self._make_SmoothFanoToricVariety(label, coordinate_names, base_ring)
      

BTF = BatyrevToricFactory()