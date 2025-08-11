r""" "
Smooth Fano toric varieties.

This module provides support for smooth Fano toric varieties, which has been classified by Batyrev in dimensions up to 4.
"""

from sage.geometry.cone import Cone
from sage.geometry.fan import Fan, FaceFan
from sage.geometry.lattice_polytope import LatticePolytope
from sage.geometry.toric_lattice import ToricLattice
from sage.rings.rational_field import QQ
from sage.arith.misc import GCD as gcd
from sage.schemes.toric.variety import (
    DEFAULT_PREFIX,
    ToricVariety_field,
    normalize_names,
)
from sage.schemes.toric.fano_variety import FanoToricVariety_field
from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall
from sage.categories.fields import Fields

_Fields = Fields()


# TODO: add the method to fix a basis of the Picard group, and express divisors in this basis
class SmoothFanoToricVariety_field(FanoToricVariety_field, metaclass=ClasscallMetaclass):
    r"""
    Construct a smooth Fano toric variety from a smooth reflexive polytope.
    
    INPUT:

    - ``Delta`` -- reflexive :class:`lattice polytope
      <sage.geometry.lattice_polytope.LatticePolytopeClass>`. The fan of the
      constructed CPR-Fano toric variety will be a crepant subdivision of the
      *normal fan* of ``Delta``. Either ``Delta`` or ``Delta_polar`` must be
      given, but not both at the same time, since one is completely determined
      by another via :meth:`polar
      <sage.geometry.lattice_polytope.LatticePolytopeClass.polar>` method.

    - ``Delta_polar`` -- reflexive :class:`lattice polytope
      <sage.geometry.lattice_polytope.LatticePolytopeClass>`. The fan of the
      constructed CPR-Fano toric variety will be a crepant subdivision of the
      *face fan* of ``Delta_polar``. Either ``Delta`` or ``Delta_polar`` must
      be given, but not both at the same time, since one is completely
      determined by another via :meth:`polar
      <sage.geometry.lattice_polytope.LatticePolytopeClass.polar>` method.

    - ``coordinate_names`` -- names of variables for the coordinate ring, see
      :func:`~sage.schemes.toric.variety.normalize_names`
      for acceptable formats. If not given, indexed variable names will be
      created automatically.

    - ``names`` -- an alias of ``coordinate_names`` for internal
      use. You may specify either ``names`` or ``coordinate_names``,
      but not both.

    - ``coordinate_name_indices`` -- list of integers, indices for indexed
      variables. If not given, the index of each variable will coincide with
      the index of the corresponding point of ``Delta_polar``.

    - ``base_ring`` -- base field of the Fano toric variety
      (default: `\QQ`)

    - ``base_field`` -- alias for ``base_ring``. Takes precedence if
      both are specified.

    - ``check`` -- by default the input data will be checked for correctness
      (e.g. that ``charts`` do form a subdivision of the normal fan of
      ``Delta``). If you know for sure that the input is valid, you may
      significantly decrease construction time using ``check=False`` option.

    OUTPUT: :class:`smooth Fano toric variety <SmoothFanoToricVariety_field>`

    Examples::

    We start with the product of two projective lines, which is smooth::

        sage: diamond = lattice_polytope.cross_polytope(2)
        sage: diamond.vertices()
        M( 1,  0),        M( 0,  1),
        M(-1,  0),        M( 0, -1)
        in 2-d lattice M
        sage: P1xP1 = SmoothFanoToricVariety(Delta_polar=diamond)
        sage: P1xP1
        2-d smooth Fano toric variety covered by 4 affine patches
        sage: P1xP1.is_smooth()
        True
    """


    @staticmethod
    def __classcall_private__(
        cls,
        Delta=None,
        Delta_polar=None,
        coordinate_names=None,
        names=None,
        coordinate_name_indices=None,
        base_ring=None,
        base_field=None,
        check=True):
      if names is not None:
          if coordinate_names is not None:
              raise ValueError("You must not specify both coordinate_names and names!")
          coordinate_names = names
      if (Delta is None) == (Delta_polar is None):
          raise ValueError("exactly one of Delta and Delta_polar must be given")
      if Delta_polar is None:
          Delta_polar = Delta.polar()

      if check and not Delta_polar.is_reflexive():
          raise ValueError("Delta_polar must be reflexive")

      fan = FaceFan(Delta_polar)

      if check and not fan.is_smooth():
          raise ValueError("The face fan of Delta_polar is not smooth")

      # Check/normalize base_field
      if base_field is not None:
          base_ring = base_field
      if base_ring is None:
          base_ring = QQ
      elif base_ring not in _Fields:
          raise TypeError(
              "need a field to construct a Fano toric variety!" "\n Got %s" % base_ring
        )
      
      return typecall(cls, 
                      Delta_polar=Delta_polar, 
                      fan=fan,
                      coordinate_names=coordinate_names, 
                      coordinate_name_indices=coordinate_name_indices, 
                      base_field=base_ring)

    def __init__(
        self, Delta_polar, fan, coordinate_names, coordinate_name_indices, base_field
    ):
        r"""
        See :class:`SmoothFanoToricVariety_field` for documentation.

        Use ``SmoothFanoToricVariety`` to construct a smooth Fano toric variety.
        """
        FanoToricVariety_field.__init__(self,
            Delta_polar, fan, coordinate_names, coordinate_name_indices, base_field
        )

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT: string

        TESTS::

            sage: P1xP1 = SmoothFanoToricVariety(
            ....:     Delta_polar=lattice_polytope.cross_polytope(2))
            sage: print(P1xP1._repr_())
            2-d smooth Fano toric variety covered by 4 affine patches
        """
        return "%d-d smooth Fano toric variety covered by %d affine patches" % (
            self.dimension_relative(),
            self.fan().ngenerating_cones(),
        )


def SmoothFanoToricVariety(*arg, **kwds):
    return SmoothFanoToricVariety_field(*arg, **kwds)
