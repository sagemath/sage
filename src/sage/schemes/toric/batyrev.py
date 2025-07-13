
from sage.structure.sage_object import SageObject

from sage.matrix.constructor import matrix
from sage.matrix.special import identity_matrix
from sage.geometry.cone import Cone
from sage.geometry.fan import Fan, FaceFan
from sage.geometry.lattice_polytope import LatticePolytope
from sage.geometry.toric_lattice import ToricLattice
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.arith.misc import GCD as gcd
from sage.schemes.toric.variety import (DEFAULT_PREFIX,
                                        ToricVariety_field,
                                        normalize_names)

from sage.categories.fields import Fields
_Fields = Fields()



# TODO: add the method to fix a basis of the Picard group, and express divisors in this basis
class SmoothFanoToricVariety_field(ToricVariety_field):
    r"""
    Construct a Fano smooth toric variety from a reflexive polytope.

    INPUT:
    - ``P`` -- reflexive polytope

    - ``fan`` -- rational polyhedral fan subdividing the face fan of
      ``Delta_polar``

    - ``coordinate_points`` -- list of indices of points of ``Delta_polar``
      used for rays of ``fan``

    - ``coordinate_names`` -- names of the variables of the coordinate ring in
      the format accepted by
      :func:`~sage.schemes.toric.variety.normalize_names`

    - ``coordinate_name_indices`` -- indices for indexed variables,
      if ``None``, will be equal to ``coordinate_points``

    - ``base_field`` -- base field of the smooth toric variety

    OUTPUT: :class:`smooth Fano toric variety <SmoothFanoToricVariety_field>`
    """

    def __init__(self, Delta_polar, fan, coordinate_points, point_to_ray,
                 coordinate_names, coordinate_name_indices, base_field):
        r"""
        See :class:`SmoothFanoToricVariety_field` for documentation.

        Use ``SmoothFanoToricVariety`` to construct a smooth Fano toric variety.

        TESTS:

            sage: 
        """
        self._Delta_polar = Delta_polar
        self._coordinate_points = tuple(coordinate_points)
        self._point_to_ray = point_to_ray # AL: I may want to drop this later, 
        if coordinate_name_indices is None:
            coordinate_name_indices = coordinate_points
        super().__init__(fan, coordinate_names, 
                         coordinate_name_indices, base_field)
        
        # TODO: the _latex_ and _repr_ methods


    def change_ring(self, F):
        r"""
        Return a smooth Fano toric variety over field ``F``, otherwise the same
        as ``self``.

        INPUT:

        - ``F`` -- field

        OUTPUT: :class:`smooth Fano toric variety <SmoothFanoToricVariety_field>` over ``F``

        .. NOTE::

            There is no need to have any relation between ``F`` and the base
            field of ``self``. If you do want to have such a relation, use
            :meth:`base_extend` instead.

        EXAMPLES::

        """
        if self.base_ring() == F:
            return self
        elif F not in _Fields:
            raise TypeError("need a field to construct a Fano toric variety!"
                            "\n Got %s" % F)
        else:
            return SmoothFanoToricVariety_field(self._P, self._fan,
                self._coordinate_points, self._point_to_ray,
                self.variable_names(), None, F)
                # coordinate_name_indices do not matter, we give explicit
                # names for all variables

    def coordinate_point_to_coordinate(self, point):
        # TODO: if this is necessary
        pass

    def coordinate_points(self):
        r"""
        Return indices of points of :meth:`Delta_polar` used for coordinates.

        OUTPUT: :class:`tuple` of integers
        """
        return self._coordinate_points

    def Delta(self):
        r"""
        Return the reflexive polytope <sage.geometry.lattice_polytope.LatticePolytopeClass>` 
        defining ``self``.

        OUTPUT: 
        
        - reflexive :class:`lattice polytope`

        EXAMPLES::

        """
        return self._Delta_polar.polar()
    
    def Delta_polar(self):
        r"""
        Return polar of :meth:`Delta`.

        OUTPUT: 
        
        - reflexive :class:`lattice polytope`

        EXAMPLES::

        """
        return self._Delta_polar

def SmoothFanoToricVariety(Delta=None,
                        Delta_polar=None,
                        coordinate_points=None,
                        charts=None,
                        coordinate_names=None,
                        names=None,
                        coordinate_name_indices=None,
                        make_simplicial=False,
                        base_ring=None,
                        base_field=None,
                        check=True):
    r"""
    Construct a CPR-Fano toric variety.
    """
    if (Delta is None) == (Delta_polar is None):
        raise ValueError("exactly one of Delta and Delta_polar must be given")
    if Delta_polar is None:
        Delta_polar = Delta.polar()

    # AL: maybe not checking reflexivity here, since we already have a database of smooth reflexive polytopes
    #     unless the reflexivity check is not expensive here, we can add it for future more general use
    # if not Delta_polar.is_reflexive():
    #     raise ValueError("Delta_polar must be reflexive")

    # Check/normalize coordinate_points and construct fan rays
    if coordinate_points is None:
        coordinate_points = list(range(Delta_polar.nvertices()))
        if charts is not None:
            for chart in charts:
                for point in chart:
                    if point not in coordinate_points:
                        coordinate_points.append(point)
    elif coordinate_points == "vertices":
        coordinate_points = list(range(Delta_polar.nvertices()))
    elif coordinate_points == "all":
        coordinate_points = list(range(Delta_polar.npoints()))
        coordinate_points.remove(Delta_polar.origin())
    elif coordinate_points == "all but facets":
        coordinate_points = Delta_polar.skeleton_points(Delta_polar.dim() - 2)
    elif isinstance(coordinate_points, str):
        raise ValueError("unrecognized description of the coordinate points!"
                         "\nGot: %s" % coordinate_points)
    elif check:
        cp_set = set(coordinate_points)
        if len(cp_set) != len(coordinate_points):
            raise ValueError(
                "no repetitions are allowed for coordinate points!\nGot: %s"
                % coordinate_points)
        if not cp_set.issuperset(list(range(Delta_polar.nvertices()))):
            raise ValueError("all %d vertices of Delta_polar must be used "
                "for coordinates!\nGot: %s"
                % (Delta_polar.nvertices(), coordinate_points))
        if Delta_polar.origin() in cp_set:
            raise ValueError("the origin (point #%d) cannot be used for a "
                "coordinate!\nGot: %s"
                % (Delta_polar.origin(), coordinate_points))
    point_to_ray = {point: n
                    for n, point in enumerate(coordinate_points)}
    
    # This can be simplified if LatticePolytopeClass is adjusted.
    rays = [Delta_polar.point(p) for p in coordinate_points]
    # Check/normalize charts and construct the fan based on them.
    if charts is None:
        # Start with the face fan
        fan = FaceFan(Delta_polar)
    else:
        # First of all, check that each chart is completely contained in a
        # single facet of Delta_polar, otherwise they do not form a
        # subdivision of the face fan of Delta_polar
        if check:
            facet_sets = [frozenset(facet.ambient_point_indices())
                          for facet in Delta_polar.facets()]
            for chart in charts:
                is_bad = True
                for fset in facet_sets:
                    if fset.issuperset(chart):
                        is_bad = False
                        break
                if is_bad:
                    raise ValueError(
                        "%s does not form a chart of a subdivision of the "
                        "face fan of %s!" % (chart, Delta_polar))
        # AL: maybe skip this logic
        cones = [Cone((rays[point_to_ray[point]] for point in chart),
                      check=check)
                 for chart in charts]
        fan = Fan(cones, check=check)
        if check and not fan.is_complete():
            raise ValueError("given charts do not form a complete fan!")

    return SmoothFanoToricVariety_field(
        Delta_polar, fan, coordinate_points,
        point_to_ray, coordinate_names, coordinate_name_indices, base_ring)

