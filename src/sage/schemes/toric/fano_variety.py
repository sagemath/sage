# sage.doctest: needs sage.geometry.polyhedron sage.graphs
r"""
Fano toric varieties

This module provides base class of (Gorenstein) Fano varieties, so we can naturally
introduce the subclass of CPR-Fano toric variety, as well as a factory 
for smooth Fano toric varieties.
"""

import re

from sage.geometry.cone import Cone
from sage.geometry.fan import FaceFan
from sage.geometry.fan import Fan
from sage.geometry.lattice_polytope import LatticePolytope
from sage.misc.latex import latex
from sage.misc.misc_c import prod
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ

from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_base
from sage.rings.polynomial.polynomial_ring import PolynomialRing_generic
from sage.rings.fraction_field import FractionField_generic

from sage.schemes.toric.toric_subscheme import AlgebraicScheme_subscheme_toric
from sage.schemes.toric.variety import (
                                            ToricVariety_field,
                                            normalize_names)
from sage.structure.all import coercion_model
from sage.categories.fields import Fields
_Fields = Fields()


# Default coefficient for anticanonical hypersurfaces
DEFAULT_COEFFICIENT = "a"
# Default coefficients for nef complete intersections
DEFAULT_COEFFICIENTS = tuple(chr(i) for i in range(ord("a"), ord("z") + 1))



def FanoToricVariety(Delta=None,
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
    Construct a Fano toric variety.

    .. NOTE::

        See documentation of the module
        :mod:`~sage.schemes.toric.fano_variety` for the used
        definitions and supported varieties.

    Due to the large number of available options, it is recommended to always
    use keyword parameters.

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

    - ``coordinate_points`` -- list of integers or string. A list will be
      interpreted as indices of (boundary) points of ``Delta_polar`` which
      should be used as rays of the underlying fan. It must include all
      vertices of ``Delta_polar`` and no repetitions are allowed. A string
      must be one of the following descriptions of points of ``Delta_polar``:

      * "vertices" (default),
      * "all" (will not include the origin),
      * "all but facets" (will not include points in the relative interior of
        facets);

    - ``charts`` -- list of lists of elements from ``coordinate_points``. Each
      of these lists must define a generating cone of a fan subdividing the
      normal fan of ``Delta``. Default ``charts`` correspond to the normal fan
      of ``Delta`` without subdivision. The fan specified by ``charts`` will
      be subdivided to include all of the requested ``coordinate_points``.

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

    - ``make_simplicial`` -- if ``True``, the underlying fan will be made
      simplicial (default: ``False``)

    - ``base_ring`` -- base field of the CPR-Fano toric variety
      (default: `\QQ`)

    - ``base_field`` -- alias for ``base_ring``. Takes precedence if
      both are specified.

    - ``check`` -- by default the input data will be checked for correctness
      (e.g. that ``charts`` do form a subdivision of the normal fan of
      ``Delta``). If you know for sure that the input is valid, you may
      significantly decrease construction time using ``check=False`` option.

    OUTPUT: :class:`CPR-Fano toric variety <CPRFanoToricVariety_field>`

    EXAMPLES:

    We start with the product of two projective lines::

        sage: diamond = lattice_polytope.cross_polytope(2)
        sage: diamond.vertices()
        M( 1,  0),        M( 0,  1),
        M(-1,  0),        M( 0, -1)
        in 2-d lattice M
        sage: P1xP1 = CPRFanoToricVariety(Delta_polar=diamond)
        sage: P1xP1
        2-d CPR-Fano toric variety covered by 4 affine patches
        sage: P1xP1.fan()
        Rational polyhedral fan in 2-d lattice M
        sage: P1xP1.fan().rays()
        M( 1,  0),        M( 0,  1),
        M(-1,  0),        M( 0, -1)
        in 2-d lattice M

    "Unfortunately," this variety is smooth to start with and we cannot
    perform any subdivisions of the underlying fan without leaving the
    category of CPR-Fano toric varieties. Our next example starts with a
    square::

        sage: square = diamond.polar()
        sage: square.vertices()
        N( 1,  1),        N( 1, -1),
        N(-1, -1),        N(-1,  1)
        in 2-d lattice N
        sage: square.points()
        N( 1,  1),        N( 1, -1),        N(-1, -1),
        N(-1,  1),        N(-1,  0),        N( 0, -1),
        N( 0,  0),        N( 0,  1),        N( 1,  0)
        in 2-d lattice N

    We will construct several varieties associated to it::

        sage: FTV = CPRFanoToricVariety(Delta_polar=square)
        sage: FTV.fan().rays()
        N( 1,  1),        N( 1, -1),
        N(-1, -1),        N(-1,  1)
        in 2-d lattice N
        sage: FTV.gens()
        (z0, z1, z2, z3)

        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ....:                           coordinate_points=[0,1,2,3,8])
        sage: FTV.fan().rays()
        N( 1,  1),        N( 1, -1),        N(-1, -1),
        N(-1,  1),        N( 1,  0)
        in 2-d lattice N
        sage: FTV.gens()
        (z0, z1, z2, z3, z8)

        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ....:                           coordinate_points=[8,0,2,1,3],
        ....:                           coordinate_names='x+')
        sage: FTV.fan().rays()
        N( 1,  0),        N( 1,  1),        N(-1, -1),
        N( 1, -1),        N(-1,  1)
        in 2-d lattice N
        sage: FTV.gens()
        (x8, x0, x2, x1, x3)

        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ....:                           coordinate_points='all',
        ....:                           coordinate_names="x y Z+")
        sage: FTV.fan().rays()
        N( 1,  1),        N( 1, -1),        N(-1, -1),        N(-1,  1),
        N(-1,  0),        N( 0, -1),        N( 0,  1),        N( 1,  0)
        in 2-d lattice N
        sage: FTV.gens()
        (x, y, Z2, Z3, Z4, Z5, Z7, Z8)

    Note that ``Z6`` is "missing". This is due to the fact that the 6-th point
    of ``square`` is the origin, and all automatically created names have the
    same indices as corresponding points of
    :meth:`~CPRFanoToricVariety_field.Delta_polar`. This is usually very
    convenient, especially if you have to work with several partial
    resolutions of the same Fano toric variety. However, you can change it, if
    you want::

        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ....:                           coordinate_points='all',
        ....:                           coordinate_names="x y Z+",
        ....:                           coordinate_name_indices=list(range(8)))
        sage: FTV.gens()
        (x, y, Z2, Z3, Z4, Z5, Z6, Z7)

    Note that you have to provide indices for *all* variables, including those
    that have "completely custom" names. Again, this is usually convenient,
    because you can add or remove "custom" variables without disturbing too
    much "automatic" ones::

        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ....:                           coordinate_points='all',
        ....:                           coordinate_names="x Z+",
        ....:                           coordinate_name_indices=list(range(8)))
        sage: FTV.gens()
        (x, Z1, Z2, Z3, Z4, Z5, Z6, Z7)

    If you prefer to always start from zero, you will have to shift indices
    accordingly::

        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ....:                           coordinate_points='all',
        ....:                           coordinate_names="x Z+",
        ....:                           coordinate_name_indices=[0] + list(range(7)))
        sage: FTV.gens()
        (x, Z0, Z1, Z2, Z3, Z4, Z5, Z6)

        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ....:                           coordinate_points='all',
        ....:                           coordinate_names="x y Z+",
        ....:                           coordinate_name_indices=[0]*2 + list(range(6)))
        sage: FTV.gens()
        (x, y, Z0, Z1, Z2, Z3, Z4, Z5)

    So you always can get any names you want, somewhat complicated default
    behaviour was designed with the hope that in most cases you will have no
    desire to provide different names.

    Now we will use the possibility to specify initial charts::

        sage: charts = [(0,1), (1,2), (2,3), (3,0)]

    (these charts actually form exactly the face fan of our square) ::

        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ....:                           coordinate_points=[0,1,2,3,4],
        ....:                           charts=charts)
        sage: FTV.fan().rays()
        N( 1,  1),        N( 1, -1),        N(-1, -1),
        N(-1,  1),        N(-1,  0)
        in 2-d lattice N
        sage: [cone.ambient_ray_indices() for cone in FTV.fan()]
        [(0, 1), (1, 2), (2, 4), (3, 4), (0, 3)]

    If charts are wrong, it should be detected::

        sage: bad_charts = charts + [(3,0)]
        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ....:                           coordinate_points=[0,1,2,3,4],
        ....:                           charts=bad_charts)
        Traceback (most recent call last):
        ...
        ValueError: you have provided 5 cones, but only 4 of them are maximal!
        Use discard_faces=True if you indeed need to construct a fan from these cones.

    These charts are technically correct, they just happened to list one of
    them twice, but it is assumed that such a situation will not happen. It is
    especially important when you try to speed up your code::

        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ....:                           coordinate_points=[0,1,2,3,4],
        ....:                           charts=bad_charts,
        ....:                           check=False)
        Traceback (most recent call last):
        ...
        IndexError: list assignment index out of range

    In this case you still get an error message, but it is harder to figure out
    what is going on. It may also happen that "everything will still work" in
    the sense of not crashing, but work with such an invalid variety may lead to
    mathematically wrong results, so use ``check=False`` carefully!

    Here are some other possible mistakes::

        sage: bad_charts = charts + [(0,2)]
        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ....:                           coordinate_points=[0,1,2,3,4],
        ....:                           charts=bad_charts)
        Traceback (most recent call last):
        ...
        ValueError: (0, 2) does not form a chart of a subdivision of
        the face fan of 2-d reflexive polytope #14 in 2-d lattice N!

        sage: bad_charts = charts[:-1]
        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ....:                           coordinate_points=[0,1,2,3,4],
        ....:                           charts=bad_charts)
        Traceback (most recent call last):
        ...
        ValueError: given charts do not form a complete fan!

        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ....:                           coordinate_points=[1,2,3,4])
        Traceback (most recent call last):
        ...
        ValueError: all 4 vertices of Delta_polar must be used for coordinates!
        Got: [1, 2, 3, 4]

        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ....:                           coordinate_points=[0,0,1,2,3,4])
        Traceback (most recent call last):
        ...
        ValueError: no repetitions are allowed for coordinate points!
        Got: [0, 0, 1, 2, 3, 4]

        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ....:                           coordinate_points=[0,1,2,3,6])
        Traceback (most recent call last):
        ...
        ValueError: the origin (point #6) cannot be used for a coordinate!
        Got: [0, 1, 2, 3, 6]

    Here is a shorthand for defining the toric variety and homogeneous
    coordinates in one go::

        sage: P1xP1.<a,b,c,d> = CPRFanoToricVariety(Delta_polar=diamond)
        sage: (a^2+b^2) * (c+d)
        a^2*c + b^2*c + a^2*d + b^2*d
    """
    if names is not None:
        if coordinate_names is not None:
            raise ValueError('You must not specify both coordinate_names and names!')
        coordinate_names = names
    # Check/normalize Delta_polar
    if Delta is None and Delta_polar is None:
        raise ValueError("either Delta or Delta_polar must be given!")
    elif Delta is not None and Delta_polar is not None:
        raise ValueError("Delta and Delta_polar cannot be given together!")
    elif Delta_polar is None:
        Delta_polar = Delta.polar()
    elif not Delta_polar.is_reflexive():
        raise ValueError("Delta_polar must be reflexive!")
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
        # We will construct the initial fan from Cone objects: since charts
        # may not use all of the necessary rays, alternative form is tedious
        # With check=False it should not be long anyway.
        cones = [Cone((rays[point_to_ray[point]] for point in chart),
                      check=check)
                 for chart in charts]
        fan = Fan(cones, check=check)
        if check and not fan.is_complete():
            raise ValueError("given charts do not form a complete fan!")
    # Subdivide this fan to use all required points
    fan = fan.subdivide(new_rays=(ray for ray in rays
                                      if ray not in fan.rays().set()),
                        make_simplicial=make_simplicial)
    # Now create yet another fan making sure that the order of the rays is
    # the same as requested (it is a bit difficult to get it from the start)
    trans = {}
    for n, ray in enumerate(fan.rays()):
        trans[n] = rays.index(ray)
    cones = tuple(tuple(sorted(trans[r] for r in cone.ambient_ray_indices()))
                  for cone in fan)
    fan = Fan(cones, rays, check=False)
    # Check/normalize base_field
    if base_field is not None:
        base_ring = base_field
    if base_ring is None:
        base_ring = QQ
    elif base_ring not in _Fields:
        raise TypeError("need a field to construct a Fano toric variety!"
                        "\n Got %s" % base_ring)
    fan._is_complete = True     # At this point it must be for sure
    return CPRFanoToricVariety_field(
        Delta_polar, fan, coordinate_points,
        point_to_ray, coordinate_names, coordinate_name_indices, base_ring)


class FanoToricVariety_field(ToricVariety_field):
    r"""
    Base class for Gorenstein Fano toric varieties over a field.
    Construct a Fano toric variety associated to a reflexive polytope.

    .. WARNING::

        This class does not perform any checks of correctness of input and it
        does assume that the internal structure of the given parameters is
        coordinated in a certain way. Use
        :func:`FanoToricVariety` to construct Fano toric varieties.

    .. NOTE::

        See documentation of the module
        :mod:`~sage.schemes.toric.fano_variety` for the used
        definitions and supported varieties.

    INPUT:

    - ``Delta_polar`` -- reflexive polytope

    - ``fan`` -- rational polyhedral fan which is the face fan of
      ``Delta_polar``

    - ``coordinate_names`` -- names of the variables of the coordinate ring in
      the format accepted by
      :func:`~sage.schemes.toric.variety.normalize_names`

    - ``coordinate_name_indices`` -- indices for indexed variables,
      if ``None``, will be equal to ``coordinate_points``

    - ``base_field`` -- base field of the Fano toric variety

    OUTPUT: :class:`Fano toric variety <FanoToricVariety_field>`

    TESTS::

        sage: P1xP1 = FanoToricVariety(
        ....:     Delta_polar=lattice_polytope.cross_polytope(2))
        sage: P1xP1
        2-d Fano toric variety covered by 4 affine patches
    """

    def __init__(self, Delta_polar, fan, coordinate_names, coordinate_name_indices, base_field):
        r"""
        See :class:`FanoToricVariety_field` for documentation.

        Use ``FanoToricVariety`` to construct Fano toric varieties.

        TESTS::

            sage: P1xP1 = FanoToricVariety(
            ....:     Delta_polar=lattice_polytope.cross_polytope(2))
            sage: P1xP1
            2-d CPR-Fano toric variety covered by 4 affine patches
        """
        self._Delta_polar = Delta_polar
        super().__init__(fan, coordinate_names,
                         coordinate_name_indices, base_field)
    
    # AL: we want to use the anticanonical surface and nef complete intersection functions
    @property
    def _coordinate_points(self):
        return list(range(self.fan().rays()))
    
    @property
    def _point_to_ray(self):
        return {i: i for i in self._coordinate_points}

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        OUTPUT: string

        TESTS::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: print(P1xP1._latex_())
            \mathbb{P}_{\Delta^{2}_{14}}
        """
        return r"\mathbb{P}_{%s}" % latex(self.Delta())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT: string

        TESTS::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: print(P1xP1._repr_())
            2-d CPR-Fano toric variety covered by 4 affine patches
        """
        return ("%d-d Fano toric variety covered by %d affine patches"
                % (self.dimension_relative(), self.fan().ngenerating_cones()))

    def anticanonical_hypersurface(self, **kwds):
        r"""
        Return an anticanonical hypersurface of ``self``.

        .. NOTE::

            The returned hypersurface may be actually a subscheme of
            **another** Fano toric variety: if the base field of ``self``
            does not include all of the required names for generic monomial
            coefficients, it will be automatically extended.

        Below `\Delta` is the reflexive polytope corresponding to ``self``,
        i.e. the fan of ``self`` is the normal fan of
        `\Delta`. This function accepts only keyword parameters.

        INPUT:

        - ``monomial_points`` -- list of integers or a string. A list will be
          interpreted as indices of points of `\Delta` which should be used
          for monomials of this hypersurface. A string must be one of the
          following descriptions of points of `\Delta`:

          * "vertices",
          * "vertices+origin",
          * "all",
          * "simplified" (default) -- all points of `\Delta` except for
            the interior points of facets, this choice corresponds to working
            with the "simplified polynomial moduli space" of anticanonical
            hypersurfaces;

        - ``coefficient_names`` -- names for the monomial coefficients, see
          :func:`~sage.schemes.toric.variety.normalize_names`
          for acceptable formats. If not given, indexed coefficient names will
          be created automatically.

        - ``coefficient_name_indices`` -- list of integers, indices for
          indexed coefficients. If not given, the index of each coefficient
          will coincide with the index of the corresponding point of `\Delta`.

        - ``coefficients`` -- as an alternative to specifying coefficient
          names and/or indices, you can give the coefficients themselves as
          arbitrary expressions and/or strings. Using strings allows you to
          easily add "parameters": the base field of ``self`` will be extended
          to include all necessary names.

        OUTPUT:

        - an :class:`anticanonical hypersurface <AnticanonicalHypersurface>` of
          ``self`` (with the extended base field, if necessary).

        EXAMPLES:

        We realize the projective plane as a Fano toric variety::

            sage: simplex = LatticePolytope([(1,0), (0,1), (-1,-1)])
            sage: P2 = FanoToricVariety(Delta_polar=simplex)

        Its anticanonical "hypersurface" is a one-dimensional Calabi-Yau
        manifold::

            sage: P2.anticanonical_hypersurface(monomial_points='all')
            Closed subscheme of 2-d Fano toric variety
             covered by 3 affine patches defined by:
              a0*z0^3 + a9*z0^2*z1 + a7*z0*z1^2 + a1*z1^3 + a8*z0^2*z2 + a6*z0*z1*z2
              + a4*z1^2*z2 + a5*z0*z2^2 + a3*z1*z2^2 + a2*z2^3
        """
        #  AL: we include this function because the construction only rely on the Delta and cox variables 
        return AnticanonicalHypersurface(self, **kwds)

    def change_ring(self, F):
        r"""
        Return a Fano toric variety over field ``F``, otherwise the same
        as ``self``.

        INPUT:

        - ``F`` -- field

        OUTPUT: :class:`Fano toric variety <FanoToricVariety_field>` over ``F``

        .. NOTE::

            There is no need to have any relation between ``F`` and the base
            field of ``self``. If you do want to have such a relation, use
            :meth:`base_extend` instead.

        EXAMPLES::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1xP1.base_ring()
            Rational Field
            sage: P1xP1_RR = P1xP1.change_ring(RR)
            sage: P1xP1_RR.base_ring()
            Real Field with 53 bits of precision
            sage: P1xP1_QQ = P1xP1_RR.change_ring(QQ)
            sage: P1xP1_QQ.base_ring()
            Rational Field
            sage: P1xP1_RR.base_extend(QQ)
            Traceback (most recent call last):
            ...
            ValueError: no natural map from the base ring
            (=Real Field with 53 bits of precision) to R (=Rational Field)!
            sage: R = PolynomialRing(QQ, 2, 'a')
            sage: P1xP1.change_ring(R)
            Traceback (most recent call last):
            ...
            TypeError: need a field to construct a Fano toric variety!
            Got Multivariate Polynomial Ring in a0, a1 over Rational Field
        """
        if self.base_ring() == F:
            return self
        elif F not in _Fields:
            raise TypeError("need a field to construct a Fano toric variety!"
                            "\n Got %s" % F)
        else:
            return FanoToricVariety_field(self._Delta_polar, self._fan,
                self.variable_names(), None, F)
        
    def Delta(self):
        r"""
        Return the reflexive polytope associated to ``self``.

        OUTPUT:

        - reflexive :class:`lattice polytope
          <sage.geometry.lattice_polytope.LatticePolytopeClass>`. The
          underlying fan of ``self`` is the
          *normal fan* of this polytope.

        EXAMPLES::

            sage: diamond = lattice_polytope.cross_polytope(2)
            sage: P1xP1 = FanoToricVariety(Delta_polar=diamond)
            sage: P1xP1.Delta()
            2-d reflexive polytope #14 in 2-d lattice N
            sage: P1xP1.Delta() is diamond.polar()
            True
        """
        return self._Delta_polar.polar()

    def Delta_polar(self):
        r"""
        Return polar of :meth:`Delta`.

        OUTPUT:

        - reflexive :class:`lattice polytope
          <sage.geometry.lattice_polytope.LatticePolytopeClass>`. The
          underlying fan of ``self`` is the *face fan* of this polytope.

        EXAMPLES::

            sage: diamond = lattice_polytope.cross_polytope(2)
            sage: P1xP1 = FanoToricVariety(Delta_polar=diamond)
            sage: P1xP1.Delta_polar()
            2-d reflexive polytope #3 in 2-d lattice M
            sage: P1xP1.Delta_polar() is diamond
            True
            sage: P1xP1.Delta_polar() is P1xP1.Delta().polar()
            True
        """
        return self._Delta_polar

    def nef_complete_intersection(self, nef_partition, **kwds):
        r"""
        Return a nef complete intersection in ``self``.

        .. NOTE::

            The returned complete intersection may be actually a subscheme of
            **another** Fano toric variety: if the base field of ``self``
            does not include all of the required names for monomial
            coefficients, it will be automatically extended.

        Below `\Delta` is the reflexive polytope corresponding to ``self``,
        i.e. the fan of ``self`` is the normal fan of
        `\Delta`. Other polytopes are described in the documentation of
        :class:`nef-partitions <sage.geometry.lattice_polytope.NefPartition>`
        of :class:`reflexive polytopes
        <sage.geometry.lattice_polytope.LatticePolytopeClass>`.

        Except for the first argument, ``nef_partition``, this method accepts
        only keyword parameters.

        INPUT:

        - ``nef_partition`` -- a `k`-part :class:`nef-partition
          <sage.geometry.lattice_polytope.NefPartition>` of `\Delta^\circ`, all
          other parameters (if given) must be lists of length `k`

        - ``monomial_points`` -- the `i`-th element of this list is either a
          list of integers or a string. A list will be interpreted as indices
          of points of `\Delta_i` which should be used for monomials of the
          `i`-th polynomial of this complete intersection. A string must be one
          of the following descriptions of points of `\Delta_i`:

          * "vertices",
          * "vertices+origin",
          * "all" (default),

          when using this description, it is also OK to pass a single string as
          ``monomial_points`` instead of repeating it `k` times.

        - ``coefficient_names`` -- the `i`-th element of this list specifies
          names for the monomial coefficients of the `i`-th polynomial, see
          :func:`~sage.schemes.toric.variety.normalize_names`
          for acceptable formats. If not given, indexed coefficient names will
          be created automatically.

        - ``coefficient_name_indices`` -- the `i`-th element of this list
          specifies indices for indexed coefficients of the `i`-th polynomial.
          If not given, the index of each coefficient will coincide with the
          index of the corresponding point of `\Delta_i`.

        - ``coefficients`` -- as an alternative to specifying coefficient
          names and/or indices, you can give the coefficients themselves as
          arbitrary expressions and/or strings. Using strings allows you to
          easily add "parameters": the base field of ``self`` will be extended
          to include all necessary names.

        OUTPUT:

        - a :class:`nef complete intersection <NefCompleteIntersection>` of
          ``self`` (with the extended base field, if necessary).

        EXAMPLES:

        We construct several complete intersections associated to the same
        nef-partition of the 3-dimensional reflexive polytope #2254::

            sage: p = ReflexivePolytope(3, 2254)
            sage: np = p.nef_partitions()[1]; np
            Nef-partition {2, 3, 4, 7, 8} ⊔ {0, 1, 5, 6}
            sage: X = FanoToricVariety(Delta_polar=p)
            sage: X.nef_complete_intersection(np)
            Closed subscheme of 3-d Fano toric variety
             covered by 10 affine patches defined by:
              a0*z1*z4^2*z5^2*z7^3 + a2*z2*z4*z5*z6*z7^2*z8^2
              + a3*z2*z3*z4*z7*z8 + a1*z0*z2,
              b3*z1*z4*z5^2*z6^2*z7^2*z8^2 + b0*z2*z5*z6^3*z7*z8^4
              + b5*z1*z3*z4*z5*z6*z7*z8 + b2*z2*z3*z6^2*z8^3
              + b1*z1*z3^2*z4 + b4*z0*z1*z5*z6

        Now we include only monomials associated to vertices of `\Delta_i`::

            sage: X.nef_complete_intersection(np, monomial_points='vertices')
            Closed subscheme of 3-d Fano toric variety
             covered by 10 affine patches defined by:
              a0*z1*z4^2*z5^2*z7^3 + a2*z2*z4*z5*z6*z7^2*z8^2
              + a3*z2*z3*z4*z7*z8 + a1*z0*z2,
              b3*z1*z4*z5^2*z6^2*z7^2*z8^2 + b0*z2*z5*z6^3*z7*z8^4
              + b2*z2*z3*z6^2*z8^3 + b1*z1*z3^2*z4 + b4*z0*z1*z5*z6

        (effectively, we set ``b5=0``). Next we provide coefficients explicitly
        instead of using default generic names::

            sage: X.nef_complete_intersection(np,
            ....:       monomial_points='vertices',
            ....:       coefficients=[("a", "a^2", "a/e", "c_i"), list(range(1,6))])
            Closed subscheme of 3-d Fano toric variety
             covered by 10 affine patches defined by:
              a*z1*z4^2*z5^2*z7^3 + a/e*z2*z4*z5*z6*z7^2*z8^2
              + (c_i)*z2*z3*z4*z7*z8 + (a^2)*z0*z2,
              4*z1*z4*z5^2*z6^2*z7^2*z8^2 + z2*z5*z6^3*z7*z8^4
              + 3*z2*z3*z6^2*z8^3 + 2*z1*z3^2*z4 + 5*z0*z1*z5*z6

        Finally, we take a look at the generic representative of these complete
        intersections in a completely resolved ambient toric variety::

            sage: X = FanoToricVariety(Delta_polar=p,
            ....:                         coordinate_points='all')
            sage: X.nef_complete_intersection(np)
            Closed subscheme of 3-d Fano toric variety
             covered by 22 affine patches defined by:
              a2*z2*z4*z5*z6*z7^2*z8^2*z9^2*z10^2*z11*z12*z13
              + a0*z1*z4^2*z5^2*z7^3*z9*z10^2*z12*z13
              + a3*z2*z3*z4*z7*z8*z9*z10*z11*z12 + a1*z0*z2,
              b0*z2*z5*z6^3*z7*z8^4*z9^3*z10^2*z11^2*z12*z13^2
              + b3*z1*z4*z5^2*z6^2*z7^2*z8^2*z9^2*z10^2*z11*z12*z13^2
              + b2*z2*z3*z6^2*z8^3*z9^2*z10*z11^2*z12*z13
              + b5*z1*z3*z4*z5*z6*z7*z8*z9*z10*z11*z12*z13
              + b1*z1*z3^2*z4*z11*z12 + b4*z0*z1*z5*z6*z13
        """
        return NefCompleteIntersection(self, nef_partition, **kwds)

    def cartesian_product(self, other,
                          coordinate_names=None, coordinate_indices=None):
        r"""
        Return the Cartesian product of ``self`` with ``other``.

        INPUT:

        - ``other`` -- a (possibly
          :class:`CPR-Fano <CPRFanoToricVariety_field>`) :class:`toric variety
          <sage.schemes.toric.variety.ToricVariety_field>`

        - ``coordinate_names`` -- names of variables for the coordinate ring,
          see :func:`normalize_names` for acceptable formats. If not given,
          indexed variable names will be created automatically.

        - ``coordinate_indices`` -- list of integers, indices for indexed
          variables. If not given, the index of each variable will coincide
          with the index of the corresponding ray of the fan.

        OUTPUT:

        - a :class:`toric variety
          <sage.schemes.toric.variety.ToricVariety_field>`, which is
          :class:`CPR-Fano <CPRFanoToricVariety_field>` if ``other`` was.

        EXAMPLES::

            sage: P1 = toric_varieties.P1()
            sage: P2 = toric_varieties.P2()
            sage: P1xP2 = P1.cartesian_product(P2); P1xP2
            3-d CPR-Fano toric variety covered by 6 affine patches
            sage: P1xP2.fan().rays()
            N+N( 1,  0,  0),        N+N(-1,  0,  0),        N+N( 0,  1,  0),
            N+N( 0,  0,  1),        N+N( 0, -1, -1)
            in 3-d lattice N+N
            sage: P1xP2.Delta_polar()
            3-d reflexive polytope in 3-d lattice N+N
        """
        if isinstance(other, FanoToricVariety_field) and not isinstance(other, CPRFanoToricVariety_field) \
            and not isinstance(self, SmoothFanoToricVariety_field):
            fan = self.fan().cartesian_product(other.fan())
            Delta_polar = LatticePolytope(fan.rays())
            return FanoToricVariety_field(Delta_polar, fan,
                                        coordinate_names, coordinate_indices,
                                        self.base_ring())
        return super().cartesian_product(other)

class AnticanonicalHypersurface(AlgebraicScheme_subscheme_toric):
    r"""
    Construct an anticanonical hypersurface of a CPR-Fano toric variety.

    INPUT:

    - ``P_Delta`` -- :class:`CPR-Fano toric variety
      <CPRFanoToricVariety_field>` associated to a reflexive polytope `\Delta`

    - see :meth:`CPRFanoToricVariety_field.anticanonical_hypersurface` for
      documentation on all other acceptable parameters

    OUTPUT:

    :class:`anticanonical hypersurface <AnticanonicalHypersurface>` of
    ``P_Delta`` (with the extended base field, if necessary).

    EXAMPLES::

        sage: P1xP1 = toric_varieties.P1xP1()
        sage: import sage.schemes.toric.fano_variety as ftv
        sage: ftv.AnticanonicalHypersurface(P1xP1)
        Closed subscheme of 2-d CPR-Fano toric variety
         covered by 4 affine patches defined by:
          a0*s^2*x^2 + a3*t^2*x^2 + a6*s*t*x*y + a1*s^2*y^2 + a2*t^2*y^2

    See :meth:`~CPRFanoToricVariety_field.anticanonical_hypersurface()` for a
    more elaborate example.
    """
    def __init__(self, P_Delta, monomial_points=None, coefficient_names=None,
                 coefficient_name_indices=None, coefficients=None):
        r"""
        See :meth:`CPRFanoToricVariety_field.anticanonical_hypersurface` for
        documentation.

        TESTS::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: import sage.schemes.toric.fano_variety as ftv
            sage: ftv.AnticanonicalHypersurface(P1xP1)
            Closed subscheme of 2-d CPR-Fano toric variety
             covered by 4 affine patches defined by:
              a0*s^2*x^2 + a3*t^2*x^2 + a6*s*t*x*y + a1*s^2*y^2 + a2*t^2*y^2

        Check that finite fields are handled correctly :issue:`14899`::

            sage: F = GF(5^2, "a")                                                      # needs sage.rings.finite_rings
            sage: X = P1xP1.change_ring(F)                                              # needs sage.rings.finite_rings
            sage: X.anticanonical_hypersurface(monomial_points='all',                   # needs sage.rings.finite_rings
            ....:                   coefficients=[1]*X.Delta().npoints())
            Closed subscheme of 2-d CPR-Fano toric variety
             covered by 4 affine patches defined by:
              s^2*x^2 + s*t*x^2 + t^2*x^2 + s^2*x*y + s*t*x*y
              + t^2*x*y + s^2*y^2 + s*t*y^2 + t^2*y^2
        """
        if not isinstance(P_Delta, CPRFanoToricVariety_field):
            raise TypeError("anticanonical hypersurfaces can only be "
                            "constructed for CPR-Fano toric varieties!"
                            "\nGot: %s" % P_Delta)
        Delta = P_Delta.Delta()
        Delta_polar = Delta.polar()
        # Monomial points normalization
        if monomial_points == "vertices":
            monomial_points = list(range(Delta.nvertices()))
        elif monomial_points == "all":
            monomial_points = list(range(Delta.npoints()))
        elif monomial_points == "vertices+origin":
            monomial_points = list(range(Delta.nvertices()))
            monomial_points.append(Delta.origin())
        elif monomial_points == "simplified" or monomial_points is None:
            monomial_points = Delta.skeleton_points(Delta.dim() - 2)
            monomial_points.append(Delta.origin())
        elif isinstance(monomial_points, str):
            raise ValueError("%s is an unsupported description of monomial "
                             "points!" % monomial_points)
        monomial_points = tuple(monomial_points)
        self._monomial_points = monomial_points
        # Make the necessary ambient space
        if coefficients is None:
            if coefficient_name_indices is None:
                coefficient_name_indices = monomial_points
            coefficient_names = normalize_names(
                                coefficient_names, len(monomial_points),
                                DEFAULT_COEFFICIENT, coefficient_name_indices)
            # We probably don't want it: the analog in else-branch is unclear.
            # self._coefficient_names = coefficient_names
            F = add_variables(P_Delta.base_ring(), coefficient_names)
            coefficients = [F(coef) for coef in coefficient_names]
        else:
            variables = set()
            nonstr = []
            regex = re.compile(r"[_A-Za-z]\w*")
            for c in coefficients:
                if isinstance(c, str):
                    variables.update(regex.findall(c))
                else:
                    nonstr.append(c)
            F = add_variables(P_Delta.base_ring(), sorted(variables))
            F = coercion_model.common_parent(F, *nonstr)
            coefficients = [F(_) for _ in coefficients]
        P_Delta = P_Delta.base_extend(F)
        if len(monomial_points) != len(coefficients):
            raise ValueError("cannot construct equation of the anticanonical"
                     " hypersurface with %d monomials and %d coefficients"
                     % (len(monomial_points), len(coefficients)))
        # Defining polynomial
        h = sum(coef * prod(P_Delta.coordinate_point_to_coordinate(n)
                            ** (Delta.point(m) * Delta_polar.point(n) + 1)
                            for n in P_Delta.coordinate_points())
            for m, coef in zip(monomial_points, coefficients))
        super().__init__(P_Delta, h)


class NefCompleteIntersection(AlgebraicScheme_subscheme_toric):
    r"""
    Construct a nef complete intersection in a Fano toric variety.

    INPUT:

    - ``P_Delta`` -- a :class:`Fano toric variety
      <FanoToricVariety_field>` associated to a reflexive polytope `\Delta`

    - see :meth:`FanoToricVariety_field.nef_complete_intersection` for
      documentation on all other acceptable parameters

    OUTPUT:

    - a :class:`nef complete intersection <NefCompleteIntersection>` of
      ``P_Delta`` (with the extended base field, if necessary).

    EXAMPLES::

        sage: o = lattice_polytope.cross_polytope(3)
        sage: np = o.nef_partitions()[0]; np
        Nef-partition {0, 1, 3} ⊔ {2, 4, 5}
        sage: X = CPRFanoToricVariety(Delta_polar=o)
        sage: X.nef_complete_intersection(np)
        Closed subscheme of 3-d Fano toric variety
         covered by 8 affine patches defined by:
          a2*z0^2*z1 + a5*z0*z1*z3 + a1*z1*z3^2 + a3*z0^2*z4 + a4*z0*z3*z4 + a0*z3^2*z4,
          b1*z1*z2^2 + b2*z2^2*z4 + b5*z1*z2*z5 + b4*z2*z4*z5 + b3*z1*z5^2 + b0*z4*z5^2

    See :meth:`FanoToricVariety_field.nef_complete_intersection` for a
    more elaborate example.
    """
    def __init__(self, P_Delta, nef_partition,
                 monomial_points='all', coefficient_names=None,
                 coefficient_name_indices=None, coefficients=None):
        r"""
        See :meth:`FanoToricVariety_field.nef_complete_intersection` for
        documentation.

        TESTS::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: np = o.nef_partitions()[0]
            sage: np
            Nef-partition {0, 1, 3} ⊔ {2, 4, 5}
            sage: X = FanoToricVariety(Delta_polar=o)
            sage: from sage.schemes.toric.fano_variety import *
            sage: NefCompleteIntersection(X, np)
            Closed subscheme of 3-d Fano toric variety
             covered by 8 affine patches defined by:
              a2*z0^2*z1 + a5*z0*z1*z3 + a1*z1*z3^2
              + a3*z0^2*z4 + a4*z0*z3*z4 + a0*z3^2*z4,
              b1*z1*z2^2 + b2*z2^2*z4 + b5*z1*z2*z5
              + b4*z2*z4*z5 + b3*z1*z5^2 + b0*z4*z5^2
        """
        if not isinstance(P_Delta, FanoToricVariety_field):
            raise TypeError("nef complete intersections can only be "
                            "constructed for CPR-Fano toric varieties!"
                            "\nGot: %s" % P_Delta)
        if nef_partition.Delta() is not P_Delta.Delta():
            raise ValueError("polytopes 'Delta' of the nef-partition and the "
                             "Fano toric variety must be the same!")
        self._nef_partition = nef_partition
        k = nef_partition.nparts()
        # Pre-normalize all parameters
        if isinstance(monomial_points, str):
            monomial_points = [monomial_points] * k
        if coefficient_names is None:
            coefficient_names = [None] * k
        if coefficient_name_indices is None:
            coefficient_name_indices = [None] * k
        if coefficients is None:
            coefficients = [None] * k

        polynomials = []
        Delta_polar = P_Delta.Delta_polar()
        for i in range(k):
            Delta_i = nef_partition.Delta(i)
            # Monomial points normalization
            if monomial_points[i] == "vertices":
                monomial_points[i] = list(range(Delta_i.nvertices()))
            elif monomial_points[i] == "all":
                monomial_points[i] = list(range(Delta_i.npoints()))
            elif monomial_points[i] == "vertices+origin":
                monomial_points[i] = list(range(Delta_i.nvertices()))
                if (Delta_i.origin() is not None
                    and Delta_i.origin() >= Delta_i.nvertices()):
                    monomial_points[i].append(Delta_i.origin())
            elif isinstance(monomial_points[i], str):
                raise ValueError("'%s' is an unsupported description of "
                                 "monomial points!" % monomial_points[i])
            monomial_points[i] = tuple(monomial_points[i])
            # Extend the base ring of the ambient space if necessary
            if coefficients[i] is None:
                if coefficient_name_indices[i] is None:
                    coefficient_name_indices[i] = monomial_points[i]
                coefficient_names[i] = normalize_names(
                        coefficient_names[i], len(monomial_points[i]),
                        DEFAULT_COEFFICIENTS[i], coefficient_name_indices[i])
                F = add_variables(P_Delta.base_ring(), coefficient_names[i])
                coefficients[i] = [F(coef) for coef in coefficient_names[i]]
            else:
                variables = set()
                nonstr = []
                regex = re.compile(r"[_A-Za-z]\w*")
                for c in coefficients[i]:
                    if isinstance(c, str):
                        variables.update(regex.findall(c))
                    else:
                        nonstr.append(c)
                F = add_variables(P_Delta.base_ring(), sorted(variables))
                F = coercion_model.common_parent(F, *nonstr)
                coefficients[i] = [F(_) for _ in coefficients[i]]
            P_Delta = P_Delta.base_extend(F)
            if len(monomial_points[i]) != len(coefficients[i]):
                raise ValueError("cannot construct equation %d of the complete"
                         " intersection with %d monomials and %d coefficients"
                         % (i, len(monomial_points[i]), len(coefficients[i])))
            # Defining polynomial
            h = sum(coef * prod(P_Delta.coordinate_point_to_coordinate(n)
                                ** (Delta_i.point(m) * Delta_polar.point(n)
                                    + (nef_partition.part_of_point(n) == i))
                                for n in P_Delta.coordinate_points())
                for m, coef in zip(monomial_points[i], coefficients[i]))
            polynomials.append(h)
        self._monomial_points = tuple(monomial_points)
        super().__init__(P_Delta, polynomials)

    def cohomology_class(self):
        r"""
        Return the class of ``self`` in the ambient space cohomology ring.

        OUTPUT: a :class:`cohomology class <sage.schemes.generic.toric_variety.CohomologyClass>`

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: np = o.nef_partitions()[0]; np
            Nef-partition {0, 1, 3} ⊔ {2, 4, 5}
            sage: X = FanoToricVariety(Delta_polar=o)
            sage: CI = X.nef_complete_intersection(np); CI
            Closed subscheme of 3-d Fano toric variety
             covered by 8 affine patches defined by:
              a2*z0^2*z1 + a5*z0*z1*z3 + a1*z1*z3^2 + a3*z0^2*z4 + a4*z0*z3*z4 + a0*z3^2*z4,
              b1*z1*z2^2 + b2*z2^2*z4 + b5*z1*z2*z5 + b4*z2*z4*z5 + b3*z1*z5^2 + b0*z4*z5^2
            sage: CI.cohomology_class()                                                 # needs sage.libs.singular
            [2*z3*z4 + 4*z3*z5 + 2*z4*z5]
        """
        X = self.ambient_space()
        H = X.cohomology_ring()
        return prod(sum(H.gen(X._point_to_ray[point])
                    for point in part if point in X._coordinate_points)
               for part in self.nef_partition().parts(all_points=True))

    def nef_partition(self):
        r"""
        Return the nef-partition associated to ``self``.

        OUTPUT: a :class:`nef-partition <sage.geometry.lattice_polytope.NefPartition>`

        EXAMPLES::

            sage: o = lattice_polytope.cross_polytope(3)
            sage: np = o.nef_partitions()[0]; np
            Nef-partition {0, 1, 3} ⊔ {2, 4, 5}
            sage: X = FanoToricVariety(Delta_polar=o)
            sage: CI = X.nef_complete_intersection(np); CI
            Closed subscheme of 3-d Fano toric variety
             covered by 8 affine patches defined by:
              a2*z0^2*z1 + a5*z0*z1*z3 + a1*z1*z3^2 + a3*z0^2*z4 + a4*z0*z3*z4 + a0*z3^2*z4,
              b1*z1*z2^2 + b2*z2^2*z4 + b5*z1*z2*z5 + b4*z2*z4*z5 + b3*z1*z5^2 + b0*z4*z5^2
            sage: CI.nef_partition()
            Nef-partition {0, 1, 3} ⊔ {2, 4, 5}
            sage: CI.nef_partition() is np
            True
        """
        return self._nef_partition


def add_variables(field, variables):
    r"""
    Extend ``field`` to include all ``variables``.

    INPUT:

    - ``field`` -- a field

    - ``variables`` -- list of strings

    OUTPUT:

    - a fraction field extending the original ``field``, which has all
      ``variables`` among its generators.

    EXAMPLES:

    We start with the rational field and slowly add more variables::

        sage: from sage.schemes.toric.fano_variety import *
        sage: F = add_variables(QQ, []); F      # No extension
        Rational Field
        sage: F = add_variables(QQ, ["a"]); F
        Fraction Field of Univariate Polynomial Ring in a over Rational Field
        sage: F = add_variables(F, ["a"]); F
        Fraction Field of Univariate Polynomial Ring in a over Rational Field
        sage: F = add_variables(F, ["b", "c"]); F
        Fraction Field of Multivariate Polynomial Ring in a, b, c over Rational Field
        sage: F = add_variables(F, ["c", "d", "b", "c", "d"]); F
        Fraction Field of Multivariate Polynomial Ring in a, b, c, d over Rational Field
    """
    if not variables:
        return field
    if isinstance(field, FractionField_generic):
        # Q(a) ---> Q(a, b) rather than Q(a)(b)
        R = field.ring()
        if isinstance(R, (PolynomialRing_generic, MPolynomialRing_base)):
            new_variables = list(R.variable_names())
            for v in variables:
                if v not in new_variables:
                    new_variables.append(v)
            if len(new_variables) > R.ngens():
                return PolynomialRing(R.base_ring(),
                                      new_variables).fraction_field()
            else:
                return field
    # "Intelligent extension" didn't work, use the "usual one."
    new_variables = []
    for v in variables:
        if v not in new_variables:
            new_variables.append(v)
    return PolynomialRing(field, new_variables).fraction_field()