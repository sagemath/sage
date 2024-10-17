r"""
Continuous Maps Between Topological Manifolds

:class:`ContinuousMap` implements continuous maps from a topological
manifold `M` to some topological manifold `N` over the same topological
field `K` as `M`.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015): initial version
- Travis Scrimshaw (2016): review tweaks

REFERENCES:

- Chap. 1 of [KN1963]_
- [Lee2011]_
"""

# ****************************************************************************
#       Copyright (C) 2015-2019 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#       Copyright (C) 2016 Travis Scrimshaw <tscrimsh@umn.edu>
#       Copyright (C) 2018 Florentin Jaffredo <florentin.jaffredo@polytechnique.edu>
#       Copyright (C) 2020 Michael Jung <micjung@uni-potsdam.de>
#       Copyright (C) 2021 Matthias Koeppe <mkoeppe@math.ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.homset import Hom
from sage.categories.morphism import Morphism

class ContinuousMap(Morphism):
    r"""
    Continuous map between two topological manifolds.

    This class implements continuous maps of the type

    .. MATH::

        \Phi: M \longrightarrow N,

    where `M` and `N` are topological manifolds over the same
    topological field `K`.

    Continuous maps are the morphisms of the category of topological
    manifolds. The set of all continuous maps from `M` to `N` is
    therefore the homset between `M` and `N`, which is denoted
    by `\mathrm{Hom}(M,N)`.

    The class :class:`ContinuousMap` is a Sage *element* class,
    whose *parent* class is
    :class:`~sage.manifolds.manifold_homset.TopologicalManifoldHomset`.

    INPUT:

    - ``parent`` -- homset `\mathrm{Hom}(M,N)` to which the continuous
      map belongs
    - ``coord_functions`` -- dictionary of the coordinate expressions
      (as lists or tuples of the coordinates of the image expressed in
      terms of the coordinates of the considered point) with the pairs
      of charts ``(chart1, chart2)`` as keys (``chart1`` being a chart
      on `M` and ``chart2`` a chart on `N`)
    - ``name`` -- (default: ``None``) name given to ``self``
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
      continuous map; if ``None``, the LaTeX symbol is set to
      ``name``
    - ``is_isomorphism`` -- boolean (default: ``False``); determines whether the
      constructed object is a isomorphism (i.e. a homeomorphism). If set to
      ``True``, then the manifolds `M` and `N` must have the same dimension.
    - ``is_identity`` -- boolean (default: ``False``); determines whether the
      constructed object is the identity map. If set to ``True``,
      then `N` must be `M` and the entry ``coord_functions`` is not used.

    .. NOTE::

        If the information passed by means of the argument
        ``coord_functions`` is not sufficient to fully specify the
        continuous map, further coordinate expressions, in other charts,
        can be subsequently added by means of the method :meth:`add_expr`.

    EXAMPLES:

    The standard embedding of the sphere `S^2` into `\RR^3`::

        sage: M = Manifold(2, 'S^2', structure='topological') # the 2-dimensional sphere S^2
        sage: U = M.open_subset('U') # complement of the North pole
        sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
        sage: V = M.open_subset('V') # complement of the South pole
        sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
        sage: M.declare_union(U,V)   # S^2 is the union of U and V
        sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
        ....:                                intersection_name='W',
        ....:                                restrictions1=x^2+y^2!=0,
        ....:                                restrictions2=u^2+v^2!=0)
        sage: uv_to_xy = xy_to_uv.inverse()
        sage: N = Manifold(3, 'R^3', latex_name=r'\RR^3', structure='topological')  # R^3
        sage: c_cart.<X,Y,Z> = N.chart()  # Cartesian coordinates on R^3
        sage: Phi = M.continuous_map(N,
        ....:   {(c_xy, c_cart): [2*x/(1+x^2+y^2), 2*y/(1+x^2+y^2), (x^2+y^2-1)/(1+x^2+y^2)],
        ....:    (c_uv, c_cart): [2*u/(1+u^2+v^2), 2*v/(1+u^2+v^2), (1-u^2-v^2)/(1+u^2+v^2)]},
        ....:   name='Phi', latex_name=r'\Phi')
        sage: Phi
        Continuous map Phi from the 2-dimensional topological manifold S^2
         to the 3-dimensional topological manifold R^3
        sage: Phi.parent()
        Set of Morphisms from 2-dimensional topological manifold S^2
         to 3-dimensional topological manifold R^3
         in Category of manifolds over Real Field with 53 bits of precision
        sage: Phi.parent() is Hom(M, N)
        True
        sage: type(Phi)
        <class 'sage.manifolds.manifold_homset.TopologicalManifoldHomset_with_category.element_class'>
        sage: Phi.display()
        Phi: S^2 → R^3
        on U: (x, y) ↦ (X, Y, Z) = (2*x/(x^2 + y^2 + 1), 2*y/(x^2 + y^2 + 1), (x^2 + y^2 - 1)/(x^2 + y^2 + 1))
        on V: (u, v) ↦ (X, Y, Z) = (2*u/(u^2 + v^2 + 1), 2*v/(u^2 + v^2 + 1), -(u^2 + v^2 - 1)/(u^2 + v^2 + 1))

    It is possible to create the map using
    :meth:`~sage.manifolds.manifold.TopologicalManifold.continuous_map`
    with only in a single pair of charts. The argument ``coord_functions``
    is then a mere list of coordinate expressions (and not a dictionary)
    and the arguments ``chart1`` and ``chart2`` have to be provided if
    the charts differ from the default ones on the domain and/or codomain::

        sage: Phi1 = M.continuous_map(N, [2*x/(1+x^2+y^2), 2*y/(1+x^2+y^2), (x^2+y^2-1)/(1+x^2+y^2)],
        ....:                         chart1=c_xy, chart2=c_cart,
        ....:                         name='Phi', latex_name=r'\Phi')

    Since ``c_xy`` and ``c_cart`` are the default charts on respectively
    ``M`` and ``N``, they can be omitted, so that the above declaration
    is equivalent to::

        sage: Phi1 = M.continuous_map(N, [2*x/(1+x^2+y^2), 2*y/(1+x^2+y^2), (x^2+y^2-1)/(1+x^2+y^2)],
        ....:                         name='Phi', latex_name=r'\Phi')

    With such a declaration, the continuous map ``Phi1`` is only partially
    defined on the manifold `S^2` as it is known in only one chart::

        sage: Phi1.display()
        Phi: S^2 → R^3
        on U: (x, y) ↦ (X, Y, Z) = (2*x/(x^2 + y^2 + 1), 2*y/(x^2 + y^2 + 1), (x^2 + y^2 - 1)/(x^2 + y^2 + 1))

    The definition can be completed by using :meth:`add_expr`::

        sage: Phi1.add_expr(c_uv, c_cart, [2*u/(1+u^2+v^2), 2*v/(1+u^2+v^2), (1-u^2-v^2)/(1+u^2+v^2)])
        sage: Phi1.display()
        Phi: S^2 → R^3
        on U: (x, y) ↦ (X, Y, Z) = (2*x/(x^2 + y^2 + 1), 2*y/(x^2 + y^2 + 1), (x^2 + y^2 - 1)/(x^2 + y^2 + 1))
        on V: (u, v) ↦ (X, Y, Z) = (2*u/(u^2 + v^2 + 1), 2*v/(u^2 + v^2 + 1), -(u^2 + v^2 - 1)/(u^2 + v^2 + 1))

    At this stage, ``Phi1`` and ``Phi`` are fully equivalent::

        sage: Phi1 == Phi
        True

    The map acts on points::

        sage: np = M.point((0,0), chart=c_uv)  # the North pole
        sage: Phi(np)
        Point on the 3-dimensional topological manifold R^3
        sage: Phi(np).coord()  # Cartesian coordinates
        (0, 0, 1)
        sage: sp = M.point((0,0), chart=c_xy)  # the South pole
        sage: Phi(sp).coord()  # Cartesian coordinates
        (0, 0, -1)

    The test suite is passed::

        sage: TestSuite(Phi).run()
        sage: TestSuite(Phi1).run()

    Continuous maps can be composed by means of the operator ``*``.
    Let us introduce the map `\RR^3 \to \RR^2` corresponding to
    the projection from the point `(X, Y, Z) = (0, 0, 1)` onto the
    equatorial plane `Z = 0`::

        sage: P = Manifold(2, 'R^2', latex_name=r'\RR^2', structure='topological') # R^2 (equatorial plane)
        sage: cP.<xP, yP> = P.chart()
        sage: Psi = N.continuous_map(P, (X/(1-Z), Y/(1-Z)), name='Psi',
        ....:                      latex_name=r'\Psi')
        sage: Psi
        Continuous map Psi from the 3-dimensional topological manifold R^3
         to the 2-dimensional topological manifold R^2
        sage: Psi.display()
        Psi: R^3 → R^2
           (X, Y, Z) ↦ (xP, yP) = (-X/(Z - 1), -Y/(Z - 1))

    Then we compose ``Psi`` with ``Phi``, thereby getting a map
    `S^2 \to \RR^2`::

        sage: ster = Psi * Phi ; ster
        Continuous map from the 2-dimensional topological manifold S^2
         to the 2-dimensional topological manifold R^2

    Let us test on the South pole (``sp``) that ``ster`` is indeed the
    composite of ``Psi`` and ``Phi``::

        sage: ster(sp) == Psi(Phi(sp))
        True

    Actually ``ster`` is the stereographic projection from the North pole,
    as its coordinate expression reveals::

        sage: ster.display()
        S^2 → R^2
        on U: (x, y) ↦ (xP, yP) = (x, y)
        on V: (u, v) ↦ (xP, yP) = (u/(u^2 + v^2), v/(u^2 + v^2))

    If the codomain of a continuous map is 1-dimensional, the map can
    be defined by a single symbolic expression for each pair of charts
    and not by a list/tuple with a single element::

        sage: N = Manifold(1, 'N', structure='topological')
        sage: c_N = N.chart('X')
        sage: Phi = M.continuous_map(N, {(c_xy, c_N): x^2+y^2,
        ....:                            (c_uv, c_N): 1/(u^2+v^2)})

        sage: Psi = M.continuous_map(N, {(c_xy, c_N): [x^2+y^2],
        ....:                            (c_uv, c_N): [1/(u^2+v^2)]})
        sage: Phi == Psi
        True

    Next we construct an example of continuous map `\RR \to \RR^2`::

        sage: R = Manifold(1, 'R', structure='topological')  # field R
        sage: T.<t> = R.chart()  # canonical chart on R
        sage: R2 = Manifold(2, 'R^2', structure='topological')  # R^2
        sage: c_xy.<x,y> = R2.chart() # Cartesian coordinates on R^2
        sage: Phi = R.continuous_map(R2, [cos(t), sin(t)], name='Phi'); Phi
        Continuous map Phi from the 1-dimensional topological manifold R
         to the 2-dimensional topological manifold R^2
        sage: Phi.parent()
        Set of Morphisms from 1-dimensional topological manifold R
         to 2-dimensional topological manifold R^2
         in Category of manifolds over Real Field with 53 bits of precision
        sage: Phi.parent() is Hom(R, R2)
        True
        sage: Phi.display()
        Phi: R → R^2
           t ↦ (x, y) = (cos(t), sin(t))

    An example of homeomorphism between the unit open disk and the
    Euclidean plane `\RR^2`::

        sage: D = R2.open_subset('D', coord_def={c_xy: x^2+y^2<1}) # the open unit disk
        sage: Phi = D.homeomorphism(R2, [x/sqrt(1-x^2-y^2), y/sqrt(1-x^2-y^2)],
        ....:                       name='Phi', latex_name=r'\Phi')
        sage: Phi
        Homeomorphism Phi from the Open subset D of the 2-dimensional
         topological manifold R^2 to the 2-dimensional topological manifold R^2
        sage: Phi.parent()
        Set of Morphisms from Open subset D of the 2-dimensional topological
         manifold R^2 to 2-dimensional topological manifold R^2 in Category of
         manifolds over Real Field with 53 bits of precision
        sage: Phi.parent() is Hom(D, R2)
        True
        sage: Phi.display()
        Phi: D → R^2
           (x, y) ↦ (x, y) = (x/sqrt(-x^2 - y^2 + 1), y/sqrt(-x^2 - y^2 + 1))

    The image of a point::

        sage: p = D.point((1/2,0))
        sage: q = Phi(p) ; q
        Point on the 2-dimensional topological manifold R^2
        sage: q.coord()
        (1/3*sqrt(3), 0)

    The inverse homeomorphism is computed by :meth:`inverse`::

        sage: Phi.inverse()
        Homeomorphism Phi^(-1) from the 2-dimensional topological manifold R^2
         to the Open subset D of the 2-dimensional topological manifold R^2
        sage: Phi.inverse().display()
        Phi^(-1): R^2 → D
           (x, y) ↦ (x, y) = (x/sqrt(x^2 + y^2 + 1), y/sqrt(x^2 + y^2 + 1))

    Equivalently, one may use the notations ``^(-1)`` or ``~`` to
    get the inverse::

        sage: Phi^(-1) is Phi.inverse()
        True
        sage: ~Phi is Phi.inverse()
        True

    Check that ``~Phi`` is indeed the inverse of ``Phi``::

        sage: (~Phi)(q) == p
        True
        sage: Phi * ~Phi == R2.identity_map()
        True
        sage: ~Phi * Phi == D.identity_map()
        True

    The coordinate expression of the inverse homeomorphism::

        sage: (~Phi).display()
        Phi^(-1): R^2 → D
           (x, y) ↦ (x, y) = (x/sqrt(x^2 + y^2 + 1), y/sqrt(x^2 + y^2 + 1))

    A special case of homeomorphism: the identity map of the open unit disk::

        sage: id = D.identity_map() ; id
        Identity map Id_D of the Open subset D of the 2-dimensional topological
         manifold R^2
        sage: latex(id)
        \mathrm{Id}_{D}
        sage: id.parent()
        Set of Morphisms from Open subset D of the 2-dimensional topological
         manifold R^2 to Open subset D of the 2-dimensional topological
         manifold R^2 in Join of Category of subobjects of sets and Category of
         manifolds over Real Field with 53 bits of precision
        sage: id.parent() is Hom(D, D)
        True
        sage: id is Hom(D,D).one()  # the identity element of the monoid Hom(D,D)
        True

    The identity map acting on a point::

        sage: id(p)
        Point on the 2-dimensional topological manifold R^2
        sage: id(p) == p
        True
        sage: id(p) is p
        True

    The coordinate expression of the identity map::

        sage: id.display()
        Id_D: D → D
           (x, y) ↦ (x, y)

    The identity map is its own inverse::

        sage: id^(-1) is id
        True
        sage: ~id is id
        True
    """
    def __init__(self, parent, coord_functions=None, name=None, latex_name=None,
                 is_isomorphism=False, is_identity=False):
        r"""
        Initialize ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: N = Manifold(3, 'N', structure='topological')
            sage: Y.<u,v,w> = N.chart()
            sage: f = Hom(M,N)({(X,Y): (x+y, x*y, x-y)}, name='f') ; f
            Continuous map f from the 2-dimensional topological manifold M
             to the 3-dimensional topological manifold N
            sage: f.display()
            f: M → N
               (x, y) ↦ (u, v, w) = (x + y, x*y, x - y)
            sage: TestSuite(f).run()

        The identity map::

            sage: f = Hom(M,M)({}, is_identity=True) ; f
            Identity map Id_M of the 2-dimensional topological manifold M
            sage: f.display()
            Id_M: M → M
               (x, y) ↦ (x, y)
            sage: TestSuite(f).run()
        """
        Morphism.__init__(self, parent)
        domain = parent.domain()
        codomain = parent.codomain()
        self._domain = domain
        self._codomain = codomain
        self._coord_expression = {}  # dict. of coordinate expressions of the
                                     # map:
                                     # - key: pair of charts
                                     # - value: instance of MultiCoordFunction
        self._is_isomorphism = False  # default value; may be redefined below
        self._is_identity = False  # default value; may be redefined below
        if is_identity:
            # Construction of the identity map
            self._is_identity = True
            self._is_isomorphism = True
            if domain != codomain:
                raise ValueError("the domain and codomain must coincide"
                                 " for the identity map")
            if name is None:
                name = 'Id_' + domain._name
            if latex_name is None:
                latex_name = r'\mathrm{Id}_{' + domain._latex_name + r'}'
            self._name = name
            self._latex_name = latex_name
            for chart in domain.atlas():
                coord_funct = chart[:]
                self._coord_expression[(chart, chart)] = chart.multifunction(
                                                                  *coord_funct)
        else:
            # Construction of a generic continuous map
            if is_isomorphism:
                self._is_isomorphism = True
                if domain.dim() != codomain.dim():
                    raise ValueError("for an isomorphism, the source"
                                     " manifold and target manifold must"
                                     " have the same dimension")
            if coord_functions is not None:
                n2 = self._codomain.dim()
                for chart_pair, expression in coord_functions.items():
                    if chart_pair[0] not in self._domain.atlas():
                        raise ValueError("{} is not a chart ".format(
                                                              chart_pair[0]) +
                                         "defined on the {}".format(self._domain))
                    if chart_pair[1] not in self._codomain.atlas():
                        raise ValueError("{} is not a chart ".format(
                                                              chart_pair[1]) +
                                         "defined on the {}".format(self._codomain))
                    if n2 == 1:
                        # a single expression entry is allowed
                        if not isinstance(expression, (tuple, list)):
                            expression = (expression,)
                    if len(expression) != n2:
                        raise ValueError("{} coordinate ".format(n2) +
                                         "functions must be provided")
                    self._coord_expression[chart_pair] = \
                                       chart_pair[0].multifunction(*expression)
            self._name = name
            if latex_name is None:
                self._latex_name = self._name
            else:
                self._latex_name = latex_name
        self._init_derived()  # initialization of derived quantities

    #
    # SageObject methods
    #

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: N = Manifold(2, 'N', structure='topological')
            sage: Y.<u,v> = N.chart()
            sage: f = Hom(M,N)({(X,Y): (x+y,x*y)})
            sage: f
            Continuous map from the 2-dimensional topological manifold M
             to the 2-dimensional topological manifold N
            sage: f = Hom(M,N)({(X,Y): (x+y,x*y)}, name='f')
            sage: f
            Continuous map f from the 2-dimensional topological manifold M
             to the 2-dimensional topological manifold N
            sage: f = Hom(M,N)({(X,Y): (x+y,x-y)}, name='f', is_isomorphism=True)
            sage: f
            Homeomorphism f from the 2-dimensional topological manifold M
             to the 2-dimensional topological manifold N
            sage: f = Hom(M,M)({(X,X): (x+y,x-y)}, name='f', is_isomorphism=True)
            sage: f
            Homeomorphism f of the 2-dimensional topological manifold M
            sage: f = Hom(M,M)({}, name='f', is_identity=True)
            sage: f
            Identity map f of the 2-dimensional topological manifold M
        """
        if self._is_identity:
            return "Identity map {} of the {}".format(self._name, self._domain)
        if self._is_isomorphism:
            description = "Homeomorphism"
        else:
            description = "Continuous map"
        if self._name is not None:
            description += " " + self._name
        if self._domain == self._codomain:
            if self._is_isomorphism:
                description += " of the {}".format(self._domain)
            else:
                description += " from the {} to itself".format(self._domain)
        else:
            description += " from the {} to the {}".format(self._domain,
                                                           self._codomain)
        return description

    def _latex_(self):
        r"""
        LaTeX representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = Hom(M,M)({(X,X): (x+y,x*y)}, name='f')
            sage: latex(f)
            f
            sage: f = Hom(M,M)({(X,X): (x+y,x*y)}, name='f', latex_name=r'\Phi')
            sage: latex(f)
            \Phi
        """
        if self._latex_name is None:
            return r'\text{' + str(self) + r'}'
        else:
            return self._latex_name

    #
    # Hash and equality
    #

    def __hash__(self):
        """
        Hash function.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: N = Manifold(2, 'N', structure='topological')
            sage: Y.<u,v> = N.chart()
            sage: f = M.continuous_map(N, {(X,Y): (x+y,x*y)})
            sage: hash(f) == f.__hash__()
            True

        Let us check that ``f`` can be used as a dictionary key::

            sage: {f: 1}[f]
            1
        """
        return hash((self._domain, self._codomain))

    def __eq__(self, other):
        r"""
        Comparison (equality) operator.

        INPUT:

        - ``other`` -- a :class:`ContinuousMap`

        OUTPUT: ``True`` if ``self`` is equal to ``other`` and ``False`` otherwise

        TESTS::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: X.<x,y,z> = M.chart()
            sage: N = Manifold(2, 'N', structure='topological')
            sage: Y.<u,v> = N.chart()
            sage: f = M.continuous_map(N, {(X,Y): [x+y+z, 2*x*y*z]}, name='f')
            sage: g = M.continuous_map(N, {(X,Y): [x+y+z, 2*x*y*z]}, name='g')
            sage: f == g
            True
            sage: g = M.continuous_map(N, {(X,Y): [x+y+z, 1]}, name='g')
            sage: f == g
            False
        """
        if other is self:
            return True
        if not isinstance(other, type(self)):
            return False
        if self.parent() != other.parent():
            return False
        if self._is_identity:
            return other.is_identity()
        if other._is_identity:
            return self.is_identity()
        for charts, coord_functions in self._coord_expression.items():
            try:
                if coord_functions.expr() != other.expr(*charts):
                    return False
            except ValueError:
                return False
        return True

    def __ne__(self, other):
        r"""
        Inequality operator.

        INPUT:

        - ``other`` -- a :class:`ContinuousMap`

        OUTPUT:

        - ``True`` if ``self`` is different from ``other`` and
          ``False`` otherwise

        TESTS::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: X.<x,y,z> = M.chart()
            sage: N = Manifold(2, 'N', structure='topological')
            sage: Y.<u,v> = N.chart()
            sage: f = M.continuous_map(N, {(X,Y): [x+y+z, 2*x*y*z]}, name='f')
            sage: g = M.continuous_map(N, {(X,Y): [x+y+z, 2*x*y*z]}, name='g')
            sage: f != g
            False
            sage: g = M.continuous_map(N, {(X,Y): [x+y+z, 1]}, name='g')
            sage: f != g
            True
        """
        return not (self == other)

    #
    # Map methods
    #

    def _call_(self, point):
        r"""
        Compute the image of a point by ``self``.

        INPUT:

        - ``point`` -- :class:`~sage.manifolds.point.TopologicalManifoldPoint`;
          point in the domain of ``self``

        OUTPUT: image of the point by ``self``

        EXAMPLES:

        Planar rotation acting on a point::

            sage: M = Manifold(2, 'R^2', latex_name=r'\RR^2', structure='topological') # Euclidean plane
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates

        A pi/3 rotation around the origin defined in Cartesian coordinates::

            sage: rot = M.continuous_map(M, ((x - sqrt(3)*y)/2, (sqrt(3)*x + y)/2),
            ....:                        name='R')
            sage: p = M.point((1,2), name='p')
            sage: q = rot(p) ; q
            Point R(p) on the 2-dimensional topological manifold R^2
            sage: q.coord()
            (-sqrt(3) + 1/2, 1/2*sqrt(3) + 1)

        Image computed after some change of coordinates::

            sage: c_spher.<r,ph> = M.chart(r'r:(0,+oo) ph:(0,2*pi):\phi') # spherical coord. on the plane
            sage: ch = c_spher.transition_map(c_cart, [r*cos(ph), r*sin(ph)])
            sage: p1 = M.point((sqrt(5), arctan(2)), chart=c_spher) # p1 is defined only in terms of c_spher
            sage: q1 = rot(p1) # but the computation of the action of rot is still possible
            sage: q1 == q
            True

        Image computed by means of spherical coordinates::

            sage: rot.add_expr(c_spher, c_spher, (r, ph+pi/3)) # now rot is known in terms of c_spher
            sage: p2 = M.point((sqrt(5), arctan(2)), chart=c_spher)
            sage: q2 = rot(p2) # computation on c_spher
            sage: q2 == q
            True
        """
        # NB: checking that ``point`` belongs to the map's domain has been
        # already performed by Map.__call__(); this check is therefore not
        # repeated here.
        if self._is_identity:
            return point
        chart1, chart2 = None, None
        for chart in point._coordinates:
            for chart_pair in self._coord_expression:
                if chart_pair[0] is chart:
                    chart1 = chart
                    chart2 = chart_pair[1]
                    break
            if chart1 is not None:
                break
        else:
            # attempt to perform a change of coordinate on the point
            for chart_pair in self._coord_expression:
                try:
                    point.coord(chart_pair[0])
                    chart1, chart2 = chart_pair
                except ValueError:
                    pass
                if chart1 is not None:
                    break
            else:
                raise ValueError("no pair of charts has been found to " +
                                 "compute the action of the {} on the {}".format(self, point))
        coord_map = self._coord_expression[(chart1, chart2)]
        y = coord_map(*(point._coordinates[chart1]))
        if point._name is None or self._name is None:
            res_name = None
        else:
            res_name = self._name + '(' + point._name + ')'
        if point._latex_name is None or self._latex_name is None:
            res_latex_name = None
        else:
            res_latex_name = (self._latex_name + r'\left(' +
                              point._latex_name + r'\right)')
        # The image point is created as an element of the domain of chart2:
        dom2 = chart2.domain()
        return dom2.element_class(dom2, coords=y, chart=chart2,
                                  name=res_name, latex_name=res_latex_name,
                                  check_coords=False)
    #
    # Morphism methods
    #

    def is_identity(self):
        r"""
        Check whether ``self`` is an identity map.

        EXAMPLES:

        Tests on continuous maps of a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: M.identity_map().is_identity()  # obviously...
            True
            sage: Hom(M, M).one().is_identity()  # a variant of the obvious
            True
            sage: a = M.continuous_map(M, coord_functions={(X,X): (x, y)})
            sage: a.is_identity()
            True
            sage: a = M.continuous_map(M, coord_functions={(X,X): (x, y+1)})
            sage: a.is_identity()
            False

        Of course, if the codomain of the map does not coincide with its
        domain, the outcome is ``False``::

            sage: N = Manifold(2, 'N', structure='topological')
            sage: Y.<u,v> = N.chart()
            sage: a = M.continuous_map(N, {(X,Y): (x, y)})
            sage: a.display()
            M → N
               (x, y) ↦ (u, v) = (x, y)
            sage: a.is_identity()
            False
        """
        if self._is_identity:
            return True
        if self._codomain != self._domain:
            return False
        for chart in self._domain._top_charts:
            try:
                if chart[:] != self.expr(chart, chart):
                    return False
            except ValueError:
                return False
        # If this point is reached, ``self`` must be the identity:
        self._is_identity = True
        return True

    def _composition_(self, other, homset):
        r"""
        Composition of ``self`` with another morphism.

        The composition is performed on the right, i.e. the returned
        morphism is ``self * other``.

        INPUT:

        - ``other`` -- a continuous map whose codomain is the domain
          of ``self``
        - ``homset`` -- the homset of the continuous map ``self*other``;
          this argument is required to follow the prototype of
          :meth:`~sage.categories.map.Map._composition_` and is determined by
          :meth:`~sage.categories.map.Map._composition` (single underscore),
          that is supposed to call the current method

        OUTPUT:

        - :class:`~sage.manifolds.continuous_map.ContinuousMap` that is
          the composite map ``self * other``

        TESTS::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: X.<x,y,z> = M.chart()
            sage: N = Manifold(2, 'N', structure='topological')
            sage: Y.<u,v> = N.chart()
            sage: Q = Manifold(4, 'Q', structure='topological')
            sage: Z.<a,b,c,d> = Q.chart()
            sage: f = N.continuous_map(Q, [u+v, u*v, 1+u, 2-v])
            sage: g = M.continuous_map(N, [x+y+z, x*y*z])
            sage: s = f._composition_(g, Hom(M,Q)); s
            Continuous map from the 3-dimensional topological manifold M
             to the 4-dimensional topological manifold Q
            sage: s.display()
            M → Q
               (x, y, z) ↦ (a, b, c, d) = ((x*y + 1)*z + x + y, x*y*z^2 + (x^2*y + x*y^2)*z, x + y + z + 1, -x*y*z + 2)
            sage: s == f*g
            True
        """
        # This method is invoked by Map._composition (single underscore),
        # which is itself invoked by Map.__mul__ . The latter performs the
        # check other._codomain == self._domain. There is therefore no need
        # to perform it here.
        if self._is_identity:
            return other
        if other._is_identity:
            return self
        resu_funct = {}
        for chart1 in other._domain._top_charts:
            for chart2 in self._domain._top_charts:
                for chart3 in self._codomain._top_charts:
                    try:
                        self23 = self.coord_functions(chart2, chart3)
                        resu_funct[(chart1, chart3)] = self23(*other.expr(chart1, chart2),
                                                              simplify=True)
                    except ValueError:
                        pass
        return homset(resu_funct)

    def image(self, subset=None, inverse=None):
        r"""
        Return the image of ``self`` or the image of ``subset`` under ``self``.

        INPUT:

        - ``inverse`` -- (default: ``None``) continuous map from
          ``map.codomain()`` to ``map.domain()``, which once restricted to the image
          of `\Phi` is the inverse of `\Phi` onto its image if the latter
          exists (NB: no check of this is performed)
        - ``subset`` -- (default: the domain of ``map``) a subset of the domain of
          ``self``

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: N = Manifold(1, 'N', ambient=M, structure='topological')
            sage: CM.<x,y> = M.chart()
            sage: CN.<u> = N.chart(coord_restrictions=lambda u: [u > -1, u < 1])
            sage: Phi = N.continuous_map(M, {(CN,CM): [u, u^2]}, name='Phi')
            sage: Phi.image()
            Image of the Continuous map Phi
              from the 1-dimensional topological submanifold N
                immersed in the 2-dimensional topological manifold M
              to the 2-dimensional topological manifold M

            sage: S = N.subset('S')
            sage: Phi_S = Phi.image(S); Phi_S
            Image of the Subset S of the
             1-dimensional topological submanifold N
              immersed in the 2-dimensional topological manifold M
             under the Continuous map Phi
              from the 1-dimensional topological submanifold N
               immersed in the 2-dimensional topological manifold M
             to the 2-dimensional topological manifold M
            sage: Phi_S.is_subset(M)
            True
        """
        from .continuous_map_image import ImageManifoldSubset
        if self._is_identity:
            if subset is None:
                return self.domain()
            else:
                return subset
        return ImageManifoldSubset(self, inverse=inverse, domain_subset=subset)

    def preimage(self, codomain_subset, name=None, latex_name=None):
        r"""
        Return the preimage of ``codomain_subset`` under ``self``.

        An alias is :meth:`pullback`.

        INPUT:

        - ``codomain_subset`` -- an instance of
          :class:`~sage.manifolds.subset.ManifoldSubset`
        - ``name`` -- string; name (symbol) given to the subset
        - ``latex_name`` -- string (default: ``None``); LaTeX symbol to
          denote the subset; if none are provided, it is set to ``name``

        OUTPUT:

        - either a :class:`~sage.manifolds.manifold.TopologicalManifold` or
          a :class:`~sage.manifolds.subsets.pullback.ManifoldSubsetPullback`

        EXAMPLES::

            sage: R = Manifold(1, 'R', structure='topological')  # field R
            sage: T.<t> = R.chart()  # canonical chart on R
            sage: R2 = Manifold(2, 'R^2', structure='topological')  # R^2
            sage: c_xy.<x,y> = R2.chart() # Cartesian coordinates on R^2
            sage: Phi = R.continuous_map(R2, [cos(t), sin(t)], name='Phi'); Phi
            Continuous map Phi
             from the 1-dimensional topological manifold R
             to the 2-dimensional topological manifold R^2
            sage: Q1 = R2.open_subset('Q1', coord_def={c_xy: [x>0, y>0]}); Q1
            Open subset Q1 of the 2-dimensional topological manifold R^2
            sage: Phi_inv_Q1 = Phi.preimage(Q1); Phi_inv_Q1
            Subset Phi_inv_Q1 of the 1-dimensional topological manifold R
            sage: R.point([pi/4]) in Phi_inv_Q1
            True
            sage: R.point([0]) in Phi_inv_Q1
            False
            sage: R.point([3*pi/4]) in Phi_inv_Q1
            False

        The identity map is handled specially::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: M.identity_map().preimage(M)
            2-dimensional topological manifold M
            sage: M.identity_map().preimage(M) is M
            True

        Another trivial case::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: D1 = M.open_subset('D1', coord_def={X: x^2+y^2<1}) # the open unit disk
            sage: D2 = M.open_subset('D2', coord_def={X: x^2+y^2<4})
            sage: f = Hom(D1,D2)({(X.restrict(D1), X.restrict(D2)): (2*x, 2*y)}, name='f')
            sage: f.preimage(D2)
            Open subset D1 of the 2-dimensional topological manifold M
            sage: f.preimage(M)
            Open subset D1 of the 2-dimensional topological manifold M
        """
        if self._is_identity:
            return codomain_subset
        if self._codomain.is_subset(codomain_subset):
            return self._domain
        from sage.manifolds.subsets.pullback import ManifoldSubsetPullback
        return ManifoldSubsetPullback(self, codomain_subset,
                                      name=name, latex_name=latex_name)

    pullback = preimage

    #
    # Monoid methods
    #

    def _mul_(self, other):
        r"""
        Composition of ``self`` with another morphism (endomorphism case).

        This applies only when the parent of ``self`` is a monoid, i.e. when
        ``self`` is an endomorphism of the category of topological manifolds,
        i.e. a continuous map `M \to M`, where `M` is a topological manifold.

        INPUT:

        - ``other`` -- a continuous map whose codomain is the domain
          of ``self``

        OUTPUT:

        - :class:`~sage.manifolds.continuous_map.ContinuousMap` that
          is the composite map ``self * other``

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.continuous_map(M, [x+y, x*y], name='f')
            sage: g = M.continuous_map(M, [1-y, 2+x], name='g')
            sage: s = f._mul_(g); s
            Continuous map from the 2-dimensional topological manifold M
             to itself
            sage: s.display()
            M → M
               (x, y) ↦ (x - y + 3, -(x + 2)*y + x + 2)
            sage: s == f*g
            True
            sage: f._mul_(M.identity_map()) == f
            True
            sage: M.identity_map()._mul_(f) == f
            True
        """
        dom = self._domain
        return self._composition_(other, Hom(dom, dom))

    #
    # Other methods
    #

    def _init_derived(self):
        r"""
        Initialize the derived quantities of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.homeomorphism(M, [x+y, x-y])
            sage: f._init_derived()
            sage: f._restrictions
            {}
            sage: f._inverse

        ``_extensions_graph`` and ``_restrictions_graph`` were not originally
        derived quantities, but this induced a bug when dealing with other
        derived quantities (see :issue:`26012`)::

            sage: M = Manifold(2, 'M')
            sage: C.<x, y> = M.chart()
            sage: U = M.open_subset('U', coord_def={C:x>0})
            sage: g = M.metric('g')
            sage: g[:] = [[1, 0], [0, 1]]
            sage: gU = g.restrict(U)
            sage: g[:] = [[1, 0], [0, 2]]
            sage: g.inverse()[:]
            [  1   0]
            [  0 1/2]
            sage: g.inverse().restrict(U)[:] # used to be wrong
            [  1   0]
            [  0 1/2]
        """
        self._restrictions = {} # dict. of restrictions to subdomains of
                                # self._domain
        self._restrictions_graph = {(self._domain, self._codomain): self}
        # dict. of known extensions of self on bigger domains,
        # including self, with pairs of domain codomain as keys.
        # Its elements can be seen as incoming edges on a graph.
        self._extensions_graph = {(self._domain, self._codomain): self}
        # dict. of known restrictions of self on smaller domains,
        # including self, with pairs of domain codomain as keys.
        # Its elements can be seen as outgoing edges on a graph.

        if self._is_identity:
            self._inverse = self
        else:
            self._inverse = None

    def _del_derived(self):
        r"""
        Delete the derived quantities of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.homeomorphism(M, [x+y, x-y])
            sage: f^(-1)
            Homeomorphism of the 2-dimensional topological manifold M
            sage: f._inverse  # was set by f^(-1)
            Homeomorphism of the 2-dimensional topological manifold M
            sage: f._del_derived()
            sage: f._inverse  # has been set to None by _del_derived()
        """
        self._restrictions.clear()
        self._restrictions_graph = {(self._domain, self._codomain): self}
        self._extensions_graph = {(self._domain, self._codomain): self}
        if not self._is_identity:
            self._inverse = None

    def display(self, chart1=None, chart2=None):
        r"""
        Display the expression of ``self`` in one or more pair of charts.

        If the expression is not known already, it is computed from some
        expression in other charts by means of change-of-coordinate formulas.

        INPUT:

        - ``chart1`` -- (default: ``None``) chart on the domain of ``self``;
          if ``None``, the display is performed on all the charts on the
          domain in which the map is known or computable via some change
          of coordinates
        - ``chart2`` -- (default: ``None``) chart on the codomain of ``self``;
          if ``None``, the display is performed on all the charts on the
          codomain in which the map is known or computable via some change
          of coordinates

        The output is either text-formatted (console mode) or LaTeX-formatted
        (notebook mode).

        EXAMPLES:

        A simple reparametrization::

            sage: R.<t> = manifolds.RealLine()
            sage: I = R.open_interval(0, 2*pi)
            sage: J = R.open_interval(2*pi, 6*pi)
            sage: h = J.continuous_map(I, ((t-2*pi)/2,), name='h')
            sage: h.display()
            h: (2*pi, 6*pi) → (0, 2*pi)
               t ↦ t = -pi + 1/2*t
            sage: latex(h.display())
            \begin{array}{llcl} h:& \left(2 \, \pi, 6 \, \pi\right) &
             \longrightarrow & \left(0, 2 \, \pi\right) \\ & t & \longmapsto &
             t = -\pi + \frac{1}{2} \, t \end{array}

        Standard embedding of the sphere `S^2` in `\RR^3`::

            sage: M = Manifold(2, 'S^2', structure='topological') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: N = Manifold(3, 'R^3', latex_name=r'\RR^3', structure='topological')  # R^3
            sage: c_cart.<X,Y,Z> = N.chart()  # Cartesian coordinates on R^3
            sage: Phi = M.continuous_map(N,
            ....:   {(c_xy, c_cart): [2*x/(1+x^2+y^2), 2*y/(1+x^2+y^2), (x^2+y^2-1)/(1+x^2+y^2)],
            ....:    (c_uv, c_cart): [2*u/(1+u^2+v^2), 2*v/(1+u^2+v^2), (1-u^2-v^2)/(1+u^2+v^2)]},
            ....:   name='Phi', latex_name=r'\Phi')
            sage: Phi.display(c_xy, c_cart)
            Phi: S^2 → R^3
            on U: (x, y) ↦ (X, Y, Z) = (2*x/(x^2 + y^2 + 1), 2*y/(x^2 + y^2 + 1), (x^2 + y^2 - 1)/(x^2 + y^2 + 1))
            sage: Phi.display(c_uv, c_cart)
            Phi: S^2 → R^3
            on V: (u, v) ↦ (X, Y, Z) = (2*u/(u^2 + v^2 + 1), 2*v/(u^2 + v^2 + 1), -(u^2 + v^2 - 1)/(u^2 + v^2 + 1))

        The LaTeX output of that embedding is::

            sage: latex(Phi.display(c_xy, c_cart))
            \begin{array}{llcl} \Phi:& S^2 & \longrightarrow & \RR^3
             \\ \text{on}\ U : & \left(x, y\right) & \longmapsto
             & \left(X, Y, Z\right) = \left(\frac{2 \, x}{x^{2} + y^{2} + 1},
               \frac{2 \, y}{x^{2} + y^{2} + 1},
               \frac{x^{2} + y^{2} - 1}{x^{2} + y^{2} + 1}\right)
             \end{array}

        If the argument ``chart2`` is not specified, the display is performed
        on all the charts on the codomain in which the map is known
        or computable via some change of coordinates (here only one chart:
        ``c_cart``)::

            sage: Phi.display(c_xy)
            Phi: S^2 → R^3
            on U: (x, y) ↦ (X, Y, Z) = (2*x/(x^2 + y^2 + 1), 2*y/(x^2 + y^2 + 1), (x^2 + y^2 - 1)/(x^2 + y^2 + 1))

        Similarly, if the argument ``chart1`` is omitted, the display is
        performed on all the charts on the domain of ``Phi`` in which the
        map is known or computable via some change of coordinates::

            sage: Phi.display(chart2=c_cart)
            Phi: S^2 → R^3
            on U: (x, y) ↦ (X, Y, Z) = (2*x/(x^2 + y^2 + 1), 2*y/(x^2 + y^2 + 1), (x^2 + y^2 - 1)/(x^2 + y^2 + 1))
            on V: (u, v) ↦ (X, Y, Z) = (2*u/(u^2 + v^2 + 1), 2*v/(u^2 + v^2 + 1), -(u^2 + v^2 - 1)/(u^2 + v^2 + 1))

        If neither ``chart1`` nor ``chart2`` is specified, the display is
        performed on all the pair of charts in which ``Phi`` is known or
        computable via some change of coordinates::

            sage: Phi.display()
            Phi: S^2 → R^3
            on U: (x, y) ↦ (X, Y, Z) = (2*x/(x^2 + y^2 + 1), 2*y/(x^2 + y^2 + 1), (x^2 + y^2 - 1)/(x^2 + y^2 + 1))
            on V: (u, v) ↦ (X, Y, Z) = (2*u/(u^2 + v^2 + 1), 2*v/(u^2 + v^2 + 1), -(u^2 + v^2 - 1)/(u^2 + v^2 + 1))

        If a chart covers entirely the map's domain, the mention "on ..."
        is omitted::

            sage: Phi.restrict(U).display()
            Phi: U → R^3
               (x, y) ↦ (X, Y, Z) = (2*x/(x^2 + y^2 + 1), 2*y/(x^2 + y^2 + 1), (x^2 + y^2 - 1)/(x^2 + y^2 + 1))

        A shortcut of ``display()`` is ``disp()``::

            sage: Phi.disp()
            Phi: S^2 → R^3
            on U: (x, y) ↦ (X, Y, Z) = (2*x/(x^2 + y^2 + 1), 2*y/(x^2 + y^2 + 1), (x^2 + y^2 - 1)/(x^2 + y^2 + 1))
            on V: (u, v) ↦ (X, Y, Z) = (2*u/(u^2 + v^2 + 1), 2*v/(u^2 + v^2 + 1), -(u^2 + v^2 - 1)/(u^2 + v^2 + 1))

        Display when SymPy is the symbolic engine::

            sage: M.set_calculus_method('sympy')
            sage: N.set_calculus_method('sympy')
            sage: Phi.display(c_xy, c_cart)
            Phi: S^2 → R^3
            on U: (x, y) ↦ (X, Y, Z) = (2*x/(x**2 + y**2 + 1),
             2*y/(x**2 + y**2 + 1), (x**2 + y**2 - 1)/(x**2 + y**2 + 1))
            sage: latex(Phi.display(c_xy, c_cart))
            \begin{array}{llcl} \Phi:& S^2 & \longrightarrow & \RR^3
             \\ \text{on}\ U : & \left(x, y\right) & \longmapsto
             & \left(X, Y, Z\right) = \left(\frac{2 x}{x^{2} + y^{2} + 1},
               \frac{2 y}{x^{2} + y^{2} + 1},
               \frac{x^{2} + y^{2} - 1}{x^{2} + y^{2} + 1}\right)
             \end{array}
        """
        from sage.misc.latex import latex
        from sage.typeset.unicode_characters import unicode_to, unicode_mapsto
        from sage.tensor.modules.format_utilities import FormattedExpansion

        def _display_expression(self, chart1, chart2, result):
            r"""
            Helper function for :meth:`display`.
            """
            try:
                # get coordinate expression
                coord_func = self.coord_functions(chart1, chart2)
                expression = coord_func.expr()
            except (TypeError, ValueError):
                return
            # if that succeeds, proceed:
            coords1 = chart1[:]
            if len(coords1) == 1:
                coords1 = coords1[0]
            coords2 = chart2[:]
            if len(coords2) == 1:
                coords2 = coords2[0]
            if len(expression) == 1:
                expression = expression[0]
                coord_func = coord_func[0]
            if chart1._domain == self._domain:
                result._txt += '   '
                result._latex += ' & '
            else:
                result._txt += 'on ' + chart1._domain._name + ': '
                result._latex += r'\text{on}\ ' + latex(chart1._domain) + \
                                r': & '
            result._txt += repr(coords1) + ' ' + unicode_mapsto + ' '
            result._latex += latex(coords1) + r'& \longmapsto & '
            if chart2 == chart1:
                result._txt += repr(expression) + '\n'
                result._latex += latex(coord_func) + r'\\'
            else:
                result._txt += repr(coords2) + ' = ' + \
                              repr(expression) + '\n'
                result._latex += latex(coords2) + ' = ' + \
                                latex(coord_func) + r'\\'

        result = FormattedExpansion()
        if self._name is None:
            symbol = ''
        else:
            symbol = self._name + ': '
        result._txt = symbol + self._domain._name + ' ' + unicode_to + ' ' \
                      + self._codomain._name + '\n'
        if self._latex_name is None:
            symbol = ''
        else:
            symbol = self._latex_name + ':'
        result._latex = r'\begin{array}{llcl} ' + symbol + r'&' + \
                       latex(self._domain) + r'& \longrightarrow & ' + \
                       latex(self._codomain) + r'\\'
        if chart1 is None:
            if chart2 is None:
                for ch1 in self._domain._top_charts:
                    for ch2 in self._codomain.atlas():
                        _display_expression(self, ch1, ch2, result)
            else:
                for ch1 in self._domain._top_charts:
                    _display_expression(self, ch1, chart2, result)
        else:
            if chart2 is None:
                for ch2 in self._codomain.atlas():
                    _display_expression(self, chart1, ch2, result)
            else:
                _display_expression(self, chart1, chart2, result)
        result._txt = result._txt[:-1]
        result._latex = result._latex[:-2] + r'\end{array}'
        return result

    disp = display

    def coord_functions(self, chart1=None, chart2=None):
        r"""
        Return the functions of the coordinates representing ``self``
        in a given pair of charts.

        If these functions are not already known, they are computed from
        known ones by means of change-of-chart formulas.

        INPUT:

        - ``chart1`` -- (default: ``None``) chart on the domain of ``self``;
          if ``None``, the domain's default chart is assumed
        - ``chart2`` -- (default: ``None``) chart on the codomain of ``self``;
          if ``None``,  the codomain's default chart is assumed

        OUTPUT:

        - a :class:`~sage.manifolds.chart_func.MultiCoordFunction`
          representing the continuous map in the above two charts

        EXAMPLES:

        Continuous map from a 2-dimensional manifold to a 3-dimensional
        one::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: N = Manifold(3, 'N', structure='topological')
            sage: c_uv.<u,v> = M.chart()
            sage: c_xyz.<x,y,z> = N.chart()
            sage: Phi = M.continuous_map(N, (u*v, u/v, u+v), name='Phi',
            ....:                        latex_name=r'\Phi')
            sage: Phi.display()
            Phi: M → N
               (u, v) ↦ (x, y, z) = (u*v, u/v, u + v)
            sage: Phi.coord_functions(c_uv, c_xyz)
            Coordinate functions (u*v, u/v, u + v) on the Chart (M, (u, v))
            sage: Phi.coord_functions() # equivalent to above since 'uv' and 'xyz' are default charts
            Coordinate functions (u*v, u/v, u + v) on the Chart (M, (u, v))
            sage: type(Phi.coord_functions())
            <class 'sage.manifolds.chart_func.MultiCoordFunction'>

        Coordinate representation in other charts::

            sage: c_UV.<U,V> = M.chart()  # new chart on M
            sage: ch_uv_UV = c_uv.transition_map(c_UV, [u-v, u+v])
            sage: ch_uv_UV.inverse()(U,V)
            (1/2*U + 1/2*V, -1/2*U + 1/2*V)
            sage: c_XYZ.<X,Y,Z> = N.chart() # new chart on N
            sage: ch_xyz_XYZ = c_xyz.transition_map(c_XYZ,
            ....:                                   [2*x-3*y+z, y+z-x, -x+2*y-z])
            sage: ch_xyz_XYZ.inverse()(X,Y,Z)
            (3*X + Y + 4*Z, 2*X + Y + 3*Z, X + Y + Z)
            sage: Phi.coord_functions(c_UV, c_xyz)
            Coordinate functions (-1/4*U^2 + 1/4*V^2, -(U + V)/(U - V), V) on
             the Chart (M, (U, V))
            sage: Phi.coord_functions(c_uv, c_XYZ)
            Coordinate functions (((2*u + 1)*v^2 + u*v - 3*u)/v,
             -((u - 1)*v^2 - u*v - u)/v, -((u + 1)*v^2 + u*v - 2*u)/v) on the
             Chart (M, (u, v))
            sage: Phi.coord_functions(c_UV, c_XYZ)
            Coordinate functions
             (-1/2*(U^3 - (U - 2)*V^2 + V^3 - (U^2 + 2*U + 6)*V - 6*U)/(U - V),
              1/4*(U^3 - (U + 4)*V^2 + V^3 - (U^2 - 4*U + 4)*V - 4*U)/(U - V),
              1/4*(U^3 - (U - 4)*V^2 + V^3 - (U^2 + 4*U + 8)*V - 8*U)/(U - V))
             on the Chart (M, (U, V))

        Coordinate representation with respect to a subchart in the domain::

            sage: A = M.open_subset('A', coord_def={c_uv: u>0})
            sage: Phi.coord_functions(c_uv.restrict(A), c_xyz)
            Coordinate functions (u*v, u/v, u + v) on the Chart (A, (u, v))

        Coordinate representation with respect to a superchart
        in the codomain::

            sage: B = N.open_subset('B', coord_def={c_xyz: x<0})
            sage: c_xyz_B = c_xyz.restrict(B)
            sage: Phi1 = M.continuous_map(B, {(c_uv, c_xyz_B): (u*v, u/v, u+v)})
            sage: Phi1.coord_functions(c_uv, c_xyz_B) # definition charts
            Coordinate functions (u*v, u/v, u + v) on the Chart (M, (u, v))
            sage: Phi1.coord_functions(c_uv, c_xyz) # c_xyz = superchart of c_xyz_B
            Coordinate functions (u*v, u/v, u + v) on the Chart (M, (u, v))

        Coordinate representation with respect to a pair
        ``(subchart, superchart)``::

            sage: Phi1.coord_functions(c_uv.restrict(A), c_xyz)
            Coordinate functions (u*v, u/v, u + v) on the Chart (A, (u, v))

        Same example with SymPy as the symbolic calculus engine::

            sage: M.set_calculus_method('sympy')
            sage: N.set_calculus_method('sympy')
            sage: Phi = M.continuous_map(N, (u*v, u/v, u+v), name='Phi',
            ....:                        latex_name=r'\Phi')
            sage: Phi.coord_functions(c_uv, c_xyz)
            Coordinate functions (u*v, u/v, u + v) on the Chart (M, (u, v))
            sage: Phi.coord_functions(c_UV, c_xyz)
            Coordinate functions (-U**2/4 + V**2/4, (-U - V)/(U - V), V) on the Chart (M, (U, V))
            sage: Phi.coord_functions(c_UV, c_XYZ)
            Coordinate functions ((-U**3 + U**2*V + U*V**2 + 2*U*V + 6*U - V**3
             - 2*V**2 + 6*V)/(2*(U - V)), (U**3/4 - U**2*V/4 - U*V**2/4 + U*V
             - U + V**3/4 - V**2 - V)/(U - V), (U**3 - U**2*V - U*V**2 - 4*U*V
             - 8*U + V**3 + 4*V**2 - 8*V)/(4*(U - V))) on the Chart (M, (U, V))
        """
        dom1 = self._domain
        dom2 = self._codomain
        def_chart1 = dom1._def_chart
        def_chart2 = dom2._def_chart
        if chart1 is None:
            chart1 = def_chart1
        if chart2 is None:
            chart2 = def_chart2
        if (chart1, chart2) not in self._coord_expression:
            # Check whether (chart1, chart2) are (subchart, superchart) of
            # a pair of charts where the expression of self is known:
            for (ochart1, ochart2) in self._coord_expression:
                if chart1 in ochart1._subcharts and ochart2 in chart2._subcharts:
                    coord_functions = self._coord_expression[(ochart1, ochart2)].expr()
                    self._coord_expression[(chart1, chart2)] = \
                                         chart1.multifunction(*coord_functions)
                    return self._coord_expression[(chart1, chart2)]
            # Special case of the identity in a single chart:
            if self._is_identity and chart1 == chart2:
                coord_functions = chart1[:]
                self._coord_expression[(chart1, chart1)] = \
                                         chart1.multifunction(*coord_functions)
                return self._coord_expression[(chart1, chart2)]
            # Some change of coordinates must be performed
            change_start = []
            change_arrival = []
            for (ochart1, ochart2) in self._coord_expression:
                if chart1 == ochart1:
                    change_arrival.append(ochart2)
                if chart2 == ochart2:
                    change_start.append(ochart1)
            # 1/ Trying to make a change of chart only on the codomain:
            # the codomain's default chart is privileged:
            sel_chart2 = None # selected chart2
            if (def_chart2 in change_arrival
                    and (def_chart2, chart2) in dom2._coord_changes):
                sel_chart2 = def_chart2
            else:
                for ochart2 in change_arrival:
                    if (ochart2, chart2) in dom2._coord_changes:
                        sel_chart2 = ochart2
                        break
            if sel_chart2 is not None:
                oexpr = self._coord_expression[(chart1, sel_chart2)]
                chg2 = dom2._coord_changes[(sel_chart2, chart2)]
                self._coord_expression[(chart1, chart2)] = \
                                 chart1.multifunction( *chg2(*oexpr.expr()) )
                return self._coord_expression[(chart1, chart2)]

            # 2/ Trying to make a change of chart only on the start domain:
            # the domain's default chart is privileged:
            sel_chart1 = None # selected chart1
            if (def_chart1 in change_start
                    and (chart1, def_chart1) in dom1._coord_changes):
                sel_chart1 = def_chart1
            else:
                for ochart1 in change_start:
                    if (chart1, ochart1) in dom1._coord_changes:
                        sel_chart1 = ochart1
                        break
            if sel_chart1 is not None:
                oexpr = self._coord_expression[(sel_chart1, chart2)]
                chg1 = dom1._coord_changes[(chart1, sel_chart1)]
                self._coord_expression[(chart1, chart2)] = \
                       chart1.multifunction( *oexpr(*chg1._transf.expr()) )
                return self._coord_expression[(chart1, chart2)]

            # 3/ If this point is reached, it is necessary to perform some
            # coordinate change both on the start domain and the arrival one
            # the default charts are privileged:
            if ((def_chart1, def_chart2) in self._coord_expression
                    and (chart1, def_chart1) in dom1._coord_changes
                    and (def_chart2, chart2) in dom2._coord_changes):
                sel_chart1 = def_chart1
                sel_chart2 = def_chart2
            else:
                for (ochart1, ochart2) in self._coord_expression:
                    if ((chart1, ochart1) in dom1._coord_changes
                            and (ochart2, chart2) in dom2._coord_changes):
                        sel_chart1 = ochart1
                        sel_chart2 = ochart2
                        break
            if sel_chart1 is not None and sel_chart2 is not None:
                oexpr = self._coord_expression[(sel_chart1, sel_chart2)]
                chg1 = dom1._coord_changes[(chart1, sel_chart1)]
                chg2 = dom2._coord_changes[(sel_chart2, chart2)]
                self._coord_expression[(chart1, chart2)] = chart1.multifunction(
                                  *chg2( *oexpr(*chg1._transf.expr()) ) )
                return self._coord_expression[(chart1, chart2)]

            # 4/ If this point is reached, the demanded value cannot be
            # computed
            raise ValueError("the expression of the map in the pair " +
                             "({}, {})".format(chart1, chart2) + " cannot " +
                             "be computed by means of known changes of charts")

        return self._coord_expression[(chart1, chart2)]

    def expr(self, chart1=None, chart2=None):
        r"""
        Return the expression of ``self`` in terms of
        specified coordinates.

        If the expression is not already known, it is computed from some
        known expression by means of change-of-chart formulas.

        INPUT:

        - ``chart1`` -- (default: ``None``) chart on the map's domain;
          if ``None``, the domain's default chart is assumed
        - ``chart2`` -- (default: ``None``) chart on the map's codomain;
          if ``None``, the codomain's default chart is assumed

        OUTPUT:

        - symbolic expression representing the continuous map in the
          above two charts

        EXAMPLES:

        Continuous map from a 2-dimensional manifold to a
        3-dimensional one::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: N = Manifold(3, 'N', structure='topological')
            sage: c_uv.<u,v> = M.chart()
            sage: c_xyz.<x,y,z> = N.chart()
            sage: Phi = M.continuous_map(N, (u*v, u/v, u+v), name='Phi',
            ....:                        latex_name=r'\Phi')
            sage: Phi.display()
            Phi: M → N
               (u, v) ↦ (x, y, z) = (u*v, u/v, u + v)
            sage: Phi.expr(c_uv, c_xyz)
            (u*v, u/v, u + v)
            sage: Phi.expr()  # equivalent to above since 'uv' and 'xyz' are default charts
            (u*v, u/v, u + v)
            sage: type(Phi.expr()[0])
            <class 'sage.symbolic.expression.Expression'>

        Expressions in other charts::

            sage: c_UV.<U,V> = M.chart()  # new chart on M
            sage: ch_uv_UV = c_uv.transition_map(c_UV, [u-v, u+v])
            sage: ch_uv_UV.inverse()(U,V)
            (1/2*U + 1/2*V, -1/2*U + 1/2*V)
            sage: c_XYZ.<X,Y,Z> = N.chart() # new chart on N
            sage: ch_xyz_XYZ = c_xyz.transition_map(c_XYZ,
            ....:                                   [2*x-3*y+z, y+z-x, -x+2*y-z])
            sage: ch_xyz_XYZ.inverse()(X,Y,Z)
            (3*X + Y + 4*Z, 2*X + Y + 3*Z, X + Y + Z)
            sage: Phi.expr(c_UV, c_xyz)
            (-1/4*U^2 + 1/4*V^2, -(U + V)/(U - V), V)
            sage: Phi.expr(c_uv, c_XYZ)
            (((2*u + 1)*v^2 + u*v - 3*u)/v,
             -((u - 1)*v^2 - u*v - u)/v,
             -((u + 1)*v^2 + u*v - 2*u)/v)
            sage: Phi.expr(c_UV, c_XYZ)
             (-1/2*(U^3 - (U - 2)*V^2 + V^3 - (U^2 + 2*U + 6)*V - 6*U)/(U - V),
              1/4*(U^3 - (U + 4)*V^2 + V^3 - (U^2 - 4*U + 4)*V - 4*U)/(U - V),
              1/4*(U^3 - (U - 4)*V^2 + V^3 - (U^2 + 4*U + 8)*V - 8*U)/(U - V))

        A rotation in some Euclidean plane::

            sage: M = Manifold(2, 'M', structure='topological') # the plane (minus a segment to have global regular spherical coordinates)
            sage: c_spher.<r,ph> = M.chart(r'r:(0,+oo) ph:(0,2*pi):\phi') # spherical coordinates on the plane
            sage: rot = M.continuous_map(M, (r, ph+pi/3), name='R') # pi/3 rotation around r=0
            sage: rot.expr()
            (r, 1/3*pi + ph)

        Expression of the rotation in terms of Cartesian coordinates::

            sage: c_cart.<x,y> = M.chart() # Declaration of Cartesian coordinates
            sage: ch_spher_cart = c_spher.transition_map(c_cart,
            ....:                 [r*cos(ph), r*sin(ph)]) # relation to spherical coordinates
            sage: ch_spher_cart.set_inverse(sqrt(x^2+y^2), atan2(y,x))
            Check of the inverse coordinate transformation:
              r == r  *passed*
              ph == arctan2(r*sin(ph), r*cos(ph))  **failed**
              x == x  *passed*
              y == y  *passed*
            NB: a failed report can reflect a mere lack of simplification.
            sage: rot.expr(c_cart, c_cart)
            (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y)
        """
        return self.coord_functions(chart1, chart2).expr()

    expression = expr

    def set_expr(self, chart1, chart2, coord_functions):
        r"""
        Set a new coordinate representation of ``self``.

        The expressions with respect to other charts are deleted, in order to
        avoid any inconsistency. To keep them, use :meth:`add_expr` instead.

        INPUT:

        - ``chart1`` -- chart for the coordinates on the domain of ``self``
        - ``chart2`` -- chart for the coordinates on the codomain of ``self``
        - ``coord_functions`` -- the coordinate symbolic expression of the
          map in the above charts: list (or tuple) of the coordinates of
          the image expressed in terms of the coordinates of the considered
          point; if the dimension of the arrival manifold is 1, a single
          coordinate expression can be passed instead of a tuple with a
          single element

        EXAMPLES:

        Polar representation of a planar rotation initially defined in
        Cartesian coordinates::

            sage: M = Manifold(2, 'R^2', latex_name=r'\RR^2', structure='topological')  # the Euclidean plane R^2
            sage: c_xy.<x,y> = M.chart() # Cartesian coordinate on R^2
            sage: U = M.open_subset('U', coord_def={c_xy: (y!=0, x<0)}) # the complement of the segment y=0 and x>0
            sage: c_cart = c_xy.restrict(U) # Cartesian coordinates on U
            sage: c_spher.<r,ph> = U.chart(r'r:(0,+oo) ph:(0,2*pi):\phi') # spherical coordinates on U

        Links between spherical coordinates and Cartesian ones::

            sage: ch_cart_spher = c_cart.transition_map(c_spher,
            ....:                                       [sqrt(x*x+y*y), atan2(y,x)])
            sage: ch_cart_spher.set_inverse(r*cos(ph), r*sin(ph))
            Check of the inverse coordinate transformation:
              x == x  *passed*
              y == y  *passed*
              r == r  *passed*
              ph == arctan2(r*sin(ph), r*cos(ph))  **failed**
            NB: a failed report can reflect a mere lack of simplification.
            sage: rot = U.continuous_map(U, ((x - sqrt(3)*y)/2, (sqrt(3)*x + y)/2),
            ....:                        name='R')
            sage: rot.display(c_cart, c_cart)
            R: U → U
               (x, y) ↦ (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y)

        Let us use the method :meth:`set_expr` to set the
        spherical-coordinate expression by hand::

            sage: rot.set_expr(c_spher, c_spher, (r, ph+pi/3))
            sage: rot.display(c_spher, c_spher)
            R: U → U
               (r, ph) ↦ (r, 1/3*pi + ph)

        The expression in Cartesian coordinates has been erased::

            sage: rot._coord_expression
            {(Chart (U, (r, ph)),
              Chart (U, (r, ph))): Coordinate functions (r, 1/3*pi + ph)
              on the Chart (U, (r, ph))}

        It is recovered (thanks to the known change of coordinates) by a call
        to :meth:`display`::

            sage: rot.display(c_cart, c_cart)
            R: U → U
               (x, y) ↦ (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y)

            sage: rot._coord_expression  # random (dictionary output)
            {(Chart (U, (x, y)),
              Chart (U, (x, y))): Coordinate functions (-1/2*sqrt(3)*y + 1/2*x,
              1/2*sqrt(3)*x + 1/2*y) on the Chart (U, (x, y)),
             (Chart (U, (r, ph)),
              Chart (U, (r, ph))): Coordinate functions (r, 1/3*pi + ph)
              on the Chart (U, (r, ph))}

        TESTS:

        We check that this does not change the equality nor the hash value::

            sage: M = Manifold(2, 'R^2', latex_name=r'\RR^2', structure='topological')
            sage: c_xy.<x,y> = M.chart()
            sage: U = M.open_subset('U', coord_def={c_xy: (y!=0, x<0)})
            sage: c_cart = c_xy.restrict(U)
            sage: c_spher.<r,ph> = U.chart(r'r:(0,+oo) ph:(0,2*pi):\phi')
            sage: ch_cart_spher = c_cart.transition_map(c_spher,
            ....:                                       [sqrt(x*x+y*y), atan2(y,x)])
            sage: ch_cart_spher.set_inverse(r*cos(ph), r*sin(ph))
            Check of the inverse coordinate transformation:
              x == x  *passed*
              y == y  *passed*
              r == r  *passed*
              ph == arctan2(r*sin(ph), r*cos(ph))  **failed**
            NB: a failed report can reflect a mere lack of simplification.
            sage: rot = U.continuous_map(U, ((x - sqrt(3)*y)/2, (sqrt(3)*x + y)/2),
            ....:                        name='R')
            sage: rot2 = copy(rot)
            sage: rot == rot2 and hash(rot) == hash(rot2)
            True
            sage: rot.set_expr(c_spher, c_spher, (r, ph+pi/3))
            sage: rot == rot2 and hash(rot) == hash(rot2)
            True
        """
        if self._is_identity:
            raise NotImplementedError("set_expr() must not be used for the identity map")
        if chart1 not in self._domain.atlas():
            raise ValueError("the {}".format(chart1) +
                             " has not been defined on the {}".format(self._domain))
        if chart2 not in self._codomain.atlas():
            raise ValueError("the {}".format(chart2) +
                             " has not been defined on the {}".format(self._codomain))
        self._coord_expression.clear()
        self._del_derived()
        n2 = self._codomain.dim()
        if n2 > 1:
            if len(coord_functions) != n2:
                raise ValueError("{} coordinate functions must ".format(n2) +
                                 "be provided.")
            self._coord_expression[(chart1, chart2)] = \
                                         chart1.multifunction(*coord_functions)
        else:
            if isinstance(coord_functions, (list, tuple)):
                coord_functions = coord_functions[0]
            self._coord_expression[(chart1, chart2)] = \
                                          chart1.multifunction(coord_functions)

    set_expression = set_expr

    def add_expr(self, chart1, chart2, coord_functions):
        r"""
        Set a new coordinate representation of ``self``.

        The previous expressions with respect to other charts are kept. To
        clear them, use :meth:`set_expr` instead.

        INPUT:

        - ``chart1`` -- chart for the coordinates on the map's domain
        - ``chart2`` -- chart for the coordinates on the map's codomain
        - ``coord_functions`` -- the coordinate symbolic expression of the
          map in the above charts: list (or tuple) of the coordinates of
          the image expressed in terms of the coordinates of the considered
          point; if the dimension of the arrival manifold is 1, a single
          coordinate expression can be passed instead of a tuple with a
          single element

        .. WARNING::

            If the map has already expressions in other charts, it
            is the user's responsibility to make sure that the expression
            to be added is consistent with them.

        EXAMPLES:

        Polar representation of a planar rotation initially defined in
        Cartesian coordinates::

            sage: M = Manifold(2, 'R^2', latex_name=r'\RR^2', structure='topological')  # the Euclidean plane R^2
            sage: c_xy.<x,y> = M.chart() # Cartesian coordinate on R^2
            sage: U = M.open_subset('U', coord_def={c_xy: (y!=0, x<0)}) # the complement of the segment y=0 and x>0
            sage: c_cart = c_xy.restrict(U) # Cartesian coordinates on U
            sage: c_spher.<r,ph> = U.chart(r'r:(0,+oo) ph:(0,2*pi):\phi') # spherical coordinates on U

        We construct the links between spherical coordinates and
        Cartesian ones::

            sage: ch_cart_spher = c_cart.transition_map(c_spher, [sqrt(x*x+y*y), atan2(y,x)])
            sage: ch_cart_spher.set_inverse(r*cos(ph), r*sin(ph))
            Check of the inverse coordinate transformation:
              x == x  *passed*
              y == y  *passed*
              r == r  *passed*
              ph == arctan2(r*sin(ph), r*cos(ph))  **failed**
            NB: a failed report can reflect a mere lack of simplification.
            sage: rot = U.continuous_map(U, ((x - sqrt(3)*y)/2, (sqrt(3)*x + y)/2),
            ....:                        name='R')
            sage: rot.display(c_cart, c_cart)
            R: U → U
               (x, y) ↦ (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y)

        If we calculate the expression in terms of spherical coordinates,
        via the method :meth:`display`, we notice some difficulties
        in ``arctan2`` simplifications::

            sage: rot.display(c_spher, c_spher)
            R: U → U
               (r, ph) ↦ (r, arctan2(1/2*(sqrt(3)*cos(ph) + sin(ph))*r, -1/2*(sqrt(3)*sin(ph) - cos(ph))*r))

        Therefore, we use the method :meth:`add_expr` to set the
        spherical-coordinate expression by hand::

            sage: rot.add_expr(c_spher, c_spher, (r, ph+pi/3))
            sage: rot.display(c_spher, c_spher)
            R: U → U
               (r, ph) ↦ (r, 1/3*pi + ph)

        The call to :meth:`add_expr` has not deleted the expression in
        terms of Cartesian coordinates, as we can check by printing the
        internal dictionary ``_coord_expression``, which stores the
        various internal representations of the continuous map::

            sage: rot._coord_expression # random (dictionary output)
            {(Chart (U, (x, y)), Chart (U, (x, y))):
             Coordinate functions (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y)
              on the Chart (U, (x, y)),
             (Chart (U, (r, ph)), Chart (U, (r, ph))):
              Coordinate functions (r, 1/3*pi + ph) on the Chart (U, (r, ph))}

        If, on the contrary, we use :meth:`set_expr`, the expression in
        Cartesian coordinates is lost::

            sage: rot.set_expr(c_spher, c_spher, (r, ph+pi/3))
            sage: rot._coord_expression
            {(Chart (U, (r, ph)), Chart (U, (r, ph))):
             Coordinate functions (r, 1/3*pi + ph) on the Chart (U, (r, ph))}

        It is recovered (thanks to the known change of coordinates) by
        a call to :meth:`display`::

            sage: rot.display(c_cart, c_cart)
            R: U → U
               (x, y) ↦ (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y)

            sage: rot._coord_expression # random (dictionary output)
            {(Chart (U, (x, y)), Chart (U, (x, y))):
             Coordinate functions (-1/2*sqrt(3)*y + 1/2*x, 1/2*sqrt(3)*x + 1/2*y)
              on the Chart (U, (x, y)),
             (Chart (U, (r, ph)), Chart (U, (r, ph))):
             Coordinate functions (r, 1/3*pi + ph) on the Chart (U, (r, ph))}

        The rotation can be applied to a point by means of either
        coordinate system::

            sage: p = M.point((1,2))  #  p defined by its Cartesian coord.
            sage: q = rot(p)  # q is computed by means of Cartesian coord.
            sage: p1 = M.point((sqrt(5), arctan(2)), chart=c_spher) # p1 is defined only in terms of c_spher
            sage: q1 = rot(p1) # computation by means of spherical coordinates
            sage: q1 == q
            True
        """
        if self._is_identity:
            raise NotImplementedError("add_expr() must not be used for the identity map")
        if chart1 not in self._domain.atlas():
            raise ValueError("the {}".format(chart1) +
                             " has not been defined on the {}".format(self._domain))
        if chart2 not in self._codomain.atlas():
            raise ValueError("the {}".format(chart2) +
                             " has not been defined on the {}".format(self._codomain))
        self._del_derived()
        n2 = self._codomain.dim()
        if n2 > 1:
            if len(coord_functions) != n2:
                raise ValueError("{} coordinate functions must be provided".format(n2))
            self._coord_expression[(chart1, chart2)] = chart1.multifunction(*coord_functions)
        else:
            if isinstance(coord_functions, (list, tuple)):
                coord_functions = coord_functions[0]
            self._coord_expression[(chart1, chart2)] = chart1.multifunction(coord_functions)

    add_expression = add_expr

    def restrict(self, subdomain, subcodomain=None):
        r"""
        Restriction of ``self`` to some open subset of its
        domain of definition.

        INPUT:

        - ``subdomain`` -- :class:`~sage.manifolds.manifold.TopologicalManifold`;
          an open subset of the domain of ``self``
        - ``subcodomain`` -- (default: ``None``) an open subset of the codomain
          of ``self``; if ``None``, the codomain of ``self`` is assumed

        OUTPUT:

        - a :class:`ContinuousMap` that is the restriction
          of ``self`` to ``subdomain``

        EXAMPLES:

        Restriction to an annulus of a homeomorphism between the open unit
        disk and `\RR^2`::

            sage: M = Manifold(2, 'R^2', structure='topological')  # R^2
            sage: c_xy.<x,y> = M.chart()  # Cartesian coord. on R^2
            sage: D = M.open_subset('D', coord_def={c_xy: x^2+y^2<1}) # the open unit disk
            sage: Phi = D.continuous_map(M, [x/sqrt(1-x^2-y^2), y/sqrt(1-x^2-y^2)],
            ....:                        name='Phi', latex_name=r'\Phi')
            sage: Phi.display()
            Phi: D → R^2
               (x, y) ↦ (x, y) = (x/sqrt(-x^2 - y^2 + 1), y/sqrt(-x^2 - y^2 + 1))
            sage: c_xy_D = c_xy.restrict(D)
            sage: U = D.open_subset('U', coord_def={c_xy_D: x^2+y^2>1/2}) # the annulus 1/2 < r < 1
            sage: Phi.restrict(U)
            Continuous map Phi
             from the Open subset U of the 2-dimensional topological manifold R^2
             to the 2-dimensional topological manifold R^2
            sage: Phi.restrict(U).parent()
            Set of Morphisms from Open subset U of the 2-dimensional topological
             manifold R^2 to 2-dimensional topological manifold R^2 in Category
             of manifolds over Real Field with 53 bits of precision
            sage: Phi.domain()
            Open subset D of the 2-dimensional topological manifold R^2
            sage: Phi.restrict(U).domain()
            Open subset U of the 2-dimensional topological manifold R^2
            sage: Phi.restrict(U).display()
            Phi: U → R^2
               (x, y) ↦ (x, y) = (x/sqrt(-x^2 - y^2 + 1), y/sqrt(-x^2 - y^2 + 1))

        The result is cached::

            sage: Phi.restrict(U) is Phi.restrict(U)
            True

        The restriction of the identity map::

            sage: id = D.identity_map() ; id
            Identity map Id_D of the Open subset D of the 2-dimensional
             topological manifold R^2
            sage: id.restrict(U)
            Identity map Id_U of the Open subset U of the 2-dimensional
             topological manifold R^2
            sage: id.restrict(U) is U.identity_map()
            True

        The codomain can be restricted (i.e. made tighter)::

            sage: Phi = D.continuous_map(M, [x/sqrt(1+x^2+y^2), y/sqrt(1+x^2+y^2)])
            sage: Phi
            Continuous map from
             the Open subset D of the 2-dimensional topological manifold R^2
             to the 2-dimensional topological manifold R^2
            sage: Phi.restrict(D, subcodomain=D)
            Continuous map from the Open subset D of the 2-dimensional
             topological manifold R^2 to itself
        """
        if subcodomain is None:
            if self._is_identity:
                subcodomain = subdomain
            else:
                subcodomain = self._codomain
        if subdomain == self._domain and subcodomain == self._codomain:
            return self
        if (subdomain, subcodomain) not in self._restrictions:
            if not subdomain.is_subset(self._domain):
                raise ValueError("the specified domain is not a subset"
                                 " of the domain of definition of the"
                                 " continuous map")
            if not subcodomain.is_subset(self._codomain):
                raise ValueError("the specified codomain is not a subset"
                                 " of the codomain of the continuous map")
            # Special case of the identity map:
            if self._is_identity:
                self._restrictions[(subdomain, subcodomain)] = subdomain.identity_map()
                return self._restrictions[(subdomain, subcodomain)]

            # First one tries to get the restriction from a tighter domain:
            for dom, rst in self._restrictions.items():
                if subdomain.is_subset(dom[0]) and (subdomain, subcodomain) in rst._restrictions:
                    res = rst._restrictions[(subdomain, subcodomain)]
                    self._restrictions[(subdomain, subcodomain)] = res
                    self._restrictions.update(res._restrictions)
                    self._restrictions_graph.update(res._restrictions_graph)
                    res._extensions_graph.update(self._extensions_graph)
                    for ext in self._extensions_graph.values():
                        ext._restrictions[subdomain] = res
                        ext._restrictions.update(res._restrictions)
                        ext._restrictions_graph.update(res._restrictions_graph)
                    return self._restrictions[(subdomain, subcodomain)]

            # Maybe it didn't exist but could have:
            for dom, rst in self._restrictions.items():
                if subdomain.is_subset(dom[0]) and subcodomain.is_subset(dom[1]):
                    res = rst.restrict(subdomain,subcodomain) # all propagation
                                                              # is done here
                    # should be useless:
                    self._restrictions[(subdomain,subcodomain)] = res
                    self._restrictions_graph[(subdomain, subcodomain)] = res
                    return self._restrictions[(subdomain,subcodomain)]

            # Secondly one tries to get the restriction from one previously
            # defined on a larger domain:
            for dom, ext in self._extensions_graph.items():
                if (subdomain,subcodomain) in ext._restrictions:
                    res = ext._restrictions[(subdomain,subcodomain)]
                    self._restrictions[(subdomain, subcodomain)] = res
                    self._restrictions.update(res._restrictions)
                    self._restrictions_graph.update(res._restrictions_graph)
                    res._extensions_graph.update(self._extensions_graph)
                    for ext in self._extensions_graph.values():
                        ext._restrictions[subdomain] = res
                        ext._restrictions.update(res._restrictions)
                        ext._restrictions_graph.update(res._restrictions_graph)
                    return self._restrictions[(subdomain, subcodomain)]

            # Generic case:
            homset = Hom(subdomain, subcodomain)
            resu = type(self)(homset, name=self._name,
                              latex_name=self._latex_name)
            for charts in self._coord_expression:
                for ch1 in charts[0]._subcharts:
                    if ch1._domain.is_subset(subdomain):
                        for ch2 in charts[1]._subcharts:
                            if ch2._domain.is_subset(subcodomain):
                                for sch2 in ch2._supercharts:
                                    if (ch1, sch2) in resu._coord_expression:
                                        break
                                else:
                                    for sch2 in ch2._subcharts:
                                        if (ch1, sch2) in resu._coord_expression:
                                            del resu._coord_expression[(ch1, sch2)]
                                    coord_functions = self._coord_expression[charts].expr()
                                    resu._coord_expression[(ch1, ch2)] = \
                                            ch1.multifunction(*coord_functions)

            # propagate extensions
            for dom, ext in self._extensions_graph.items(): # includes self
                ext._restrictions[(subdomain, subcodomain)] = resu
                ext._restrictions_graph[(subdomain, subcodomain)] = resu

            # propagate restrictions
            for dom, rst in self._restrictions.items():
                if dom[0].is_subset(subdomain) and dom[1].is_subset(subcodomain):
                    if rst is not resu:
                        resu._restrictions.update(rst._restrictions_graph)
                        resu._restrictions_graph.update(rst._restrictions_graph)
                        rst._extensions_graph.update(resu._extensions_graph)

            self._restrictions[(subdomain, subcodomain)] = resu
            self._restrictions_graph[(subdomain, subcodomain)] = resu
            resu._extensions_graph.update(self._extensions_graph)

        return self._restrictions[(subdomain, subcodomain)]

    def __invert__(self):
        r"""
        Return the inverse of ``self`` if it is an isomorphism.

        OUTPUT: the inverse isomorphism

        EXAMPLES:

        The inverse of a rotation in the Euclidean plane::

            sage: M = Manifold(2, 'R^2', latex_name=r'\RR^2', structure='topological')
            sage: c_cart.<x,y> = M.chart()

        A pi/3 rotation around the origin::

            sage: rot = M.homeomorphism(M, ((x - sqrt(3)*y)/2, (sqrt(3)*x + y)/2),
            ....:                       name='R')
            sage: rot.inverse()
            Homeomorphism R^(-1) of the 2-dimensional topological manifold R^2
            sage: rot.inverse().display()
            R^(-1): R^2 → R^2
               (x, y) ↦ (1/2*sqrt(3)*y + 1/2*x, -1/2*sqrt(3)*x + 1/2*y)

        Checking that applying successively the homeomorphism and its
        inverse results in the identity::

            sage: (a, b) = var('a b')
            sage: p = M.point((a,b)) # a generic point on M
            sage: q = rot(p)
            sage: p1 = rot.inverse()(q)
            sage: p1 == p
            True

        The result is cached::

            sage: rot.inverse() is rot.inverse()
            True

        The notations ``^(-1)`` or ``~`` can also be used for the inverse::

            sage: rot^(-1) is rot.inverse()
            True
            sage: ~rot is rot.inverse()
            True

        An example with multiple charts: the equatorial symmetry on the
        2-sphere::

            sage: M = Manifold(2, 'M', structure='topological') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                                intersection_name='W',
            ....:                                restrictions1=x^2+y^2!=0,
            ....:                                restrictions2=u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: s = M.homeomorphism(M, {(c_xy, c_uv): [x, y], (c_uv, c_xy): [u, v]},
            ....:                     name='s')
            sage: s.display()
            s: M → M
            on U: (x, y) ↦ (u, v) = (x, y)
            on V: (u, v) ↦ (x, y) = (u, v)
            sage: si = s.inverse(); si
            Homeomorphism s^(-1) of the 2-dimensional topological manifold M
            sage: si.display()
            s^(-1): M → M
            on U: (x, y) ↦ (u, v) = (x, y)
            on V: (u, v) ↦ (x, y) = (u, v)

        The equatorial symmetry is of course an involution::

            sage: si == s
            True
        """
        from sage.symbolic.ring import SR
        from sage.symbolic.relation import solve
        if self._inverse is not None:
            return self._inverse
        if not self._is_isomorphism:
            raise ValueError("the {} is not an isomorphism".format(self))
        coord_functions = {} # coordinate expressions of the result
        for (chart1, chart2) in self._coord_expression:
            coord_map = self._coord_expression[(chart1, chart2)]
            n1 = len(chart1._xx)
            n2 = len(chart2._xx)
            # New symbolic variables (different from chart2._xx to allow for a
            #  correct solution even when chart2 = chart1):
            x2 = SR.temp_var(n=n2)
            equations = [x2[i] == coord_map._functions[i].expr()
                         for i in range(n2)]
            solutions = solve(equations, chart1._xx, solution_dict=True)
            if not solutions:
                raise ValueError("no solution found")
            if len(solutions) > 1:
                raise ValueError("non-unique solution found")
            substitutions = dict(zip(x2, chart2._xx))
            sol = solutions[0]
            inv_functions = [sol[chart1._xx[i]].subs(substitutions)
                             for i in range(n1)]
            for i in range(n1):
                x = inv_functions[i]
                try:
                    # simplify derived from calculus_method
                    inv_functions[i] = chart2.simplify(x)
                except AttributeError:
                    pass
            coord_functions[(chart2, chart1)] = inv_functions
            SR.cleanup_var(x2)
        if self._name is None:
            name = None
        else:
            name = self._name + '^(-1)'
        if self._latex_name is None:
            latex_name = None
        else:
            latex_name = self._latex_name + r'^{-1}'
        homset = Hom(self._codomain, self._domain)
        self._inverse = type(self)(homset, coord_functions=coord_functions,
                                   name=name, latex_name=latex_name,
                                   is_isomorphism=True)
        return self._inverse

    inverse = __invert__
