r"""
Differentiable Maps between Differentiable Manifolds

The class :class:`DiffMap` implements differentiable maps from a differentiable
manifold `M` to a differentiable manifold `N` over the same topological field
`K` as `M` (in most applications, `K = \RR` or `K = \CC`):

.. MATH::

    \Phi: M \longrightarrow N


AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015): initial version
- Marco Mancini (2018): pullback parallelization

REFERENCES:

- Chap. 1 of [KN1963]_
- Chaps. 2 and 3 of [Lee2013]_
"""

# ****************************************************************************
#       Copyright (C) 2015-2021 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#       Copyright (C) 2018 Marco Mancini <marco.mancini@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import annotations

from typing import TYPE_CHECKING

from sage.manifolds.continuous_map import ContinuousMap
from sage.parallel.decorate import parallel
from sage.parallel.parallelism import Parallelism

if TYPE_CHECKING:
    from sage.manifolds.point import ManifoldPoint
    from sage.tensor.modules.free_module_morphism import FiniteRankFreeModuleMorphism


class DiffMap(ContinuousMap):
    r"""
    Differentiable map between two differentiable manifolds.

    This class implements differentiable maps of the type

    .. MATH::

        \Phi: M \longrightarrow  N

    where `M` and `N` are differentiable manifolds over the same topological
    field `K` (in most applications, `K = \RR` or `K = \CC`).

    Differentiable maps are the *morphisms* of the *category* of
    differentiable manifolds. The set of all differentiable maps from `M` to
    `N` is therefore the homset between `M` and `N`, which is denoted by
    `\mathrm{Hom}(M,N)`.

    The class :class:`DiffMap` is a Sage *element* class, whose *parent*
    class is
    :class:`~sage.manifolds.differentiable.manifold_homset.DifferentiableManifoldHomset`.
    It inherits from the class
    :class:`~sage.manifolds.continuous_map.ContinuousMap` since a
    differentiable map is obviously a continuous one.

    INPUT:

    - ``parent`` -- homset `\mathrm{Hom}(M,N)` to which the differentiable
      map belongs
    - ``coord_functions`` -- (default: ``None``) if not ``None``, must be
      a dictionary of the coordinate expressions (as lists (or tuples) of the
      coordinates of the image expressed in terms of the coordinates of
      the considered point) with the pairs of charts (chart1, chart2)
      as keys (chart1 being a chart on `M` and chart2 a chart on `N`).
      If the dimension of the map's codomain is 1, a single coordinate
      expression can be passed instead of a tuple with a single element
    - ``name`` -- (default: ``None``) name given to the differentiable map
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
      differentiable map; if ``None``, the LaTeX symbol is set to ``name``
    - ``is_isomorphism`` -- boolean (default: ``False``); determines whether the
      constructed object is a isomorphism (i.e. a diffeomorphism). If set to
      ``True``, then the manifolds `M` and `N` must have the same dimension.
    - ``is_identity`` -- boolean (default: ``False``); determines whether the
      constructed object is the identity map. If set to ``True``,
      then `N` must be `M` and the entry ``coord_functions`` is not used.

    .. NOTE::

        If the information passed by means of the argument ``coord_functions``
        is not sufficient to fully specify the differentiable map,
        further coordinate expressions, in other charts, can be subsequently
        added by means of the method
        :meth:`~sage.manifolds.continuous_map.ContinuousMap.add_expr`

    EXAMPLES:

    The standard embedding of the sphere `S^2` into `\RR^3`::

        sage: M = Manifold(2, 'S^2') # the 2-dimensional sphere S^2
        sage: U = M.open_subset('U') # complement of the North pole
        sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
        sage: V = M.open_subset('V') # complement of the South pole
        sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
        sage: M.declare_union(U,V)   # S^2 is the union of U and V
        sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
        ....:                 intersection_name='W', restrictions1= x^2+y^2!=0,
        ....:                 restrictions2= u^2+v^2!=0)
        sage: uv_to_xy = xy_to_uv.inverse()
        sage: N = Manifold(3, 'R^3', r'\RR^3')  # R^3
        sage: c_cart.<X,Y,Z> = N.chart()  # Cartesian coordinates on R^3
        sage: Phi = M.diff_map(N,
        ....: {(c_xy, c_cart): [2*x/(1+x^2+y^2), 2*y/(1+x^2+y^2), (x^2+y^2-1)/(1+x^2+y^2)],
        ....:  (c_uv, c_cart): [2*u/(1+u^2+v^2), 2*v/(1+u^2+v^2), (1-u^2-v^2)/(1+u^2+v^2)]},
        ....: name='Phi', latex_name=r'\Phi')
        sage: Phi
        Differentiable map Phi from the 2-dimensional differentiable manifold
         S^2 to the 3-dimensional differentiable manifold R^3
        sage: Phi.parent()
        Set of Morphisms from 2-dimensional differentiable manifold S^2 to
         3-dimensional differentiable manifold R^3 in Category of smooth
         manifolds over Real Field with 53 bits of precision
        sage: Phi.parent() is Hom(M, N)
        True
        sage: type(Phi)
        <class 'sage.manifolds.differentiable.manifold_homset.DifferentiableManifoldHomset_with_category.element_class'>
        sage: Phi.display()
        Phi: S^2 → R^3
        on U: (x, y) ↦ (X, Y, Z) = (2*x/(x^2 + y^2 + 1), 2*y/(x^2 + y^2 + 1),
                                    (x^2 + y^2 - 1)/(x^2 + y^2 + 1))
        on V: (u, v) ↦ (X, Y, Z) = (2*u/(u^2 + v^2 + 1), 2*v/(u^2 + v^2 + 1),
                                    -(u^2 + v^2 - 1)/(u^2 + v^2 + 1))

    It is possible to create the map via the method
    :meth:`~sage.manifolds.differentiable.manifold.DifferentiableManifold.diff_map`
    only in a single pair of charts: the argument ``coord_functions`` is then
    a mere list of coordinate expressions (and not a dictionary) and the
    arguments ``chart1`` and ``chart2`` have to be provided if the charts
    differ from the default ones on the domain and/or the codomain::

        sage: Phi1 = M.diff_map(N, [2*x/(1+x^2+y^2), 2*y/(1+x^2+y^2),
        ....:                       (x^2+y^2-1)/(1+x^2+y^2)],
        ....:                   chart1=c_xy, chart2=c_cart, name='Phi',
        ....:                   latex_name=r'\Phi')

    Since ``c_xy`` and ``c_cart`` are the default charts on respectively ``M``
    and ``N``, they can be omitted, so that the above declaration is equivalent
    to::

        sage: Phi1 = M.diff_map(N, [2*x/(1+x^2+y^2), 2*y/(1+x^2+y^2),
        ....:                       (x^2+y^2-1)/(1+x^2+y^2)],
        ....:                   name='Phi', latex_name=r'\Phi')

    With such a declaration, the differentiable map is only partially defined
    on the manifold `S^2`, being known in only one chart::

        sage: Phi1.display()
        Phi: S^2 → R^3
        on U: (x, y) ↦ (X, Y, Z) = (2*x/(x^2 + y^2 + 1), 2*y/(x^2 + y^2 + 1),
                                    (x^2 + y^2 - 1)/(x^2 + y^2 + 1))

    The definition can be completed by means of the method
    :meth:`~sage.manifolds.continuous_map.ContinuousMap.add_expr`::

        sage: Phi1.add_expr(c_uv, c_cart, [2*u/(1+u^2+v^2), 2*v/(1+u^2+v^2),
        ....:                              (1-u^2-v^2)/(1+u^2+v^2)])
        sage: Phi1.display()
        Phi: S^2 → R^3
        on U: (x, y) ↦ (X, Y, Z) = (2*x/(x^2 + y^2 + 1), 2*y/(x^2 + y^2 + 1),
                                    (x^2 + y^2 - 1)/(x^2 + y^2 + 1))
        on V: (u, v) ↦ (X, Y, Z) = (2*u/(u^2 + v^2 + 1), 2*v/(u^2 + v^2 + 1),
                                    -(u^2 + v^2 - 1)/(u^2 + v^2 + 1))

    At this stage, ``Phi1`` and ``Phi`` are fully equivalent::

        sage: Phi1 == Phi
        True

    The test suite is passed::

        sage: TestSuite(Phi).run()
        sage: TestSuite(Phi1).run()

    The map acts on points::

        sage: np = M.point((0,0), chart=c_uv, name='N')  # the North pole
        sage: Phi(np)
        Point Phi(N) on the 3-dimensional differentiable manifold R^3
        sage: Phi(np).coord() # Cartesian coordinates
        (0, 0, 1)
        sage: sp = M.point((0,0), chart=c_xy, name='S')  # the South pole
        sage: Phi(sp).coord() # Cartesian coordinates
        (0, 0, -1)

    The differential `\mathrm{d}\Phi` of the map `\Phi` at the North pole and
    at the South pole::

        sage: Phi.differential(np)
        Generic morphism:
          From: Tangent space at Point N on the 2-dimensional differentiable manifold S^2
          To:   Tangent space at Point Phi(N) on the 3-dimensional differentiable manifold R^3
        sage: Phi.differential(sp)
        Generic morphism:
          From: Tangent space at Point S on the 2-dimensional differentiable manifold S^2
          To:   Tangent space at Point Phi(S) on the 3-dimensional differentiable manifold R^3

    The matrix of the linear map `\mathrm{d}\Phi_N` with respect to the default
    bases of `T_N S^2` and `T_{\Phi(N)} \RR^3`::

        sage: Phi.differential(np).matrix()
        [2 0]
        [0 2]
        [0 0]

    the default bases being::

        sage: Phi.differential(np).domain().default_basis()
        Basis (∂/∂u,∂/∂v) on the Tangent space at Point N on the 2-dimensional
         differentiable manifold S^2
        sage: Phi.differential(np).codomain().default_basis()
        Basis (∂/∂X,∂/∂Y,∂/∂Z) on the Tangent space at Point Phi(N) on the
         3-dimensional differentiable manifold R^3

    A convenient way to display the matrix of the differential::

        sage: Phi.differential(np).display()
             ∂/∂u ∂/∂v
        ∂/∂X⎛   2    0⎞
        ∂/∂Y⎜   0    2⎟
        ∂/∂Z⎝   0    0⎠

    Differentiable maps can be composed by means of the operator ``*``: let
    us introduce the map `\RR^3\rightarrow \RR^2` corresponding to
    the projection from the point `(X,Y,Z)=(0,0,1)` onto the equatorial plane
    `Z=0`::

        sage: P = Manifold(2, 'R^2', r'\RR^2') # R^2 (equatorial plane)
        sage: cP.<xP, yP> = P.chart()
        sage: Psi = N.diff_map(P, (X/(1-Z), Y/(1-Z)), name='Psi',
        ....:                  latex_name=r'\Psi')
        sage: Psi
        Differentiable map Psi from the 3-dimensional differentiable manifold
         R^3 to the 2-dimensional differentiable manifold R^2
        sage: Psi.display()
        Psi: R^3 → R^2
           (X, Y, Z) ↦ (xP, yP) = (-X/(Z - 1), -Y/(Z - 1))

    Then we compose ``Psi`` with ``Phi``, thereby getting a map
    `S^2\rightarrow \RR^2`::

        sage: ster = Psi*Phi ; ster
        Differentiable map from the 2-dimensional differentiable manifold S^2
         to the 2-dimensional differentiable manifold R^2

    Let us test on the South pole (``sp``) that ``ster`` is indeed the
    composite of ``Psi`` and ``Phi``::

        sage: ster(sp) == Psi(Phi(sp))
        True

    Actually ``ster`` is the stereographic projection from the North pole, as
    its coordinate expression reveals::

        sage: ster.display()
        S^2 → R^2
        on U: (x, y) ↦ (xP, yP) = (x, y)
        on V: (u, v) ↦ (xP, yP) = (u/(u^2 + v^2), v/(u^2 + v^2))

    If its codomain is 1-dimensional, a differentiable map must be
    defined by a single symbolic expression for each pair of charts, and not
    by a list/tuple with a single element::

        sage: N = Manifold(1, 'N')
        sage: c_N = N.chart('X')
        sage: Phi = M.diff_map(N, {(c_xy, c_N): x^2+y^2,
        ....: (c_uv, c_N): 1/(u^2+v^2)})  # not ...[1/(u^2+v^2)] or (1/(u^2+v^2),)

    An example of differentiable map `\RR \rightarrow \RR^2`::

        sage: R = Manifold(1, 'R')  # field R
        sage: T.<t> = R.chart()  # canonical chart on R
        sage: R2 = Manifold(2, 'R^2')  # R^2
        sage: c_xy.<x,y> = R2.chart() # Cartesian coordinates on R^2
        sage: Phi = R.diff_map(R2, [cos(t), sin(t)], name='Phi') ; Phi
        Differentiable map Phi from the 1-dimensional differentiable manifold R
         to the 2-dimensional differentiable manifold R^2
        sage: Phi.parent()
        Set of Morphisms from 1-dimensional differentiable manifold R to
         2-dimensional differentiable manifold R^2 in Category of smooth
         manifolds over Real Field with 53 bits of precision
        sage: Phi.parent() is Hom(R, R2)
        True
        sage: Phi.display()
        Phi: R → R^2
           t ↦ (x, y) = (cos(t), sin(t))

    An example of diffeomorphism between the unit open disk and the Euclidean
    plane `\RR^2`::

        sage: D = R2.open_subset('D', coord_def={c_xy: x^2+y^2<1}) # the open unit disk
        sage: Phi = D.diffeomorphism(R2, [x/sqrt(1-x^2-y^2), y/sqrt(1-x^2-y^2)],
        ....:                        name='Phi', latex_name=r'\Phi')
        sage: Phi
        Diffeomorphism Phi from the Open subset D of the 2-dimensional
         differentiable manifold R^2 to the 2-dimensional differentiable
         manifold R^2
        sage: Phi.parent()
        Set of Morphisms from Open subset D of the 2-dimensional differentiable
         manifold R^2 to 2-dimensional differentiable manifold R^2 in Category
         of smooth manifolds over Real Field with 53 bits of precision
        sage: Phi.parent() is Hom(D, R2)
        True
        sage: Phi.display()
        Phi: D → R^2
           (x, y) ↦ (x, y) = (x/sqrt(-x^2 - y^2 + 1), y/sqrt(-x^2 - y^2 + 1))

    The image of a point::

        sage: p = D.point((1/2,0))
        sage: q = Phi(p) ; q
        Point on the 2-dimensional differentiable manifold R^2
        sage: q.coord()
        (1/3*sqrt(3), 0)

    The inverse diffeomorphism is computed by means of the method
    :meth:`~sage.manifolds.continuous_map.ContinuousMap.inverse`::

        sage: Phi.inverse()
        Diffeomorphism Phi^(-1) from the 2-dimensional differentiable manifold R^2
         to the Open subset D of the 2-dimensional differentiable manifold R^2
        sage: Phi.inverse().display()
        Phi^(-1): R^2 → D
           (x, y) ↦ (x, y) = (x/sqrt(x^2 + y^2 + 1), y/sqrt(x^2 + y^2 + 1))

    Equivalently, one may use the notations ``^(-1)`` or ``~`` to get the
    inverse::

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

    The coordinate expression of the inverse diffeomorphism::

        sage: (~Phi).display()
        Phi^(-1): R^2 → D
           (x, y) ↦ (x, y) = (x/sqrt(x^2 + y^2 + 1), y/sqrt(x^2 + y^2 + 1))

    A special case of diffeomorphism: the identity map of the open unit disk::

        sage: id = D.identity_map() ; id
        Identity map Id_D of the Open subset D of the 2-dimensional
         differentiable manifold R^2
        sage: latex(id)
        \mathrm{Id}_{D}
        sage: id.parent()
        Set of Morphisms from Open subset D of the 2-dimensional differentiable
         manifold R^2 to Open subset D of the 2-dimensional differentiable
         manifold R^2 in Join of Category of subobjects of sets and Category of
         smooth manifolds over Real Field with 53 bits of precision
        sage: id.parent() is Hom(D, D)
        True
        sage: id is Hom(D,D).one()  # the identity element of the monoid Hom(D,D)
        True

    The identity map acting on a point::

        sage: id(p)
        Point on the 2-dimensional differentiable manifold R^2
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
    def __init__(self, parent, coord_functions=None, name=None,
                 latex_name=None, is_isomorphism=False, is_identity=False):
        r"""
        Construct a differentiable map.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: N = Manifold(3, 'N')
            sage: Y.<u,v,w> = N.chart()
            sage: f = Hom(M,N)({(X,Y): (x+y, x*y, x-y)}, name='f') ; f
            Differentiable map f from the 2-dimensional differentiable manifold
             M to the 3-dimensional differentiable manifold N
            sage: f.display()
            f: M → N
               (x, y) ↦ (u, v, w) = (x + y, x*y, x - y)
            sage: TestSuite(f).run()

        The identity map::

            sage: f = Hom(M,M)({}, is_identity=True) ; f
            Identity map Id_M of the 2-dimensional differentiable manifold M
            sage: f.display()
            Id_M: M → M
               (x, y) ↦ (x, y)
            sage: TestSuite(f).run()
        """
        ContinuousMap.__init__(self, parent, coord_functions=coord_functions,
                               name=name, latex_name=latex_name,
                               is_isomorphism=is_isomorphism,
                               is_identity=is_identity)

    #
    # SageObject methods
    #

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: N = Manifold(2, 'N')
            sage: Y.<u,v> = N.chart()
            sage: f = Hom(M,N)({(X,Y): (x+y,x*y)})
            sage: f._repr_()
            'Differentiable map from the 2-dimensional differentiable manifold M to
             the 2-dimensional differentiable manifold N'
            sage: f = Hom(M,N)({(X,Y): (x+y,x*y)}, name='f')
            sage: f._repr_()
            'Differentiable map f from the 2-dimensional differentiable manifold M to
             the 2-dimensional differentiable manifold N'
            sage: f = Hom(M,N)({(X,Y): (x+y,x-y)}, name='f', is_isomorphism=True)
            sage: f._repr_()
            'Diffeomorphism f from the 2-dimensional differentiable manifold M to
             the 2-dimensional differentiable manifold N'
            sage: f = Hom(M,M)({(X,X): (x+y,x-y)}, name='f', is_isomorphism=True)
            sage: f._repr_()
            'Diffeomorphism f of the 2-dimensional differentiable manifold M'
            sage: f = Hom(M,M)({}, name='f', is_identity=True)
            sage: f._repr_()
            'Identity map f of the 2-dimensional differentiable manifold M'
        """
        if self._is_identity:
            return "Identity map " + self._name + \
                   " of the {}".format(self._domain)
        if self._is_isomorphism:
            description = "Diffeomorphism"
        else:
            description = "Differentiable map"
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

    def _init_derived(self):
        r"""
        Initialize the derived quantities.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = M.diffeomorphism(M, [x+y, x-y])
            sage: f._init_derived()
            sage: f._restrictions
            {}
            sage: f._inverse
        """
        ContinuousMap._init_derived(self)
        # derived quantities of the mother class
        self._diff = {}
        # dict. of the coord. expressions of the differential
        # keys: pair of charts

    def _del_derived(self):
        r"""
        Delete the derived quantities.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: f = M.diffeomorphism(M, [x+y, x-y])
            sage: f^(-1)
            Diffeomorphism of the 2-dimensional differentiable manifold M
            sage: f._inverse  # was set by f^(-1)
            Diffeomorphism of the 2-dimensional differentiable manifold M
            sage: f._del_derived()
            sage: f._inverse  # has been set to None by _del_derived()
        """
        ContinuousMap._del_derived(self)  # derived quantities of the mother
                                          # class
        self._diff.clear()

    def differential(self, point: ManifoldPoint) -> FiniteRankFreeModuleMorphism:
        r"""
        Return the differential of ``self`` at a given point.

        If the differentiable map ``self`` is

        .. MATH::

            \Phi: M \longrightarrow N,

        where `M` and `N` are differentiable manifolds, the *differential*
        of `\Phi` at a point `p \in M` is the tangent space linear map:

        .. MATH::

            \mathrm{d}\Phi_p: T_p M \longrightarrow T_{\Phi(p)} N

        defined by

        .. MATH::

            \begin{array}{rccc}
            \forall v\in T_p M,\quad \mathrm{d}\Phi_p(v) :
                                & C^k(N) & \longrightarrow & \mathbb{R} \\
                                & f & \longmapsto & v(f\circ \Phi)
            \end{array}

        INPUT:

        - ``point`` -- point `p` in the domain `M` of the differentiable
          map `\Phi`

        OUTPUT:

        - `\mathrm{d}\Phi_p`, the differential of `\Phi` at `p`, as a
          :class:`~sage.tensor.modules.free_module_morphism.FiniteRankFreeModuleMorphism`

        EXAMPLES:

        Differential of a differentiable map between a 2-dimensional manifold
        and a 3-dimensional one::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: N = Manifold(3, 'N')
            sage: Y.<u,v,w> = N.chart()
            sage: Phi = M.diff_map(N, {(X,Y): (x-2*y, x*y, x^2-y^3)}, name='Phi',
            ....:                  latex_name = r'\Phi')
            sage: p = M.point((2,-1), name='p')
            sage: dPhip = Phi.differential(p) ; dPhip
            Generic morphism:
              From: Tangent space at Point p on the 2-dimensional differentiable manifold M
              To:   Tangent space at Point Phi(p) on the 3-dimensional differentiable manifold N
            sage: latex(dPhip)
            {\mathrm{d}\Phi}_{p}
            sage: dPhip.parent()
            Set of Morphisms from Tangent space at Point p on the 2-dimensional
             differentiable manifold M to Tangent space at Point Phi(p) on the
             3-dimensional differentiable manifold N in Category of finite
             dimensional vector spaces over Symbolic Ring

        The matrix of `\mathrm{d}\Phi_p` w.r.t. to the default bases of
        `T_p M` and `T_{\Phi(p)} N`::

            sage: dPhip.matrix()
            [ 1 -2]
            [-1  2]
            [ 4 -3]
        """
        image_point = self(point)
        tsp_image = image_point._manifold.tangent_space(image_point)
        tsp_source = point._manifold.tangent_space(point)
        # Search for a common chart to perform the computation
        chartp = None

        # 1/ Search without any extra computation
        for chart in point._coordinates:
            for chart_pair in self._diff:
                if chart == chart_pair[0]:
                    chartp = chart_pair
                    break
            if chartp is not None:
                break
        else:
            # 2/ Search with a coordinate transformation on the point
            for chart_pair in self._diff:
                try:
                    point.coord(chart_pair[0])
                    chartp = chart_pair
                except ValueError:
                    pass

        if chartp is None:
            # 3/ Search with a coordinate evaluation of self
            for chart1 in point._coordinates:
                for chart2 in self._codomain.atlas():
                    try:
                        self.differential_functions(chart1, chart2)
                        chartp = (chart1, chart2)
                        break
                    except ValueError:
                        pass
                if chartp is not None:
                    break

        if chartp is None:
            raise ValueError("no common chart have been found for the " +
                     "coordinate expressions of {} and {}".format(self, point))

        diff_funct = self.differential_functions(*chartp)
        chart1 = chartp[0]
        chart2 = chartp[1]
        coord_point = point.coord(chart1)
        n1 = self._domain.dim()
        n2 = self._codomain.dim()
        matrix = [[diff_funct[i][j](*coord_point) for j in range(n1)]
                  for i in range(n2)]
        bases = (chart1.frame().at(point), chart2.frame().at(image_point))
        if self._name is not None and point._name is not None:
            name = 'd%s_%s' % (self._name, point._name)
        else:
            name = None
        if self._latex_name is not None and point._latex_name is not None:
            latex_name = r'{\mathrm{d}%s}_{%s}' % (self._latex_name,
                                                 point._latex_name)
        else:
            latex_name = None
        return tsp_source.hom(tsp_image, matrix, bases=bases,
                              name=name, latex_name=latex_name)

    def differential_functions(self, chart1=None, chart2=None):
        r"""
        Return the coordinate expression of the differential of the
        differentiable map with respect to a pair of charts.

        If the differentiable map is

        .. MATH::

            \Phi: M \longrightarrow  N,

        where `M` and `N` are differentiable manifolds, the *differential*
        of `\Phi` at a point `p \in M` is the tangent space linear map:

        .. MATH::

            \mathrm{d}\Phi_p: T_p M \longrightarrow T_{\Phi(p)} N

        defined by

        .. MATH::

            \begin{array}{rccc}
            \forall v\in T_p M,\quad \mathrm{d}\Phi_p(v) : & C^k(N) &
                                              \longrightarrow & \mathbb{R}, \\
                                & f & \longmapsto & v(f\circ \Phi).
            \end{array}

        If the coordinate expression of `\Phi` is

        .. MATH::

            y^i = Y^i(x^1, \ldots, x^n), \quad 1 \leq i \leq m,

        where `(x^1, \ldots, x^n)` are coordinates of a chart on `M` and
        `(y^1, \ldots, y^m)` are coordinates of a chart on `\Phi(M)`, the
        expression of the differential of `\Phi` with respect to these
        coordinates is

        .. MATH::

            J_{ij} = \frac{\partial Y^i}{\partial x^j} \quad 1\leq i \leq m,
                            \qquad 1 \leq j \leq n.

        `\left. J_{ij} \right|_p` is then the matrix of the linear map
        `\mathrm{d}\Phi_p` with respect to the bases of `T_p M` and
        `T_{\Phi(p)} N` associated to the above charts:

        .. MATH::

            \mathrm{d}\Phi_p\left( \left. \frac{\partial}{\partial x^j}
               \right|_p \right) = \left. J_{ij} \right|_p \;
             \left. \frac{\partial}{\partial y^i} \right| _{\Phi(p)}.

        INPUT:

        - ``chart1`` -- (default: ``None``) chart on the domain `M` of `\Phi`
          (coordinates denoted by `(x^j)` above); if ``None``, the domain's
          default chart is assumed
        - ``chart2`` -- (default: ``None``) chart on the codomain of `\Phi`
          (coordinates denoted by `(y^i)` above); if ``None``, the codomain's
          default chart is assumed

        OUTPUT:

        - the functions `J_{ij}` as a double array, `J_{ij}` being
          the element ``[i][j]`` represented by a
          :class:`~sage.manifolds.chart_func.ChartFunction`

        To get symbolic expressions, use the method
        :meth:`jacobian_matrix` instead.

        EXAMPLES:

        Differential functions of a map between a 2-dimensional manifold
        and a 3-dimensional one::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: N = Manifold(3, 'N')
            sage: Y.<u,v,w> = N.chart()
            sage: Phi = M.diff_map(N, {(X,Y): (x-2*y, x*y, x^2-y^3)}, name='Phi',
            ....:                  latex_name = r'\Phi')
            sage: J = Phi.differential_functions(X, Y) ; J
            [     1     -2]
            [     y      x]
            [   2*x -3*y^2]

        The result is cached::

            sage: Phi.differential_functions(X, Y) is J
            True

        The elements of ``J`` are functions of the coordinates of
        the chart ``X``::

            sage: J[2][0]
            2*x
            sage: type(J[2][0])
            <class 'sage.manifolds.chart_func.ChartFunctionRing_with_category.element_class'>

            sage: J[2][0].display()
            (x, y) ↦ 2*x

        In contrast, the method :meth:`jacobian_matrix` leads directly to
        symbolic expressions::

            sage: JJ = Phi.jacobian_matrix(X,Y) ; JJ
            [     1     -2]
            [     y      x]
            [   2*x -3*y^2]
            sage: JJ[2,0]
            2*x
            sage: type(JJ[2,0])
            <class 'sage.symbolic.expression.Expression'>
            sage: bool( JJ[2,0] == J[2][0].expr() )
            True
        """
        dom1 = self._domain
        dom2 = self._codomain
        if chart1 is None:
            chart1 = dom1._def_chart
        if chart2 is None:
            chart2 = dom2._def_chart
        if (chart1, chart2) not in self._diff:
            funct = self.coord_functions(chart1, chart2)
            self._diff[(chart1, chart2)] = funct.jacobian()
        return self._diff[(chart1, chart2)]

    def jacobian_matrix(self, chart1=None, chart2=None):
        r"""
        Return the Jacobian matrix resulting from the coordinate expression
        of the differentiable map with respect to a pair of charts.

        If `\Phi` is the current differentiable map and its coordinate
        expression is

        .. MATH::

            y^i = Y^i(x^1, \ldots, x^n),  \quad 1 \leq i \leq m,

        where `(x^1, \ldots, x^n)` are coordinates of a chart `X` on the
        domain of `\Phi` and `(y^1, \ldots, y^m)` are coordinates of a chart
        `Y` on the codomain of `\Phi`, the *Jacobian matrix* of the
        differentiable map `\Phi` w.r.t. to charts `X` and `Y` is

        .. MATH::

            J = \left( \frac{\partial Y^i}{\partial x^j}
              \right)_{\substack{1 \leq i \leq m \\ 1 \leq j \leq n}},

        where `i` is the row index and `j` the column one.

        INPUT:

        - ``chart1`` -- (default: ``None``) chart `X` on the domain of
          `\Phi`; if none is provided, the domain's default chart is assumed
        - ``chart2`` -- (default: ``None``) chart `Y` on the codomain of
          `\Phi`; if none is provided, the codomain's default chart is
          assumed

        OUTPUT: the matrix `J` defined above

        EXAMPLES:

        Jacobian matrix of a map between a 2-dimensional manifold and a
        3-dimensional one::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: N = Manifold(3, 'N')
            sage: Y.<u,v,w> = N.chart()
            sage: Phi = M.diff_map(N, {(X,Y): (x-2*y, x*y, x^2-y^3)}, name='Phi',
            ....:                  latex_name = r'\Phi')
            sage: Phi.display()
            Phi: M → N
               (x, y) ↦ (u, v, w) = (x - 2*y, x*y, -y^3 + x^2)
            sage: J = Phi.jacobian_matrix(X, Y) ; J
            [     1     -2]
            [     y      x]
            [   2*x -3*y^2]
            sage: J.parent()
            Full MatrixSpace of 3 by 2 dense matrices over Symbolic Ring
        """
        from sage.matrix.constructor import matrix
        diff_funct = self.differential_functions(chart1, chart2)
        n1 = self._domain.dim()
        n2 = self._codomain.dim()
        return matrix([[diff_funct[i][j].expr() for j in range(n1)]
                       for i in range(n2)])

    def pullback(self, tensor_or_codomain_subset, name=None, latex_name=None):
        r"""
        Pullback operator associated with ``self``.

        In what follows, let `\Phi` denote a differentiable map ``self``,
        `M` its domain and `N` its codomain.

        INPUT:

        One of the following:

        - ``tensor_or_codomain_subset`` -- one of the following:

          - a :class:`~sage.manifolds.differentiable.tensorfield.TensorField`;
            a fully covariant tensor field `T` on `N`, i.e. a tensor
            field of type `(0, p)`, with `p` a positive or zero integer; the
            case `p = 0` corresponds to a scalar field
          - a :class:`~sage.manifolds.subset.ManifoldSubset` `S`

        OUTPUT:

        - (if the input is a tensor field `T`)
          a :class:`~sage.manifolds.differentiable.tensorfield.TensorField`
          representing a fully covariant tensor field on `M` that is the
          pullback of `T` by `\Phi`
        - (if the input is a manifold subset `S`)
          a :class:`~sage.manifolds.subset.ManifoldSubset` that is the
          preimage `\Phi^{-1}(S)`; same as :meth:`preimage`

        EXAMPLES:

        Pullback on `S^2` of a scalar field defined on `R^3`::

            sage: M = Manifold(2, 'S^2', start_index=1)
            sage: U = M.open_subset('U') # the complement of a meridian (domain of spherical coordinates)
            sage: c_spher.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi') # spherical coord. on U
            sage: N = Manifold(3, 'R^3', r'\RR^3', start_index=1)
            sage: c_cart.<x,y,z> = N.chart() # Cartesian coord. on R^3
            sage: Phi = U.diff_map(N, (sin(th)*cos(ph), sin(th)*sin(ph), cos(th)),
            ....:                  name='Phi', latex_name=r'\Phi')
            sage: f = N.scalar_field(x*y*z, name='f') ; f
            Scalar field f on the 3-dimensional differentiable manifold R^3
            sage: f.display()
            f: R^3 → ℝ
               (x, y, z) ↦ x*y*z
            sage: pf = Phi.pullback(f) ; pf
            Scalar field Phi^*(f) on the Open subset U of the 2-dimensional
             differentiable manifold S^2
            sage: pf.display()
            Phi^*(f): U → ℝ
               (th, ph) ↦ cos(ph)*cos(th)*sin(ph)*sin(th)^2

        Pullback on `S^2` of the standard Euclidean metric on `R^3`::

            sage: g = N.sym_bilin_form_field(name='g')
            sage: g[1,1], g[2,2], g[3,3] = 1, 1, 1
            sage: g.display()
            g = dx⊗dx + dy⊗dy + dz⊗dz
            sage: pg = Phi.pullback(g) ; pg
            Field of symmetric bilinear forms Phi^*(g) on the Open subset U of
             the 2-dimensional differentiable manifold S^2
            sage: pg.display()
            Phi^*(g) = dth⊗dth + sin(th)^2 dph⊗dph

        Parallel computation::

           sage: Parallelism().set('tensor', nproc=2)
           sage: pg = Phi.pullback(g) ; pg
           Field of symmetric bilinear forms Phi^*(g) on the Open subset U of
            the 2-dimensional differentiable manifold S^2
           sage: pg.display()
           Phi^*(g) = dth⊗dth + sin(th)^2 dph⊗dph
           sage: Parallelism().set('tensor', nproc=1)  # switch off parallelization


        Pullback on `S^2` of a 3-form on `R^3`::

            sage: a = N.diff_form(3, name='A')
            sage: a[1,2,3] = f
            sage: a.display()
            A = x*y*z dx∧dy∧dz
            sage: pa = Phi.pullback(a) ; pa
            3-form Phi^*(A) on the Open subset U of the 2-dimensional
             differentiable manifold S^2
            sage: pa.display() # should be zero (as any 3-form on a 2-dimensional manifold)
            Phi^*(A) = 0

        TESTS:

        Check that :issue:`31904` is fixed::

            sage: E.<x,y> = EuclideanSpace()
            sage: polar.<r,ph> = E.polar_coordinates()
            sage: g = E.metric()
            sage: M = Manifold(1, 'M')
            sage: Ct.<t> = M.chart()
            sage: F = M.diff_map(E, coord_functions={(Ct, polar): (1 + cos(t), t)})
            sage: gM = F.pullback(g)
            sage: gM.display()
            (2*cos(t) + 2) dt⊗dt
        """
        if not hasattr(tensor_or_codomain_subset, '_domain'):
            return super().pullback(tensor_or_codomain_subset,
                                    name=name, latex_name=latex_name)
        tensor = tensor_or_codomain_subset

        from sage.manifolds.differentiable.tensorfield_paral import TensorFieldParal
        from sage.tensor.modules.comp import (
            CompFullyAntiSym,
            CompFullySym,
            Components,
            CompWithSym,
        )

        def _pullback_chart(diff_map, tensor, chart1, chart2):
            r"""
            Helper function performing the pullback on chart domains
            only.

            INPUT:

            - ``diff_map`` -- a restriction of ``self``, whose both
              domain and codomain are chart domains, corresponding
              respectively to ``chart1`` and ``chart2``
            - ``tensor`` -- a covariant tensor field, whose domain is
              the codomain of ``diff_map`` and whose components are known
              in ``chart2.frame()``
            - ``chart1`` -- chart on the domain of ``diff_map``
            - ``chart2`` -- chart on the codomain of ``diff_map``

            OUTPUT: the pull back of ``tensor`` by ``diff_map``
            """
            dom1 = diff_map._domain
            dom2 = diff_map._codomain
            ncov = tensor._tensor_type[1]
            resu_name = None
            resu_latex_name = None
            if diff_map._name is not None and tensor._name is not None:
                resu_name = diff_map._name + '^*(' + tensor._name + ')'
            if (diff_map._latex_name is not None and
                tensor._latex_name is not None):
                resu_latex_name = '{' + diff_map._latex_name + '}^*' \
                                  + tensor._latex_name
            fmodule1 = dom1.vector_field_module()
            ring1 = fmodule1._ring
            si1 = fmodule1._sindex
            of1 = fmodule1._output_formatter
            si2 = dom2._sindex
            resu = fmodule1.tensor((0, ncov), name=resu_name,
                                   latex_name=resu_latex_name, sym=tensor._sym,
                                   antisym=tensor._antisym)

            nproc = Parallelism().get('tensor')
            ind_old_list = list(dom2.manifold().index_generator(ncov))

            frame1 = chart1.frame()
            frame2 = chart2.frame()

            tcomp = tensor._components[frame2]
            if isinstance(tcomp, CompFullySym):
                ptcomp = CompFullySym(ring1, frame1, ncov, start_index=si1,
                                      output_formatter=of1)
            elif isinstance(tcomp, CompFullyAntiSym):
                ptcomp = CompFullyAntiSym(ring1, frame1, ncov, start_index=si1,
                                          output_formatter=of1)
            elif isinstance(tcomp, CompWithSym):
                ptcomp = CompWithSym(ring1, frame1, ncov, start_index=si1,
                                     output_formatter=of1, sym=tcomp.sym,
                                     antisym=tcomp.antisym)
            else:
                ptcomp = Components(ring1, frame1, ncov, start_index=si1,
                                    output_formatter=of1)
            phi = diff_map._coord_expression[(chart1, chart2)]
            jacob = phi.jacobian()
            # X2 coordinates expressed in terms of X1 ones via the diff. map:
            coord2_1 = phi(*(chart1._xx))

            if nproc != 1:
                # Parallel computation
                lol = lambda lst, sz: [lst[i:i+sz] for i in range(0, len(lst), sz)]
                ind_list = list(ptcomp.non_redundant_index_generator())
                ind_step = max(1, (len(ind_list) // nproc) // 2)
                local_list = lol(ind_list, ind_step)
                # list of input parameters
                listParalInput = [(tcomp, chart1, chart2, coord2_1, jacob,
                                   ind_old_list, si1, si2, ncov, ind_part)
                                  for ind_part in local_list]

                @parallel(p_iter='multiprocessing', ncpus=nproc)
                def paral_comp(tcomp, chart1, chart2, coord2_1, jacob,
                               ind_old_list, si1, si2, ncov, local_list_ind):
                    partial = []
                    for ind_new in local_list_ind:
                        res = 0
                        for ind_old in ind_old_list:
                            ff = tcomp[[ind_old]].coord_function(chart2)
                            t = chart1.function(ff(*coord2_1))
                            for i in range(ncov):
                                t *= jacob[ind_old[i]-si2, ind_new[i]-si1]
                            res += t
                        partial.append([ind_new, res])
                    return partial

                for ii, val in paral_comp(listParalInput):
                    for jj in val:
                        ptcomp[[jj[0]]] = jj[1]

            else:
                # Sequential computation
                for ind_new in ptcomp.non_redundant_index_generator():
                    res = 0
                    for ind_old in ind_old_list:
                        ff = tcomp[[ind_old]].coord_function(chart2)
                        t = chart1.function(ff(*coord2_1))
                        for i in range(ncov):
                            t *= jacob[ind_old[i]-si2, ind_new[i]-si1]
                        res += t
                    ptcomp[ind_new] = res

            resu._components[frame1] = ptcomp
            return resu
        # End of function _pullback_chart

        # Special case of the identity map:
        if self._is_identity:
            return tensor  # no test for efficiency
        # Generic case:
        dom1 = self._domain
        dom2 = self._codomain
        tdom = tensor._domain
        if not tdom.is_subset(dom2):
            raise ValueError("the tensor field is not defined on the map codomain")
        (ncon, ncov) = tensor._tensor_type
        if ncon != 0:
            raise TypeError("the pullback cannot be taken on a tensor " +
                            "with some contravariant part")
        resu_name = None
        resu_latex_name = None
        if self._name is not None and tensor._name is not None:
            resu_name = self._name + '^*(' + tensor._name + ')'
        if self._latex_name is not None and tensor._latex_name is not None:
            resu_latex_name = '{' + self._latex_name + '}^*' \
                              + tensor._latex_name
        if ncov == 0:
            # Case of a scalar field
            resu_fc = []
            for chart2 in tensor._express:
                for chart1 in dom1._atlas:
                    if (chart1, chart2) in self._coord_expression:
                        phi = self._coord_expression[(chart1, chart2)]
                        coord1 = chart1._xx
                        ff = tensor._express[chart2]
                        resu_fc.append( chart1.function(ff(*(phi(*coord1)))) )
            dom_resu = resu_fc[0].parent()._chart.domain()
            for fc in resu_fc[1:]:
                dom_resu = dom_resu.union(fc.parent()._chart.domain())
            resu = dom_resu.scalar_field(name=resu_name,
                                         latex_name=resu_latex_name)
            for fc in resu_fc:
                resu._express[fc.parent()._chart] = fc
        else:
            # Case of tensor field of rank >= 1
            if tensor._vmodule._dest_map is not tdom.identity_map():
                raise TypeError("the pullback is defined only for tensor " +
                                "fields on {}".format(dom2))
            resu_rst = []
            for chart_pair in self._coord_expression:
                chart1 = chart_pair[0]
                chart2 = chart_pair[1]
                ch2dom = chart2._domain
                if ch2dom.is_subset(tdom):
                    self_r = self.restrict(chart1._domain, subcodomain=ch2dom)
                    tensor_r = tensor.restrict(ch2dom)
                    if chart2.frame() in tensor_r._components:
                        resu_rst.append(_pullback_chart(self_r, tensor_r,
                                                        chart1, chart2))
            dom_resu = resu_rst[0]._domain
            for rst in resu_rst[1:]:
                dom_resu = dom_resu.union(rst._domain)
            resu = dom_resu.tensor_field(0, ncov, name=resu_name,
                                         latex_name=resu_latex_name,
                                         sym=resu_rst[0]._sym,
                                         antisym=resu_rst[0]._antisym)
            for rst in resu_rst:
                if rst._domain is not resu._domain:
                    resu._restrictions[rst._domain] = rst
            if isinstance(resu, TensorFieldParal):
                for rst in resu_rst:
                    if rst._domain is resu._domain:
                        for frame, comp in rst._components.items():
                            resu._components[frame] = comp
        return resu

    def pushforward(self, tensor):
        r"""
        Pushforward operator associated with ``self``.

        In what follows, let `\Phi` denote the differentiable map, `M` its
        domain and `N` its codomain.

        INPUT:

        - ``tensor`` --
          :class:`~sage.manifolds.differentiable.tensorfield.TensorField`;
          a fully contrariant tensor field `T` on `M`, i.e. a tensor
          field of type `(p, 0)`, with `p` a positive integer

        OUTPUT:

        - a :class:`~sage.manifolds.differentiable.tensorfield.TensorField`
          representing a fully contravariant tensor field along `M` with
          values in `N`, which is the pushforward of `T` by `\Phi`

        EXAMPLES:

        Pushforward of a vector field on the 2-sphere `S^2` to the Euclidean
        3-space `\RR^3`, via the standard embedding of `S^2`::

            sage: S2 = Manifold(2, 'S^2', start_index=1)
            sage: U = S2.open_subset('U')  # domain of spherical coordinates
            sage: spher.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
            sage: R3 = Manifold(3, 'R^3', start_index=1)
            sage: cart.<x,y,z> = R3.chart()
            sage: Phi = U.diff_map(R3, {(spher, cart): [sin(th)*cos(ph),
            ....:   sin(th)*sin(ph), cos(th)]}, name='Phi', latex_name=r'\Phi')
            sage: v = U.vector_field(name='v')
            sage: v[:] = 0, 1
            sage: v.display()
            v = ∂/∂ph
            sage: pv = Phi.pushforward(v); pv
            Vector field Phi_*(v) along the Open subset U of the 2-dimensional
             differentiable manifold S^2 with values on the 3-dimensional
             differentiable manifold R^3
            sage: pv.display()
            Phi_*(v) = -sin(ph)*sin(th) ∂/∂x + cos(ph)*sin(th) ∂/∂y

        Pushforward of a vector field on the real line to the `\RR^3`, via a
        helix embedding::

            sage: R.<t> = manifolds.RealLine()
            sage: Psi = R.diff_map(R3, [cos(t), sin(t), t], name='Psi',
            ....:                  latex_name=r'\Psi')
            sage: u = R.vector_field(name='u')
            sage: u[0] = 1
            sage: u.display()
            u = ∂/∂t
            sage: pu = Psi.pushforward(u); pu
            Vector field Psi_*(u) along the Real number line ℝ with values on
             the 3-dimensional differentiable manifold R^3
            sage: pu.display()
            Psi_*(u) = -sin(t) ∂/∂x + cos(t) ∂/∂y + ∂/∂z
        """
        from sage.manifolds.differentiable.tensorfield_paral import TensorFieldParal
        from sage.tensor.modules.comp import (
            CompFullyAntiSym,
            CompFullySym,
            Components,
            CompWithSym,
        )
        vmodule = tensor.base_module()
        dest_map = vmodule.destination_map()
        dom1 = tensor.domain()
        ambient_dom1 = dest_map.codomain()
        if not ambient_dom1.is_subset(self._domain):
            raise ValueError("the {} does not take its ".format(tensor) +
                             "values on the domain of the {}".format(self))
        (ncon, ncov) = tensor.tensor_type()
        if ncov != 0:
            raise ValueError("the pushforward cannot be taken on a tensor " +
                             "with some covariant part")
        if ncon == 0:
            raise ValueError("the pushforward cannot be taken on a scalar " +
                             "field")
        if dest_map != dom1.identity_map():
            raise NotImplementedError("the case of a non-trivial destination" +
                                      " map is not implemented yet")
        if not isinstance(tensor, TensorFieldParal):
            raise NotImplementedError("the case of a non-parallelizable " +
                                      "domain is not implemented yet")
        # A pair of charts (chart1, chart2) where the computation
        # is feasible is searched, privileging the default chart of the
        # map's domain for chart1
        chart1 = None
        chart2 = None
        def_chart1 = dom1.default_chart()
        def_chart2 = self._codomain.default_chart()
        if (def_chart1._frame in tensor._components
                and (def_chart1, def_chart2) in self._coord_expression):
            chart1 = def_chart1
            chart2 = def_chart2
        else:
            for (chart1n, chart2n) in self._coord_expression:
                if (chart2n == def_chart2
                         and chart1n._frame in tensor._components):
                    chart1 = chart1n
                    chart2 = def_chart2
                    break
        if chart2 is None:
            # It is not possible to have def_chart2 as chart for
            # expressing the result; any other chart is then looked for:
            for (chart1n, chart2n) in self._coord_expression:
                if chart1n._frame in tensor._components:
                    chart1 = chart1n
                    chart2 = chart2n
                    break
        if chart1 is None:
            # Computations of components are attempted
            chart1 = dom1.default_chart()
            tensor.comp(chart1.frame())
            for chart2n in self._codomain.atlas():
                try:
                    self.coord_functions(chart1, chart2n)
                    chart2 = chart2n
                    break
                except ValueError:
                    pass
            else:
                raise ValueError("no pair of charts could be found to " +
                                 "compute the pushforward of " +
                                 "the {} by the {}".format(tensor, self))
        # Vector field module for the result:
        fmodule2 = dom1.vector_field_module(dest_map=self)

        frame2 = fmodule2.basis(from_frame=chart2.frame())
        si1 = dom1.start_index()
        si2 = fmodule2._sindex
        ring2 = fmodule2._ring
        of2 = fmodule2._output_formatter
        # Computation at the component level:
        tcomp = tensor._components[chart1.frame()]
        # Construction of the pushforward components (ptcomp):
        if isinstance(tcomp, CompFullySym):
            ptcomp = CompFullySym(ring2, frame2, ncon, start_index=si2,
                                  output_formatter=of2)
        elif isinstance(tcomp, CompFullyAntiSym):
            ptcomp = CompFullyAntiSym(ring2, frame2, ncon, start_index=si2,
                                      output_formatter=of2)
        elif isinstance(tcomp, CompWithSym):
            ptcomp = CompWithSym(ring2, frame2, ncon, start_index=si2,
                                 output_formatter=of2, sym=tcomp._sym,
                                 antisym=tcomp._antisym)
        else:
            ptcomp = Components(ring2, frame2, ncon, start_index=si2,
                                output_formatter=of2)
        # Computation of the pushforward components:
        jacob = self.differential_functions(chart1=chart1, chart2=chart2)
        si2 = chart2.domain().start_index()
        for ind_new in ptcomp.non_redundant_index_generator():
            res = 0
            for ind_old in dom1.index_generator(ncon):
                t = tcomp[[ind_old]].coord_function(chart1)
                for i in range(ncon):
                    t *= jacob[ind_new[i]-si2, ind_old[i]-si1]
                res += t
            ptcomp[ind_new] = res
        # Name of the result:
        resu_name = None
        resu_latex_name = None
        if self._name is not None and tensor._name is not None:
            resu_name = self._name + '_*(' + tensor._name + ')'
        if self._latex_name is not None and tensor._latex_name is not None:
            resu_latex_name = '{' + self._latex_name + '}_*' \
                              + tensor._latex_name
        # Creation of the result with the components obtained above:
        resu = fmodule2.tensor_from_comp((ncon, 0), ptcomp, name=resu_name,
                                         latex_name=resu_latex_name)
        return resu
