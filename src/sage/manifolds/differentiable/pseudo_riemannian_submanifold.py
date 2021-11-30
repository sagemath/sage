r"""
Pseudo-Riemannian submanifolds

An *embedded (resp. immersed) submanifold of a pseudo-Riemannian manifold*
`(M,g)` is an embedded (resp. immersed) submanifold `N` of `M` as a
differentiable manifold (see
:mod:`~sage.manifolds.differentiable.differentiable_submanifold`) such that
pull back of the metric tensor `g` via the embedding (resp. immersion) endows
`N` with the structure of a pseudo-Riemannian manifold.

The following example shows how to compute the various quantities related
to the intrinsic and extrinsic geometries of a hyperbolic slicing of the
3-dimensional Minkowski space.

We start by declaring the ambient manifold `M` and the submanifold `N`::

    sage: M = Manifold(3, 'M', structure="Lorentzian")
    sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian", start_index=1)

The considered slices being spacelike hypersurfaces, they are Riemannian
manifolds.

Let us introduce the Minkowskian coordinates `(w,x,y)` on `M` and the polar
coordinates `(\rho, \theta)` on the submanifold `N`::

    sage: E.<w,x,y> = M.chart()
    sage: C.<rh,th> = N.chart(r'rh:(0,+oo):\rho th:(0,2*pi):\theta')

Let `b` be the hyperbola semi-major axis and `t` the parameter of the
foliation::

    sage: b = var('b', domain='real')
    sage: assume(b>0)
    sage: t = var('t', domain='real')

One can then define the embedding `\phi_t`::

    sage: phi = N.diff_map(M, {(C,E): [b*cosh(rh)+t,
    ....:                              b*sinh(rh)*cos(th),
    ....:                              b*sinh(rh)*sin(th)]})
    sage: phi.display()
    N → M
       (rh, th) ↦ (w, x, y) = (b*cosh(rh) + t, b*cos(th)*sinh(rh),
                                  b*sin(th)*sinh(rh))

as well as its inverse (when considered as a diffeomorphism onto its image)::

    sage: phi_inv = M.diff_map(N, {(E,C): [log(sqrt(x^2+y^2+b^2)/b+
    ....:                                  sqrt((x^2+y^2+b^2)/b^2-1)),
    ....:                                  atan2(y,x)]})
    sage: phi_inv.display()
    M → N
       (w, x, y) ↦ (rh, th) = (log(sqrt((b^2 + x^2 + y^2)/b^2 - 1)
                                  + sqrt(b^2 + x^2 + y^2)/b), arctan2(y, x))

and the partial inverse expressing the foliation parameter `t` as a scalar
field on `M`::

    sage: phi_inv_t = M.scalar_field({E: w-sqrt(x^2+y^2+b^2)})
    sage: phi_inv_t.display()
    M → ℝ
    (w, x, y) ↦ w - sqrt(b^2 + x^2 + y^2)

One can check that the inverse is correct with::

    sage: (phi*phi_inv).display()
    M → M
       (w, x, y) ↦ ((b^2 + x^2 + y^2 + sqrt(b^2 + x^2 + y^2)*(t + sqrt(x^2 +
     y^2)) + sqrt(x^2 + y^2)*t)/(sqrt(b^2 + x^2 + y^2) + sqrt(x^2 + y^2)), x, y)

The first item of the 3-uple in the right-hand does not appear as `w` because
`t` has not been replaced by its value provided by ``phi_inv_t``. Once this is
done, we do get `w`::

    sage: (phi*phi_inv).expr()[0].subs({t: phi_inv_t.expr()}).simplify_full()
    w

The embedding can then be declared::

    sage: N.set_embedding(phi, inverse=phi_inv, var=t,
    ....:                 t_inverse = {t: phi_inv_t})

This line does not perform any calculation yet. It just check the coherence of
the arguments, but not the inverse, the user is trusted on this point.

Finally, we initialize the metric of `M` to be that of Minkowski space::

    sage: g = M.metric()
    sage: g[0,0], g[1,1], g[2,2] = -1, 1, 1
    sage: g.display()
    g = -dw⊗dw + dx⊗dx + dy⊗dy

With this, the declaration the ambient manifold and its foliation parametrized
by `t` is finished, and calculations can be performed.

The first step is always to find a chart adapted to the foliation. This is done
by the method "adapted_chart"::

    sage: T = N.adapted_chart(); T
    [Chart (M, (rh_M, th_M, t_M))]

``T`` contains a new chart defined on `M`. By default, the coordinate names
are constructed from the names of the submanifold coordinates and the foliation
parameter indexed by the name of the ambient manifold. By this can be
customized, see
:meth:`~sage.manifolds.topological_submanifold.TopologicalSubmanifold.adapted_chart`.

One can check that the adapted chart has been added to `M`'s atlas, along with
some coordinates changes::

    sage: M.atlas()
    [Chart (M, (w, x, y)), Chart (M, (rh_M, th_M, t_M))]
    sage: len(M.coord_changes())
    2

Let us compute the induced metric (or first fundamental form)::

    sage: gamma = N.induced_metric()  # long time
    sage: gamma.display()  # long time
    gamma = b^2 drh⊗drh + b^2*sinh(rh)^2 dth⊗dth
    sage: gamma[:]  # long time
    [           b^2              0]
    [             0 b^2*sinh(rh)^2]
    sage: gamma[1,1]  # long time
    b^2

the normal vector::

    sage: N.normal().display()  # long time
    n = sqrt(b^2 + x^2 + y^2)/b ∂/∂w + x/b ∂/∂x + y/b ∂/∂y

Check that the hypersurface is indeed spacelike, i.e. that its normal is
timelike::

    sage: N.ambient_metric()(N.normal(), N.normal()).display()  # long time
    g(n,n): M → ℝ
       (w, x, y) ↦ -1
       (rh_M, th_M, t_M) ↦ -1

The lapse function is::

    sage: N.lapse().display()  # long time
    N: M → ℝ
       (w, x, y) ↦ sqrt(b^2 + x^2 + y^2)/b
       (rh_M, th_M, t_M) ↦ cosh(rh_M)

while the shift vector is::

    sage: N.shift().display()  # long time
    beta = -(x^2 + y^2)/b^2 ∂/∂w - sqrt(b^2 + x^2 + y^2)*x/b^2 ∂/∂x
     - sqrt(b^2 + x^2 + y^2)*y/b^2 ∂/∂y

The extrinsic curvature (or second fundamental form) as a tensor field on the
ambient manifold::

    sage: N.ambient_extrinsic_curvature()[:] # long time
    [                                 -(x^2 + y^2)/b^3 (b^2*x + x^3 + x*y^2)/(sqrt(b^2 + x^2 + y^2)*b^3) (y^3 + (b^2 + x^2)*y)/(sqrt(b^2 + x^2 + y^2)*b^3)]
    [                      sqrt(b^2 + x^2 + y^2)*x/b^3                                  -(b^2 + x^2)/b^3                                          -x*y/b^3]
    [                      sqrt(b^2 + x^2 + y^2)*y/b^3                                          -x*y/b^3                                  -(b^2 + y^2)/b^3]

The extrinsic curvature as a tensor field on the submanifold::

    sage: N.extrinsic_curvature()[:] # long time
    [           -b             0]
    [            0 -b*sinh(rh)^2]


AUTHORS:

- Florentin Jaffredo (2018): initial version
- Eric Gourgoulhon (2018-2019): add documentation
- Matthias Koeppe (2021): open subsets of submanifolds

REFERENCES:

- \B. O'Neill : *Semi-Riemannian Geometry* [ONe1983]_
- \J. M. Lee : *Riemannian Manifolds* [Lee1997]_

"""

# *****************************************************************************
#  Copyright (C) 2018      Florentin Jaffredo <florentin.jaffredo@polytechnique.edu>
#  Copyright (C) 2018-2019 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#  Copyright (C) 2021      Matthias Koeppe <mkoeppe@math.ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from sage.manifolds.differentiable.pseudo_riemannian import \
    PseudoRiemannianManifold
from sage.manifolds.differentiable.differentiable_submanifold import \
    DifferentiableSubmanifold
from sage.rings.infinity import infinity
from sage.matrix.constructor import matrix
from sage.functions.other import factorial
from sage.symbolic.ring import SR
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer
from queue import Queue


class PseudoRiemannianSubmanifold(PseudoRiemannianManifold,
                                  DifferentiableSubmanifold):
    r"""
    Pseudo-Riemannian submanifold.

    An *embedded (resp. immersed) submanifold of a pseudo-Riemannian manifold*
    `(M,g)` is an embedded (resp. immersed) submanifold `N` of `M` as a
    differentiable manifold such that pull back of the metric tensor `g` via
    the embedding (resp. immersion) endows `N` with the structure of a
    pseudo-Riemannian manifold.

    INPUT:

    - ``n`` -- positive integer; dimension of the submanifold
    - ``name`` -- string; name (symbol) given to the submanifold
    - ``ambient`` -- (default: ``None``) pseudo-Riemannian manifold `M` in
      which the submanifold is embedded (or immersed). If ``None``, it is set
      to ``self``
    - ``metric_name`` -- (default: ``None``) string; name (symbol) given to the
      metric; if ``None``, ``'gamma'`` is used
    - ``signature`` -- (default: ``None``) signature `S` of the metric as a
      single integer: `S = n_+ - n_-`, where `n_+` (resp. `n_-`) is the
      number of positive terms (resp. number of negative terms) in any
      diagonal writing of the metric components; if ``signature`` is not
      provided, `S` is set to the submanifold's dimension (Riemannian
      signature)
    - ``base_manifold`` -- (default: ``None``) if not ``None``, must be a
      differentiable manifold; the created object is then an open subset of
      ``base_manifold``
    - ``diff_degree`` -- (default: ``infinity``) degree of differentiability
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the submanifold; if none is provided, it is set to ``name``
    - ``metric_latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the metric; if none is provided, it is set to ``metric_name`` if
      the latter is not ``None`` and to ``r'\gamma'`` otherwise
    - ``start_index`` -- (default: 0) integer; lower value of the range of
      indices used for "indexed objects" on the submanifold, e.g. coordinates
      in a chart
    - ``category`` -- (default: ``None``) to specify the category; if ``None``,
      ``Manifolds(RR).Differentiable()`` (or ``Manifolds(RR).Smooth()``
      if ``diff_degree`` = ``infinity``) is assumed (see the category
      :class:`~sage.categories.manifolds.Manifolds`)
    - ``unique_tag`` -- (default: ``None``) tag used to force the construction
      of a new object when all the other arguments have been used previously
      (without ``unique_tag``, the
      :class:`~sage.structure.unique_representation.UniqueRepresentation`
      behavior inherited from
      :class:`~sage.manifolds.subset.ManifoldSubset`, via
      :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`
      and :class:`~sage.manifolds.manifold.TopologicalManifold`,
      would return the previously constructed object corresponding to these
      arguments).

    EXAMPLES:

    Let `N` be a 2-dimensional submanifold of a 3-dimensional Riemannian
    manifold `M`::

        sage: M = Manifold(3, 'M', structure ="Riemannian")
        sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
        sage: N
        2-dimensional Riemannian submanifold N immersed in the 3-dimensional
         Riemannian manifold M
        sage: CM.<x,y,z> = M.chart()
        sage: CN.<u,v> = N.chart()

    Let us define a 1-dimension foliation indexed by `t`. The inverse map is
    needed in order to compute the adapted chart in the ambient manifold::

        sage: t = var('t')
        sage: phi = N.diff_map(M, {(CN,CM):[u, v, t+u^2+v^2]}); phi
        Differentiable map from the 2-dimensional Riemannian submanifold N
         immersed in the 3-dimensional Riemannian manifold M to the
         3-dimensional Riemannian manifold M
        sage: phi_inv = M.diff_map(N,{(CM, CN): [x,y]})
        sage: phi_inv_t = M.scalar_field({CM: z-x^2-y^2})

    `\phi` can then be declared as an embedding `N\to M`::

        sage: N.set_embedding(phi, inverse=phi_inv, var=t,
        ....:                 t_inverse={t: phi_inv_t})

    The foliation can also be used to find new charts on the ambient manifold
    that are adapted to the foliation, ie in which the expression of the
    immersion is trivial. At the same time, the appropriate coordinate changes
    are computed::

        sage: N.adapted_chart()
        [Chart (M, (u_M, v_M, t_M))]
        sage: len(M.coord_changes())
        2

    .. SEEALSO::

        :mod:`~sage.manifolds.manifold` and
        :mod:`~sage.manifolds.differentiable.differentiable_submanifold`

   """
    def __init__(self, n, name, ambient=None, metric_name=None,
                 signature=None, base_manifold=None, diff_degree=infinity,
                 latex_name=None, metric_latex_name=None, start_index=0,
                 category=None, unique_tag=None):
        r"""
        Construct a pseudo-Riemannian submanifold.

        TESTS::

            sage: M = Manifold(3, 'M', structure='Lorentzian')
            sage: N = Manifold(2, 'N', ambient=M, structure='Riemannian')
            sage: N
            2-dimensional Riemannian submanifold N immersed in the
             3-dimensional Lorentzian manifold M
            sage: phi = N.diff_map(M)
            sage: N.set_embedding(phi)
            sage: N
            2-dimensional Riemannian submanifold N embedded in the
             3-dimensional Lorentzian manifold M
            sage: S = Manifold(2, 'S', latex_name=r"\Sigma", ambient=M,
            ....:              structure="Riemannian", start_index=1)
            sage: latex(S)
            \Sigma
            sage: S.start_index()
            1

        ::

            sage: M = Manifold(5, 'M', structure='pseudo-Riemannian',
            ....:              signature=1)
            sage: N = Manifold(4, 'N', ambient=M,
            ....:              structure='pseudo-Riemannian', signature=0)
            sage: N
            4-dimensional pseudo-Riemannian submanifold N immersed in the
             5-dimensional pseudo-Riemannian manifold M

        """
        if metric_name is None:
            metric_name = 'gamma'
            metric_latex_name = r'\gamma'
        PseudoRiemannianManifold.__init__(self, n, name=name,
                                          metric_name=metric_name,
                                          signature=signature,
                                          base_manifold=base_manifold,
                                          diff_degree=diff_degree,
                                          latex_name=latex_name,
                                          metric_latex_name=metric_latex_name,
                                          start_index=start_index,
                                          category=category)
        self._init_immersion(ambient=ambient)
        self._difft = None
        self._gradt = None
        self._normal = None
        self._lapse = None
        self._shift = None
        self._first_fundamental_form = None
        self._ambient_first_fundamental_form = None
        self._second_fundamental_form = None
        self._ambient_second_fundamental_form = None
        self._ambient_metric = None
        self._projector = None
        self._gauss_curvature = None
        self._principal_directions = {}
        self._principal_curvatures = {}
        self._mean_curvature = None
        self._shape_operator = None
        self._sgn = 1 if ambient._structure.name == "Riemannian" else -1

    def _repr_(self):
        r"""
        Return a string representation of the submanifold.

        If no ambient manifold is specified, the submanifold is considered
        as a manifold.

        TESTS::

            sage: M = Manifold(3, 'M', structure="Lorentzian")
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: N
            2-dimensional Riemannian submanifold N immersed in the
             3-dimensional Lorentzian manifold M
            sage: phi = N.diff_map(M)
            sage: N.set_embedding(phi)
            sage: N
            2-dimensional Riemannian submanifold N embedded in the
             3-dimensional Lorentzian manifold M

        """
        if self is not self._manifold:
            return "Open subset {} of the {}".format(self._name, self._manifold)
        if self._ambient is None:
            return super(PseudoRiemannianManifold, self).__repr__()
        if self._embedded:
            return "{}-dimensional {} submanifold {} embedded in the {}".format(
                self._dim, self._structure.name, self._name, self._ambient)
        return "{}-dimensional {} submanifold {} immersed in the {}".format(
                self._dim, self._structure.name, self._name, self._ambient)

    def open_subset(self, name, latex_name=None, coord_def={}, supersets=None):
        r"""
        Create an open subset of ``self``.

        An open subset is a set that is (i) included in the manifold and (ii)
        open with respect to the manifold's topology. It is a differentiable
        manifold by itself. Moreover, equipped with the restriction of the
        manifold metric to itself, it is a pseudo-Riemannian manifold.

        As ``self`` is a submanifold of its ambient manifold,
        the new open subset is also considered a submanifold of that.
        Hence the returned object is an instance of
        :class:`PseudoRiemannianSubmanifold`.

        INPUT:

        - ``name`` -- name given to the open subset
        - ``latex_name`` --  (default: ``None``) LaTeX symbol to denote the
          subset; if none is provided, it is set to ``name``
        - ``coord_def`` -- (default: {}) definition of the subset in
          terms of coordinates; ``coord_def`` must a be dictionary with keys
          charts in the manifold's atlas and values the symbolic expressions
          formed by the coordinates to define the subset.
        - ``supersets`` -- (default: only ``self``) list of sets that the
          new open subset is a subset of

        OUTPUT:

        - instance of :class:`PseudoRiemannianSubmanifold` representing the
          created open subset

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure="Riemannian")
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian"); N
            2-dimensional Riemannian submanifold N immersed in the
             3-dimensional Riemannian manifold M
            sage: S = N.subset('S'); S
            Subset S of the
             2-dimensional Riemannian submanifold N immersed in the
              3-dimensional Riemannian manifold M
            sage: O = N.subset('O', is_open=True); O  # indirect doctest
            Open subset O of the
             2-dimensional Riemannian submanifold N immersed in the
              3-dimensional Riemannian manifold M

            sage: phi = N.diff_map(M)
            sage: N.set_embedding(phi)
            sage: N
            2-dimensional Riemannian submanifold N embedded in the
             3-dimensional Riemannian manifold M
            sage: S = N.subset('S'); S
            Subset S of the
             2-dimensional Riemannian submanifold N embedded in the
              3-dimensional Riemannian manifold M
            sage: O = N.subset('O', is_open=True); O  # indirect doctest
            Open subset O of the
             2-dimensional Riemannian submanifold N embedded in the
              3-dimensional Riemannian manifold M

        """
        resu = PseudoRiemannianSubmanifold(self._dim, name,
                                           ambient=self._ambient,
                                           metric_name=self._metric_name,
                                           signature=self._metric_signature,
                                           base_manifold=self._manifold,
                                           diff_degree=self._diff_degree,
                                           latex_name=latex_name,
                                           metric_latex_name=self._metric_latex_name,
                                           start_index=self._sindex)
        if supersets is None:
            supersets = [self]
        for superset in supersets:
            superset._init_open_subset(resu, coord_def=coord_def)
        return resu

    def ambient_metric(self):
        r"""
        Return the metric of the ambient manifold.

        OUTPUT:

        - the metric of the ambient manifold

        EXAMPLES::

            sage: M.<x,y,z> = EuclideanSpace()
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: N.ambient_metric()
            Riemannian metric g on the Euclidean space E^3
            sage: N.ambient_metric().display()
            g = dx⊗dx + dy⊗dy + dz⊗dz
            sage: N.ambient_metric() is M.metric()
            True

        """
        if self._ambient_metric is None:
            self._ambient_metric = self._ambient.metric()
        return self._ambient_metric

    def first_fundamental_form(self):
        r"""
        Return the first fundamental form of the submanifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        OUTPUT:

        - the first fundamental form, as an instance of
          :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`

        EXAMPLES:

        A sphere embedded in Euclidean space::

            sage: M.<x,y,z> = EuclideanSpace()
            sage: N = Manifold(2, 'N', ambient=M, structure='Riemannian')
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r', domain='real')
            sage: assume(r>0)
            sage: E = M.cartesian_coordinates()
            sage: phi = N.diff_map(M, {(C,E): [r*sin(th)*cos(ph),
            ....:                              r*sin(th)*sin(ph),
            ....:                              r*cos(th)]})
            sage: N.set_embedding(phi)
            sage: N.first_fundamental_form()  # long time
            Riemannian metric gamma on the 2-dimensional Riemannian
             submanifold N embedded in the Euclidean space E^3
            sage: N.first_fundamental_form()[:]  # long time
            [          r^2             0]
            [            0 r^2*sin(th)^2]

        An alias is ``induced_metric``::

            sage: N.induced_metric()[:]  # long time
            [          r^2             0]
            [            0 r^2*sin(th)^2]

        By default, the first fundamental form is named ``gamma``, but this
        can be customized by means of the argument ``metric_name`` when
        declaring the submanifold::

            sage: P = Manifold(1, 'P', ambient=M, structure='Riemannian',
            ....:              metric_name='g')
            sage: CP.<t> = P.chart()
            sage: F = P.diff_map(M, [t, 2*t, 3*t])
            sage: P.set_embedding(F)
            sage: P.induced_metric()
            Riemannian metric g on the 1-dimensional Riemannian submanifold P
             embedded in the Euclidean space E^3
            sage: P.induced_metric().display()
            g = 14 dt⊗dt

        """
        if self._first_fundamental_form is None:
            self._first_fundamental_form = super().metric()
            self._first_fundamental_form.set(
                               self._immersion.pullback(self.ambient_metric()))
        return self._first_fundamental_form

    induced_metric = first_fundamental_form

    def metric(self, name=None, signature=None, latex_name=None,
               dest_map=None):
        r"""
        Return the induced metric (first fundamental form) or define a new
        metric tensor on the submanifold.

        A new (uninitialzed) metric is returned only if the argument ``name``
        is provided and differs from the metric name declared at the
        construction of the submanifold; otherwise, the first fundamental
        form is returned.

        INPUT:

        - ``name`` -- (default: ``None``) name given to the metric; if ``name``
          is ``None`` or equals the metric name declared when constructing
          the submanifold, the first fundamental form is returned (see
          :meth:`first_fundamental_form`)
        - ``signature`` -- (default: ``None``; ignored if ``name`` is ``None``)
          signature `S` of the metric as a single integer: `S = n_+ - n_-`,
          where `n_+` (resp. `n_-`) is the number of positive terms (resp.
          number of negative terms) in any diagonal writing of the metric
          components; if ``signature`` is not provided, `S` is set to the
          submanifold's dimension (Riemannian signature)
        - ``latex_name`` -- (default: ``None``; ignored if ``name`` is ``None``)
          LaTeX symbol to denote the metric; if ``None``, it is formed from
          ``name``
        - ``dest_map`` -- (default: ``None``; ignored if ``name`` is ``None``)
          instance of
          class :class:`~sage.manifolds.differentiable.diff_map.DiffMap`
          representing the destination map `\Phi:\ U \rightarrow M`, where `U`
          is the current submanifold; if ``None``, the identity map is assumed
          (case of a metric tensor field *on* `U`)

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`

        EXAMPLES:

        Induced metric on a straight line of the Euclidean plane::

            sage: M.<x,y> = EuclideanSpace()
            sage: N = Manifold(1, 'N', ambient=M, structure='Riemannian')
            sage: CN.<t> = N.chart()
            sage: F = N.diff_map(M, [t, 2*t])
            sage: N.set_embedding(F)
            sage: N.metric()
            Riemannian metric gamma on the 1-dimensional Riemannian
             submanifold N embedded in the Euclidean plane E^2
            sage: N.metric().display()
            gamma = 5 dt⊗dt

        Setting the argument ``name`` to that declared while constructing
        the submanifold (default = ``'gamma'``) yields the same result::

            sage: N.metric(name='gamma') is N.metric()
            True

        while using a different name allows one to define a new metric on the
        submanifold::

            sage: h = N.metric(name='h'); h
            Riemannian metric h on the 1-dimensional Riemannian submanifold N
             embedded in the Euclidean plane E^2
            sage: h[0, 0] = 1  # initialization
            sage: h.display()
            h = dt⊗dt

        """
        if name is None or name == self._metric_name:
            return self.first_fundamental_form()
        return super().metric(name=name, signature=signature,
                              latex_name=latex_name, dest_map=dest_map)

    @cached_method
    def difft(self):
        r"""
        Return the differential of the scalar field on the ambient manifold
        representing the first parameter of the foliation associated to
        ``self``.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        OUTPUT:

        - 1-form field on the ambient manifold

        EXAMPLES:

        Foliation of the Euclidean 3-space by 2-spheres parametrized by their
        radii::

            sage: M.<x,y,z> = EuclideanSpace()
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r', domain='real')
            sage: assume(r>0)
            sage: E = M.cartesian_coordinates()
            sage: phi = N.diff_map(M, {(C,E): [r*sin(th)*cos(ph),
            ....:                              r*sin(th)*sin(ph),
            ....:                              r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E: sqrt(x^2+y^2+z^2)})
            sage: N.set_embedding(phi, inverse=phi_inv, var=r,
            ....:                 t_inverse={r: phi_inv_r})
            sage: N.difft()
            1-form dr on the Euclidean space E^3
            sage: N.difft().display()
            dr = x/sqrt(x^2 + y^2 + z^2) dx + y/sqrt(x^2 + y^2 + z^2) dy +
             z/sqrt(x^2 + y^2 + z^2) dz

        """
        if self._dim_foliation == 0:
            raise ValueError("A foliation is needed to "
                             "perform this calculation")
        self._difft = self._t_inverse[self._var[0]].differential()
        self._difft.set_name("d" + self._var[0]._repr_(),
                             r"\mathrm{d}" + self._var[0]._latex_())
        return self._difft

    @cached_method
    def gradt(self):
        r"""
        Return the gradient of the scalar field on the ambient manifold
        representing the first parameter of the foliation associated to
        ``self``.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        OUTPUT:

        - vector field on the ambient manifold

        EXAMPLES:

        Foliation of the Euclidean 3-space by 2-spheres parametrized by their
        radii::

            sage: M.<x,y,z> = EuclideanSpace()
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r', domain='real')
            sage: assume(r>0)
            sage: E = M.cartesian_coordinates()
            sage: phi = N.diff_map(M, {(C,E): [r*sin(th)*cos(ph),
            ....:                              r*sin(th)*sin(ph),
            ....:                              r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E: sqrt(x^2+y^2+z^2)})
            sage: N.set_embedding(phi, inverse=phi_inv, var=r,
            ....:                 t_inverse={r: phi_inv_r})
            sage: N.gradt()
            Vector field grad(r) on the Euclidean space E^3
            sage: N.gradt().display()
            grad(r) = x/sqrt(x^2 + y^2 + z^2) e_x + y/sqrt(x^2 + y^2 + z^2) e_y
             + z/sqrt(x^2 + y^2 + z^2) e_z

        """
        if self._dim_foliation == 0:
            raise ValueError("A foliation is needed to perform "
                             "this calculation")
        param = self._var[0]
        self._gradt = self._t_inverse[param].gradient()
        self._gradt.set_name("grad({})".format(param),
                             r"\mathrm{grad}\left(" + param._latex_()
                             + r"\right)")
        return self._gradt

    @cached_method
    def normal(self):
        r"""
        Return a normal unit vector to the submanifold.

        If a foliation is defined, it is used to compute the gradient of the
        foliation parameter and then the normal vector. If not, the normal
        vector is computed using the following formula:

        .. MATH::

            n = \vec{*}(\mathrm{d}x_0\wedge\mathrm{d}x_1\wedge\cdots
            \wedge\mathrm{d}x_{n-1})

        where the star stands for the Hodge dual operator and the wedge for the
        exterior product.

        This formula does not always define a proper vector field when
        multiple charts overlap, because of the arbitrariness of the direction
        of the normal vector. To avoid this problem, the method ``normal()``
        considers the graph defined by the atlas of the submanifold and the
        changes of coordinates, and only calculate the normal vector once by
        connected component. The expression is then propagate by restriction,
        continuation, or change of coordinates using a breadth-first
        exploration of the graph.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        OUTPUT:

        - vector field on the ambient manifold (case of a foliation) or along
          the submanifold with values in the ambient manifold (case of a
          single submanifold)

        EXAMPLES:

        Foliation of the Euclidean 3-space by 2-spheres parametrized by their
        radii::

            sage: M.<x,y,z> = EuclideanSpace()
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r', domain='real') # foliation parameter
            sage: assume(r>0)
            sage: E = M.cartesian_coordinates()
            sage: phi = N.diff_map(M, {(C,E): [r*sin(th)*cos(ph),
            ....:                              r*sin(th)*sin(ph),
            ....:                              r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E: sqrt(x^2+y^2+z^2)})
            sage: N.set_embedding(phi, inverse=phi_inv, var=r,
            ....:                 t_inverse={r: phi_inv_r})
            sage: T = N.adapted_chart()
            sage: N.normal()  # long time
            Vector field n on the Euclidean space E^3
            sage: N.normal().display()  # long time
            n = x/sqrt(x^2 + y^2 + z^2) e_x + y/sqrt(x^2 + y^2 + z^2) e_y
             + z/sqrt(x^2 + y^2 + z^2) e_z

        Or in spherical coordinates::

            sage: N.normal().display(T[0].frame(),T[0])  # long time
            n = ∂/∂r_E3

        Let us now consider a sphere of constant radius, i.e. not assumed to be
        part of a foliation, in stereographic coordinates::

            sage: M.<X,Y,Z> = EuclideanSpace()
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: U = N.open_subset('U')
            sage: V = N.open_subset('V')
            sage: N.declare_union(U, V)
            sage: stereoN.<x,y> = U.chart()
            sage: stereoS.<xp,yp> = V.chart("xp:x' yp:y'")
            sage: stereoN_to_S = stereoN.transition_map(stereoS,
            ....:                                 (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                                 intersection_name='W',
            ....:                                 restrictions1= x^2+y^2!=0,
            ....:                                 restrictions2= xp^2+yp^2!=0)
            sage: stereoS_to_N = stereoN_to_S.inverse()
            sage: W = U.intersection(V)
            sage: stereoN_W = stereoN.restrict(W)
            sage: stereoS_W = stereoS.restrict(W)
            sage: A = W.open_subset('A', coord_def={stereoN_W: (y!=0, x<0),
            ....:                                   stereoS_W: (yp!=0, xp<0)})
            sage: spher.<the,phi> = A.chart(r'the:(0,pi):\theta phi:(0,2*pi):\phi')
            sage: stereoN_A = stereoN_W.restrict(A)
            sage: spher_to_stereoN = spher.transition_map(stereoN_A,
            ....:                              (sin(the)*cos(phi)/(1-cos(the)),
            ....:                               sin(the)*sin(phi)/(1-cos(the))))
            sage: spher_to_stereoN.set_inverse(2*atan(1/sqrt(x^2+y^2)),
            ....:                              atan2(-y,-x)+pi)
            Check of the inverse coordinate transformation:
              the == 2*arctan(sqrt(-cos(the) + 1)/sqrt(cos(the) + 1))  **failed**
              phi == pi + arctan2(sin(phi)*sin(the)/(cos(the) - 1),
                                  cos(phi)*sin(the)/(cos(the) - 1))  **failed**
              x == x  *passed*
              y == y  *passed*
            NB: a failed report can reflect a mere lack of simplification.
            sage: stereoN_to_S_A = stereoN_to_S.restrict(A)
            sage: spher_to_stereoS = stereoN_to_S_A * spher_to_stereoN
            sage: stereoS_to_N_A = stereoN_to_S.inverse().restrict(A)
            sage: stereoS_to_spher = spher_to_stereoN.inverse() * stereoS_to_N_A
            sage: E = M.cartesian_coordinates()
            sage: phi = N.diff_map(M, {(stereoN, E): [2*x/(1+x^2+y^2),
            ....:                                    2*y/(1+x^2+y^2),
            ....:                                    (x^2+y^2-1)/(1+x^2+y^2)],
            ....:                   (stereoS, E): [2*xp/(1+xp^2+yp^2),
            ....:                                  2*yp/(1+xp^2+yp^2),
            ....:                                 (1-xp^2-yp^2)/(1+xp^2+yp^2)]},
            ....:                  name='Phi', latex_name=r'\Phi')
            sage: N.set_embedding(phi)

        The method ``normal()`` now returns a tensor field along ``N``::

            sage: n = N.normal()  # long time
            sage: n  # long time
            Vector field n along the 2-dimensional Riemannian submanifold N
             embedded in the Euclidean space E^3 with values on the Euclidean
             space E^3

        Let us check that the choice of orientation is coherent on the two top
        frames::

            sage: n.restrict(V).display(format_spec=spher)  # long time
            n = -cos(phi)*sin(the) e_X - sin(phi)*sin(the) e_Y - cos(the) e_Z
            sage: n.restrict(U).display(format_spec=spher)  # long time
            n = -cos(phi)*sin(the) e_X - sin(phi)*sin(the) e_Y - cos(the) e_Z

        """
        if self._dim_foliation != 0:    # case of a foliation
            self._normal = self._sgn * self.lapse() * self.gradt()
            self._normal.set_name("n")
            return self._normal
        # case of no foliation:
        max_frame = self._ambient.default_frame().along(self._immersion)
        self._normal = self.multivector_field(self._ambient._dim - self._dim,
                                              name="n",
                                              dest_map=self._immersion)

        # an auxiliary function:
        def calc_normal(chart):
            """
            Calculate the normal vector field according to the formula in the
            documentation in a given chart.
            """
            eps = self.ambient_metric().volume_form(self._dim).along(
                self._immersion).restrict(chart.domain())
            args = list(range(self._dim)) + [eps] + list(range(self._dim))
            r = self.irange()
            n_form = self._immersion.restrict(chart.domain()).pushforward(
                chart.frame()[next(r)]).down(
                self.ambient_metric().along(self._immersion).restrict(
                    chart.domain()))
            for i in r:
                n_form = n_form.wedge(
                    self._immersion.restrict(chart.domain()).pushforward(
                        chart.frame()[i]).down(
                        self.ambient_metric().along(
                            self._immersion).restrict(
                            chart.domain())))
            n_comp = (n_form.contract(*args) / factorial(self._dim)).contract(
                self.ambient_metric().inverse().along(self._immersion))
            if self._ambient._dim - self._dim == 1:
                n_comp = n_comp / n_comp.norm(self.ambient_metric())

            norm_rst = self._normal.restrict(chart.domain())
            norm_rst.add_comp(max_frame.restrict(chart.domain()))[:] = n_comp[:]
            self._normal.add_comp_by_continuation(max_frame, chart.domain(),
                                                  chart)

        # start breadth-first graph exploration
        marked = set()
        f = Queue()

        for v in self.top_charts():
            if v not in marked:
                f.put(v)
                calc_normal(v)  # initial calculus
                marked.add(v)
                while not f.empty():
                    v = f.get()
                    # for each neighbor:
                    for vp in self.atlas():
                        # case restriction
                        if vp in v._subcharts and vp not in marked:
                            f.put(vp)
                            self._normal.restrict(vp.domain())
                            marked.add(vp)

                        # case continuation
                        if vp in v._supercharts and vp not in marked:
                            f.put(vp)
                            self._normal.add_comp_by_continuation(
                                max_frame.restrict(vp.domain()), v.domain(), vp)
                            marked.add(vp)

                        # case coordinates change
                        if (v, vp) in self.coord_changes() and vp not in marked:
                            f.put(vp)
                            self._normal.comp(max_frame, vp)
                            marked.add(vp)

        # Going up from each top_chart to the full manifold :
        for v in self.top_charts():
            self._normal.add_expr_from_subdomain(max_frame, v.domain())

        return self._normal

    def ambient_first_fundamental_form(self):
        r"""
        Return the first fundamental form of the submanifold as a tensor of the
        ambient manifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        OUTPUT:

        - (0,2) tensor field on the ambient manifold describing the induced
          metric before projection on the submanifold

        EXAMPLES:

        A unit circle embedded in the Euclidean plane::

            sage: M.<X,Y> = EuclideanSpace()
            sage: N = Manifold(1, 'N', ambient=M, structure="Riemannian")
            sage: U = N.open_subset('U')
            sage: V = N.open_subset('V')
            sage: N.declare_union(U,V)
            sage: stereoN.<x> = U.chart()
            sage: stereoS.<y> = V.chart()
            sage: stereoN_to_S = stereoN.transition_map(stereoS, (4/x),
            ....:                   intersection_name='W',
            ....:                   restrictions1=x!=0, restrictions2=y!=0)
            sage: stereoS_to_N = stereoN_to_S.inverse()
            sage: E = M.cartesian_coordinates()
            sage: phi = N.diff_map(M,
            ....:         {(stereoN, E): [1/sqrt(1+x^2/4), x/2/sqrt(1+x^2/4)],
            ....:          (stereoS, E): [1/sqrt(1+4/y^2), 2/y/sqrt(1+4/y^2)]})
            sage: N.set_embedding(phi)
            sage: N.ambient_first_fundamental_form()
            Tensor field gamma of type (0,2) along the 1-dimensional Riemannian
             submanifold N embedded in the Euclidean plane E^2 with values on
             the Euclidean plane E^2
            sage: N.ambient_first_fundamental_form()[:]
            [ x^2/(x^2 + 4) -2*x/(x^2 + 4)]
            [-2*x/(x^2 + 4)    4/(x^2 + 4)]

        An alias is ``ambient_induced_metric``::

            sage: N.ambient_induced_metric()[:]
            [ x^2/(x^2 + 4) -2*x/(x^2 + 4)]
            [-2*x/(x^2 + 4)    4/(x^2 + 4)]

        """
        if self._ambient._dim - self._dim != 1:
            raise NotImplementedError("ambient_first_fundamental_form() is "
                                      "implemented only for hypersurfaces")
        if self._ambient_first_fundamental_form is None:
            g = self.ambient_metric()
            if self._dim_foliation == 0:  # case no foliation
                g = g.along(self._immersion)
            self._ambient_first_fundamental_form = g - self._sgn * g.contract(
                self.normal()) * g.contract(self.normal())
            self._ambient_first_fundamental_form.set_name("gamma", r"\gamma")
        return self._ambient_first_fundamental_form

    ambient_induced_metric = ambient_first_fundamental_form

    @cached_method
    def lapse(self):
        r"""
        Return the lapse function of the foliation.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        OUTPUT:

        - the lapse function, as a scalar field on the ambient manifold

        EXAMPLES:

        Foliation of the Euclidean 3-space by 2-spheres parametrized by their
        radii::

            sage: M.<x,y,z> = EuclideanSpace()
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r', domain='real') # foliation parameter
            sage: assume(r>0)
            sage: E = M.cartesian_coordinates()
            sage: phi = N.diff_map(M, {(C,E): [r*sin(th)*cos(ph),
            ....:                              r*sin(th)*sin(ph),
            ....:                              r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E: sqrt(x^2+y^2+z^2)})
            sage: N.set_embedding(phi, inverse=phi_inv, var=r,
            ....:                 t_inverse={r: phi_inv_r})
            sage: T = N.adapted_chart()
            sage: N.lapse()
            Scalar field N on the Euclidean space E^3
            sage: N.lapse().display()
            N: E^3 → ℝ
               (x, y, z) ↦ 1
               (th_E3, ph_E3, r_E3) ↦ 1

        """
        if self._dim_foliation == 0:
            raise ValueError("A foliation is needed "
                             "to perform this calculation")
        self._lapse = 1 / (self._sgn * self.ambient_metric()(
            self.gradt(), self.gradt())).sqrt()
        self._lapse.set_name("N")
        return self._lapse

    @cached_method
    def shift(self):
        r"""
        Return the shift vector associated with the first adapted chart of the
        foliation.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        OUTPUT:

        - shift vector field on the ambient manifold

        EXAMPLES:

        Foliation of the Euclidean 3-space by 2-spheres parametrized by their
        radii::

            sage: M.<x,y,z> = EuclideanSpace()
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r', domain='real') # foliation parameter
            sage: assume(r>0)
            sage: E = M.cartesian_coordinates()
            sage: phi = N.diff_map(M, {(C,E): [r*sin(th)*cos(ph),
            ....:                              r*sin(th)*sin(ph),
            ....:                              r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E: sqrt(x^2+y^2+z^2)})
            sage: N.set_embedding(phi, inverse=phi_inv, var=r,
            ....:                 t_inverse={r: phi_inv_r})
            sage: T = N.adapted_chart()
            sage: N.shift()  # long time
            Vector field beta on the Euclidean space E^3
            sage: N.shift().display()  # long time
            beta = 0

        """
        if self._dim_foliation == 0:
            raise ValueError("A foliation is needed "
                             "to perform this calculation")
        sia = self._ambient._sindex
        self._shift = self._adapted_charts[0].frame()[self._dim + sia]\
            - self.lapse() * self.normal()
        self._shift.set_name("beta", r"\beta")
        return self._shift

    def ambient_second_fundamental_form(self):
        r"""
        Return the second fundamental form of the submanifold as a tensor field
        on the ambient manifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        OUTPUT:

        - (0,2) tensor field on the ambient manifold equal to the second
          fundamental form once orthogonally projected onto the submanifold

        EXAMPLES:

        A unit circle embedded in the Euclidean plane::

            sage: M.<X,Y> = EuclideanSpace()
            sage: N = Manifold(1, 'N', ambient=M, structure="Riemannian")
            sage: U = N.open_subset('U')
            sage: V = N.open_subset('V')
            sage: N.declare_union(U,V)
            sage: stereoN.<x> = U.chart()
            sage: stereoS.<y> = V.chart()
            sage: stereoN_to_S = stereoN.transition_map(stereoS, (4/x),
            ....:                   intersection_name='W',
            ....:                   restrictions1=x!=0, restrictions2=y!=0)
            sage: stereoS_to_N = stereoN_to_S.inverse()
            sage: E = M.cartesian_coordinates()
            sage: phi = N.diff_map(M,
            ....:         {(stereoN, E): [1/sqrt(1+x^2/4), x/2/sqrt(1+x^2/4)],
            ....:          (stereoS, E): [1/sqrt(1+4/y^2), 2/y/sqrt(1+4/y^2)]})
            sage: N.set_embedding(phi)
            sage: N.ambient_second_fundamental_form()  # long time
            Field of symmetric bilinear forms K along the 1-dimensional
             Riemannian submanifold N embedded in the Euclidean plane E^2 with
             values on the Euclidean plane E^2
            sage: N.ambient_second_fundamental_form()[:] # long time
            [-x^2/(x^2 + 4)  2*x/(x^2 + 4)]
            [ 2*x/(x^2 + 4)   -4/(x^2 + 4)]

        An alias is ``ambient_extrinsic_curvature``::

            sage: N.ambient_extrinsic_curvature()[:]  # long time
            [-x^2/(x^2 + 4)  2*x/(x^2 + 4)]
            [ 2*x/(x^2 + 4)   -4/(x^2 + 4)]

        """
        if self._ambient._dim - self._dim != 1:
            raise ValueError("ambient_second_fundamental_form is defined only "
                             "for hypersurfaces")
        if self._ambient_second_fundamental_form is None:
            if self._dim_foliation == 0:
                self._ambient_second_fundamental_form = self.tensor_field(0, 2,
                                        sym=[(0, 1)], dest_map=self._immersion)
                k = self.second_fundamental_form()
                g = self.ambient_metric().along(self._immersion)
                max_frame = self._ambient.default_frame().along(self._immersion)
                for chart in self.top_charts():
                    pf = [self._immersion.restrict(chart.domain()).pushforward(
                        chart.frame()[i]) for i in self.irange()]
                    for i in range(self._dim):
                        pf[i] = pf[i] / g(pf[i], pf[i])
                    gam_rst = sum(
                        g.restrict(chart.domain()).contract(pf[i]) *
                        g.restrict(chart.domain()).contract(pf[j]) *
                        self.scalar_field({chart: k.comp(chart.frame())[:][i, j]})
                        for i in range(self._dim) for j in range(self._dim))
                    gam_rst._sym = [(0, 1)]
                    self._ambient_second_fundamental_form.set_restriction(gam_rst)

                charts = iter(self.top_charts())
                self._ambient_second_fundamental_form.add_comp_by_continuation(
                    max_frame, next(charts).domain())
                for chart in charts:
                    self._ambient_second_fundamental_form.add_expr_from_subdomain(
                        max_frame, chart.domain())
            else:
                nab = self.ambient_metric().connection('nabla', r'\nabla')
                self._ambient_second_fundamental_form = \
                    -self.ambient_metric().contract(nab(self.normal())) \
                    - nab(self.normal()).contract(self.normal())\
                    .contract(self.ambient_metric())\
                    * self.normal().contract(self.ambient_metric())
            self._ambient_second_fundamental_form.set_name("K")
        return self._ambient_second_fundamental_form

    ambient_extrinsic_curvature = ambient_second_fundamental_form

    def second_fundamental_form(self):
        r"""
        Return the second fundamental form of the submanifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        OUTPUT:

        - the second fundamental form, as a symmetric tensor field of type
          (0,2) on the submanifold

        EXAMPLES:

        A unit circle embedded in the Euclidean plane::

            sage: M.<X,Y> = EuclideanSpace()
            sage: N = Manifold(1, 'N', ambient=M, structure="Riemannian")
            sage: U = N.open_subset('U')
            sage: V = N.open_subset('V')
            sage: N.declare_union(U,V)
            sage: stereoN.<x> = U.chart()
            sage: stereoS.<y> = V.chart()
            sage: stereoN_to_S = stereoN.transition_map(stereoS, (4/x),
            ....:                   intersection_name='W',
            ....:                   restrictions1=x!=0, restrictions2=y!=0)
            sage: stereoS_to_N = stereoN_to_S.inverse()
            sage: E = M.cartesian_coordinates()
            sage: phi = N.diff_map(M,
            ....:         {(stereoN, E): [1/sqrt(1+x^2/4), x/2/sqrt(1+x^2/4)],
            ....:          (stereoS, E): [1/sqrt(1+4/y^2), 2/y/sqrt(1+4/y^2)]})
            sage: N.set_embedding(phi)
            sage: N.second_fundamental_form()  # long time
            Field of symmetric bilinear forms K on the 1-dimensional Riemannian
             submanifold N embedded in the Euclidean plane E^2
            sage: N.second_fundamental_form().display()  # long time
            K = -4/(x^4 + 8*x^2 + 16) dx⊗dx

        An alias is ``extrinsic_curvature``::

            sage: N.extrinsic_curvature().display()  # long time
            K = -4/(x^4 + 8*x^2 + 16) dx⊗dx

        An example with a non-Euclidean ambient metric::

            sage: M = Manifold(2, 'M', structure='Riemannian')
            sage: N = Manifold(1, 'N', ambient=M, structure='Riemannian',
            ....:              start_index=1)
            sage: CM.<x,y> = M.chart()
            sage: CN.<u> = N.chart()
            sage: g = M.metric()
            sage: g[0, 0], g[1, 1] = 1, 1/(1 + y^2)^2
            sage: phi = N.diff_map(M, (u, u))
            sage: N.set_embedding(phi)
            sage: N.second_fundamental_form()
            Field of symmetric bilinear forms K on the 1-dimensional Riemannian
             submanifold N embedded in the 2-dimensional Riemannian manifold M
            sage: N.second_fundamental_form().display()
            K = 2*sqrt(u^4 + 2*u^2 + 2)*u/(u^6 + 3*u^4 + 4*u^2 + 2) du⊗du

        """
        if self._ambient._dim - self._dim != 1:
            raise ValueError("second_fundamental_form is defined only for"
                             + " hypersurfaces")
        if self._second_fundamental_form is None:
            resu = self.vector_field_module().tensor((0, 2), name='K',
                                                      sym=[(0, 1)])
            if self._dim_foliation != 0:
                inverse_subs = {v: k for k, v in self._subs[0].items()}
                asff = self.ambient_second_fundamental_form()
                adapted_chart = self._adapted_charts[0]
                dsi = self._ambient._sindex - self._sindex
                for i in self.irange():
                    for j in self.irange(start=i):
                        resu[i, j] = asff[adapted_chart.frame(),
                                          i + dsi, j + dsi,
                                          adapted_chart].expr().subs(inverse_subs)
            else:
                nab = self.ambient_metric().connection('nabla', r'\nabla')
                n = self.normal()

                for chart in self.atlas():
                    gamma_n = matrix(SR, self._dim + 1, self._dim + 1)
                    subs = dict(zip(self._ambient.default_chart()[:],
                                    self._immersion.expression(chart)))
                    for i in range(self._dim + 1):
                        for j in range(self._dim + 1):
                            Gam_ij = [nab[self._ambient.frames()[0],
                                          :][i][j][k].expr().subs(subs)
                                      for k in range(self._dim + 1)]
                            gamma_n[i, j] = chart.simplify(sum(
                                Gam_ij[k] *
                                n.restrict(chart.domain()).comp(
                                  n.restrict(chart.domain())._fmodule.bases()[0])
                                [:][k].expr() for k in range(self._dim + 1)))
                    dXdu = self._immersion.differential_functions(chart)
                    dNdu = matrix(SR, self._dim + 1, self._dim)
                    for i in range(self._dim + 1):
                        for j in range(self._dim):
                            dNdu[i, j] = n.restrict(chart.domain()).comp(
                                n.restrict(chart.domain())._fmodule.bases()[0])[:,
                                chart][i].diff(chart[:][j]).expr()
                    g = self.ambient_metric().along(
                        self._immersion.restrict(chart.domain())).restrict(
                        chart.domain())[:, chart]
                    K = dXdu.transpose() * g * (dNdu + gamma_n * dXdu)
                    si = self._sindex
                    for i in self.irange():
                        for j in self.irange(i):  # since K is symmetric
                            resu[chart.frame(), i, j, chart] = chart.simplify(
                                                      K[i - si, j - si].expr())

            self._second_fundamental_form = resu
        return self._second_fundamental_form

    extrinsic_curvature = second_fundamental_form

    @cached_method
    def projector(self):
        r"""
        Return the orthogonal projector onto the submanifold.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        OUTPUT:

        - the orthogonal projector onto the submanifold, as tensor field of
          type (1,1) on the ambient manifold

        EXAMPLES:

        Foliation of the Euclidean 3-space by 2-spheres parametrized by their
        radii::

            sage: M.<x,y,z> = EuclideanSpace()
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r', domain='real') # foliation parameter
            sage: assume(r>0)
            sage: E = M.cartesian_coordinates()
            sage: phi = N.diff_map(M, {(C,E): [r*sin(th)*cos(ph),
            ....:                              r*sin(th)*sin(ph),
            ....:                              r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E: sqrt(x^2+y^2+z^2)})
            sage: N.set_embedding(phi, inverse=phi_inv, var=r,
            ....:                 t_inverse={r: phi_inv_r})
            sage: T = N.adapted_chart()

        The orthogonal projector onto ``N`` is a type-(1,1) tensor field on
        ``M``::

            sage: N.projector()  # long time
            Tensor field gamma of type (1,1) on the Euclidean space E^3

        Check that the orthogonal projector applied to the normal vector is
        zero::

            sage: N.projector().contract(N.normal()).display()  # long time
            0

        """
        if self._ambient._dim - self._dim != 1:
            raise NotImplementedError("projector() is implemented only for "
                                      "hypersurfaces")
        g = self.ambient_metric().inverse()
        if self._dim_foliation == 0:
            g = g.along(self._immersion)

        self._projector = self.ambient_first_fundamental_form().contract(0, g)
        self._projector.set_name("gamma", r"\vec{\gamma}")
        return self._projector

    def project(self, tensor):
        r"""
        Return the orthogonal projection of a tensor field onto the submanifold.

        INPUT:

        - ``tensor`` -- any tensor field to be projected onto the submanifold.
          If no foliation is provided, must be a tensor field along the
          submanifold.

        OUTPUT:

        - orthogonal projection of ``tensor`` onto the submanifold, as a
          tensor field of the *ambient* manifold

        EXAMPLES:

        Foliation of the Euclidean 3-space by 2-spheres parametrized by their
        radii::

            sage: M.<x,y,z> = EuclideanSpace()
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r', domain='real') # foliation parameter
            sage: assume(r>0)
            sage: E = M.cartesian_coordinates()
            sage: phi = N.diff_map(M, {(C,E): [r*sin(th)*cos(ph),
            ....:                              r*sin(th)*sin(ph),
            ....:                              r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E: sqrt(x^2+y^2+z^2)})
            sage: N.set_embedding(phi, inverse=phi_inv, var=r,
            ....:                 t_inverse={r: phi_inv_r})
            sage: T = N.adapted_chart()

        Let us perform the projection of the ambient metric and check that it
        is equal to the first fundamental form::

            sage: pg = N.project(M.metric()); pg  # long time
            Tensor field of type (0,2) on the Euclidean space E^3
            sage: pg == N.ambient_first_fundamental_form()  # long time
            True

        Note that the output of ``project()`` is not cached.

        """
        if self._ambient._dim - self._dim != 1:
            raise NotImplementedError("project() is implemented only for "
                                      "hypersurfaces")
        resu = tensor.copy()
        resu.set_name(tensor._name + "_" + self._name,
                      r"{" + tensor._latex_() + r"}_{" + self._latex_() + r"}")
        for i in range(tensor.tensor_type()[0]):
            resu = self.projector().contract(1, resu, i)
        for i in range(tensor.tensor_type()[1]):
            resu = self.projector().contract(0, resu, i)
        return resu

    def mixed_projection(self, tensor, indices=0):
        r"""
        Return de n+1 decomposition of a tensor on the submanifold and the
        normal vector.

        The n+1 decomposition of a tensor of rank `k` can be obtained by
        contracting each index either with the normal vector or the projection
        operator of the submanifold (see
        :meth:`~sage.manifolds.differentiable.pseudo_riemannian_submanifold.PseudoRiemannianSubmanifold.projector`).

        INPUT:

        - ``tensor`` -- any tensor field, eventually along the submanifold if
          no foliation is provided.
        - ``indices`` -- (default: ``0``) list of integers containing the
          indices on which the projection is made on the normal vector.
          By default, all projections are made on the submanifold. If
          an integer `n` is provided, the `n` first contractions are made with
          the normal vector, all the other ones with the orthogonal projection
          operator.

        OUTPUT:

        - tensor field of rank `k`-``len(indices)``

        EXAMPLES:

        Foliation of the Euclidean 3-space by 2-spheres parametrized by their
        radii::

            sage: M.<x,y,z> = EuclideanSpace()
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r', domain='real') # foliation parameter
            sage: assume(r>0)
            sage: E = M.cartesian_coordinates()
            sage: phi = N.diff_map(M, {(C,E): [r*sin(th)*cos(ph),
            ....:                              r*sin(th)*sin(ph),
            ....:                              r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E: sqrt(x^2+y^2+z^2)})
            sage: N.set_embedding(phi, inverse=phi_inv, var=r,
            ....:                 t_inverse={r: phi_inv_r})
            sage: T = N.adapted_chart()

        If ``indices`` is not specified, the mixed projection of the ambient
        metric coincides with the first fundamental form::

            sage: g = M.metric()
            sage: gpp = N.mixed_projection(g); gpp  # long time
            Tensor field of type (0,2) on the Euclidean space E^3
            sage: gpp == N.ambient_first_fundamental_form()  # long time
            True

        The other non-redundant projections are::

            sage: gnp = N.mixed_projection(g, [0]); gnp  # long time
            1-form on the Euclidean space E^3

        and::

            sage: gnn = N.mixed_projection(g, [0,1]); gnn
            Scalar field on the Euclidean space E^3

        which is constant and equal to 1 (the norm of the unit normal vector)::

            sage: gnn.display()
            E^3 → ℝ
            (x, y, z) ↦ 1
            (th_E3, ph_E3, r_E3) ↦ 1

        """
        if self._ambient._dim - self._dim != 1:
            raise NotImplementedError("mixed_projection() is implemented only "
                                      "for hypersurfaces")
        if isinstance(indices, (Integer, int)):
            indices = list(range(indices))

        if len(indices) > tensor.tensor_rank():
            raise ValueError("Too much contractions")

        g = self.ambient_metric()
        if self._dim_foliation == 0:
            g = g.along(self._immersion)

        multiprojector = 1
        k = tensor.tensor_rank()      # order of the tensor
        kp = 2 * k - len(indices)       # order of the multiprojector
        for i in range(tensor.tensor_type()[1]):
            if i in indices:
                multiprojector = multiprojector * self.normal()
            else:
                multiprojector = multiprojector * self.projector()
        for i in range(tensor.tensor_type()[0]):
            if i in indices:
                multiprojector = multiprojector * self.normal().contract(g)
            else:
                multiprojector = multiprojector * self.projector()
        args = list(range(kp - tensor.tensor_type()[0], kp)) + list(range(
                tensor.tensor_type()[1])) + [tensor] + list(range(k))
        return multiprojector.contract(*args)

    @cached_method
    def gauss_curvature(self):
        r"""
        Return the Gauss curvature of the submanifold.

        The *Gauss curvature* is the product or the principal curvatures, or
        equivalently the determinant of the projection operator.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        OUTPUT:

        - the Gauss curvature as a scalar field on the submanifold

        EXAMPLES:

        A unit circle embedded in the Euclidean plane::

            sage: M.<X,Y> = EuclideanSpace()
            sage: N = Manifold(1, 'N', ambient=M, structure="Riemannian")
            sage: U = N.open_subset('U')
            sage: V = N.open_subset('V')
            sage: N.declare_union(U,V)
            sage: stereoN.<x> = U.chart()
            sage: stereoS.<y> = V.chart()
            sage: stereoN_to_S = stereoN.transition_map(stereoS, (4/x),
            ....:                   intersection_name='W',
            ....:                   restrictions1=x!=0, restrictions2=y!=0)
            sage: stereoS_to_N = stereoN_to_S.inverse()
            sage: E = M.cartesian_coordinates()
            sage: phi = N.diff_map(M,
            ....:         {(stereoN, E): [1/sqrt(1+x^2/4), x/2/sqrt(1+x^2/4)],
            ....:          (stereoS, E): [1/sqrt(1+4/y^2), 2/y/sqrt(1+4/y^2)]})
            sage: N.set_embedding(phi)
            sage: N.gauss_curvature()  # long time
            Scalar field on the 1-dimensional Riemannian submanifold N embedded
             in the Euclidean plane E^2
            sage: N.gauss_curvature().display()  # long time
            N → ℝ
            on U: x ↦ -1
            on V: y ↦ -1

        """
        if self._ambient._dim - self._dim != 1:
            raise ValueError("gauss_curvature is defined only for "
                             "hypersurfaces")
        a = self.shape_operator()
        self._gauss_curvature = self.scalar_field(
            {chart: a[chart.frame(), :, chart].determinant()
             for chart in self.top_charts()})
        return self._gauss_curvature

    @cached_method
    def principal_directions(self, chart):
        r"""
        Return the principal directions of the submanifold.

        The *principal directions* are the eigenvectors of the projection
        operator. The result is formatted as a list of pairs
        (eigenvector, eigenvalue).

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        INPUT:

        - ``chart`` --  chart in which the principal directions are to be
          computed

        OUTPUT:

        - list of pairs (vector field, scalar field) representing the
          principal directions and the associated principal curvatures

        EXAMPLES:

        A unit circle embedded in the Euclidean plane::

            sage: M.<X,Y> = EuclideanSpace()
            sage: N = Manifold(1, 'N', ambient=M, structure="Riemannian")
            sage: U = N.open_subset('U')
            sage: V = N.open_subset('V')
            sage: N.declare_union(U,V)
            sage: stereoN.<x> = U.chart()
            sage: stereoS.<y> = V.chart()
            sage: stereoN_to_S = stereoN.transition_map(stereoS, (4/x),
            ....:                   intersection_name='W',
            ....:                   restrictions1=x!=0, restrictions2=y!=0)
            sage: stereoS_to_N = stereoN_to_S.inverse()
            sage: E = M.cartesian_coordinates()
            sage: phi = N.diff_map(M,
            ....:         {(stereoN, E): [1/sqrt(1+x^2/4), x/2/sqrt(1+x^2/4)],
            ....:          (stereoS, E): [1/sqrt(1+4/y^2), 2/y/sqrt(1+4/y^2)]})
            sage: N.set_embedding(phi)
            sage: N.principal_directions(stereoN)  # long time
            [(Vector field e_0 on the 1-dimensional Riemannian submanifold N
              embedded in the Euclidean plane E^2, -1)]
            sage: N.principal_directions(stereoN)[0][0].display()  # long time
            e_0 = ∂/∂x

        """
        if self._ambient._dim - self._dim != 1:
            raise ValueError("principal directions is defined only for "
                             "hypersurfaces")
        a = self.shape_operator()
        pr_d = matrix(
            [[a[chart.frame(), :, chart][i, j].expr() for i in self.irange()]
             for j in self.irange()]).eigenvectors_right()
        res = []
        v = self.vector_field()
        counter = self.irange()
        for eigen_space in pr_d:
            for eigen_vector in eigen_space[1]:
                v[chart.frame(), :] = eigen_vector
                res.append((v.copy(), eigen_space[0]))
                res[-1][0].set_name("e_{}".format(next(counter)))
        self._principal_directions[chart] = res
        return res

    @cached_method
    def principal_curvatures(self, chart):
        r"""
        Return the principal curvatures of the submanifold.

        The *principal curvatures* are the eigenvalues of the projection
        operator. The resulting scalar fields are named ``k_i`` with the
        index ``i`` ranging from 0 to the submanifold dimension minus one.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        INPUT:

        - ``chart`` --  chart in which the principal curvatures are to be
          computed

        OUTPUT:

        - the principal curvatures, as a list of scalar fields on the
          submanifold

        EXAMPLES:

        A unit circle embedded in the Euclidean plane::

            sage: M.<X,Y> = EuclideanSpace()
            sage: N = Manifold(1, 'N', ambient=M, structure="Riemannian")
            sage: U = N.open_subset('U')
            sage: V = N.open_subset('V')
            sage: N.declare_union(U,V)
            sage: stereoN.<x> = U.chart()
            sage: stereoS.<y> = V.chart()
            sage: stereoN_to_S = stereoN.transition_map(stereoS, (4/x),
            ....:                   intersection_name='W',
            ....:                   restrictions1=x!=0, restrictions2=y!=0)
            sage: stereoS_to_N = stereoN_to_S.inverse()
            sage: E = M.cartesian_coordinates()
            sage: phi = N.diff_map(M,
            ....:         {(stereoN, E): [1/sqrt(1+x^2/4), x/2/sqrt(1+x^2/4)],
            ....:          (stereoS, E): [1/sqrt(1+4/y^2), 2/y/sqrt(1+4/y^2)]})
            sage: N.set_embedding(phi)
            sage: N.principal_curvatures(stereoN)  # long time
            [Scalar field k_0 on the 1-dimensional Riemannian submanifold N
             embedded in the Euclidean plane E^2]
            sage: N.principal_curvatures(stereoN)[0].display()  # long time
            k_0: N → ℝ
            on U: x ↦ -1
            on W: y ↦ -1

        """
        if self._ambient._dim - self._dim != 1:
            raise ValueError("principal_curvatures is defined only for "
                             "hypersurfaces")
        a = self.shape_operator()
        res = matrix(
            [[a[chart.frame(), :, chart][i, j].expr() for i in self.irange()]
             for j in self.irange()]).eigenvalues()
        counter = self.irange()
        for i in range(self._dim):
            res[i] = self.scalar_field({chart: res[i]},
                                       name="k_{}".format(next(counter)))
        self._principal_curvatures[chart] = res
        return res

    @cached_method
    def mean_curvature(self):
        r"""
        Return the mean curvature of the submanifold.

        The *mean curvature* is the arithmetic mean of the principal curvatures,
        or equivalently the trace of the projection operator.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        OUTPUT:

        - the mean curvature, as a scalar field on the submanifold

        EXAMPLES:

       A unit circle embedded in the Euclidean plane::

            sage: M.<X,Y> = EuclideanSpace()
            sage: N = Manifold(1, 'N', ambient=M, structure="Riemannian")
            sage: U = N.open_subset('U')
            sage: V = N.open_subset('V')
            sage: N.declare_union(U,V)
            sage: stereoN.<x> = U.chart()
            sage: stereoS.<y> = V.chart()
            sage: stereoN_to_S = stereoN.transition_map(stereoS, (4/x),
            ....:                   intersection_name='W',
            ....:                   restrictions1=x!=0, restrictions2=y!=0)
            sage: stereoS_to_N = stereoN_to_S.inverse()
            sage: E = M.cartesian_coordinates()
            sage: phi = N.diff_map(M,
            ....:         {(stereoN, E): [1/sqrt(1+x^2/4), x/2/sqrt(1+x^2/4)],
            ....:          (stereoS, E): [1/sqrt(1+4/y^2), 2/y/sqrt(1+4/y^2)]})
            sage: N.set_embedding(phi)
            sage: N.mean_curvature()  # long time
            Scalar field on the 1-dimensional Riemannian submanifold N
             embedded in the Euclidean plane E^2
            sage: N.mean_curvature().display()  # long time
            N → ℝ
            on U: x ↦ -1
            on V: y ↦ -1

        """
        if self._ambient._dim - self._dim != 1:
            raise ValueError("mean_curvature is defined only for "
                             "hypersurfaces")
        self._shape_operator = self.scalar_field({chart: self._sgn * sum(
            self.principal_curvatures(chart)).expr(chart) / self._dim
                                                  for chart in
                                                  self.top_charts()})
        return self._shape_operator

    @cached_method
    def shape_operator(self):
        r"""
        Return the shape operator of the submanifold.

        The shape operator is equal to the second fundamental form with one of
        the indices upped.

        The result is cached, so calling this method multiple times always
        returns the same result at no additional cost.

        OUTPUT:

        - the shape operator, as a tensor field of type (1,1) on the
          submanifold

        EXAMPLES:

        A unit circle embedded in the Euclidean plane::

            sage: M.<X,Y> = EuclideanSpace()
            sage: N = Manifold(1, 'N', ambient=M, structure="Riemannian")
            sage: U = N.open_subset('U')
            sage: V = N.open_subset('V')
            sage: N.declare_union(U,V)
            sage: stereoN.<x> = U.chart()
            sage: stereoS.<y> = V.chart()
            sage: stereoN_to_S = stereoN.transition_map(stereoS, (4/x),
            ....:                   intersection_name='W',
            ....:                   restrictions1=x!=0, restrictions2=y!=0)
            sage: stereoS_to_N = stereoN_to_S.inverse()
            sage: E = M.cartesian_coordinates()
            sage: phi = N.diff_map(M,
            ....:         {(stereoN, E): [1/sqrt(1+x^2/4), x/2/sqrt(1+x^2/4)],
            ....:          (stereoS, E): [1/sqrt(1+4/y^2), 2/y/sqrt(1+4/y^2)]})
            sage: N.set_embedding(phi)
            sage: N.shape_operator()  # long time
            Tensor field of type (1,1) on the 1-dimensional Riemannian
             submanifold N embedded in the Euclidean plane E^2
            sage: N.shape_operator().display()  # long time
            -∂/∂x⊗dx

        """
        if self._ambient._dim - self._dim != 1:
            raise ValueError("shape_operator is defined only for "
                             "hypersurfaces")
        self._shape_operator = self.second_fundamental_form().contract(
                                               self.induced_metric().inverse())
        return self._shape_operator

    def clear_cache(self):
        r"""
        Reset all the cached functions and the derived quantities.

        Use this function if you modified the immersion (or embedding) of the
        submanifold. Note that when calling a calculus function after clearing,
        new Python objects will be created.

        EXAMPLES::

            sage: M.<x,y,z> = EuclideanSpace()
            sage: N = Manifold(2, 'N', ambient=M, structure="Riemannian")
            sage: C.<th,ph> = N.chart(r'th:(0,pi):\theta ph:(-pi,pi):\phi')
            sage: r = var('r', domain='real') # foliation parameter
            sage: assume(r>0)
            sage: E = M.cartesian_coordinates()
            sage: phi = N.diff_map(M, {(C,E): [r*sin(th)*cos(ph),
            ....:                              r*sin(th)*sin(ph),
            ....:                              r*cos(th)]})
            sage: phi_inv = M.diff_map(N, {(E,C): [arccos(z/r), atan2(y,x)]})
            sage: phi_inv_r = M.scalar_field({E: sqrt(x^2+y^2+z^2)})
            sage: N.set_embedding(phi, inverse=phi_inv, var=r,
            ....:                 t_inverse={r: phi_inv_r})
            sage: T = N.adapted_chart()
            sage: n = N.normal()
            sage: n is N.normal()
            True
            sage: N.clear_cache()
            sage: n is N.normal()
            False
            sage: n == N.normal()
            True

        """
        self.difft.clear_cache()
        self.gradt.clear_cache()
        self.normal.clear_cache()
        self.lapse.clear_cache()
        self.shift.clear_cache()
        self.projector.clear_cache()
        self.gauss_curvature.clear_cache()
        self.principal_directions.clear_cache()
        self.principal_curvatures.clear_cache()
        self.shape_operator.clear_cache()
        self._difft = None
        self._gradt = None
        self._normal = None
        self._lapse = None
        self._shift = None
        self._first_fundamental_form = None
        self._ambient_first_fundamental_form = None
        self._second_fundamental_form = None
        self._ambient_second_fundamental_form = None
        self._ambient_metric = None
        self._projector = None
        self._gauss_curvature = None
        self._principal_directions = {}
        self._principal_curvatures = {}
        self._mean_curvature = None
        self._shape_operator = None
