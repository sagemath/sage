r"""
Pseudo-Riemannian Manifolds

A *pseudo-Riemannian manifold* is a pair `(M,g)` where `M` is a real
differentiable manifold `M` (see
:class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`)
and `g` is a field of non-degenerate symmetric bilinear forms on `M`, which is
called the *metric tensor*, or simply the *metric* (see
:class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`).

Two important subcases are

- *Riemannian manifold*: the metric `g` is positive definite, i.e. its
  signature is `n = \dim M`;
- *Lorentzian manifold*: the metric `g` has signature `n-2` (positive
  convention) or `2-n` (negative convention).

On a pseudo-Riemannian manifold, one may use various standard
:mod:`~sage.manifolds.operators` acting on scalar and tensor
fields, like :func:`~sage.manifolds.operators.grad` or
:func:`~sage.manifolds.operators.div`.

All pseudo-Riemannian manifolds, whatever the metric signature, are implemented
via the class :class:`PseudoRiemannianManifold`.

.. RUBRIC:: Example 1: the sphere as a Riemannian manifold of dimension 2

We start by declaring `S^2` as a 2-dimensional Riemannian manifold::

    sage: M = Manifold(2, 'S^2', structure='Riemannian')
    sage: M
    2-dimensional Riemannian manifold S^2

We then cover `S^2` by two stereographic charts, from the North pole and from
the South pole respectively::

    sage: U = M.open_subset('U')
    sage: stereoN.<x,y> = U.chart()
    sage: V = M.open_subset('V')
    sage: stereoS.<u,v> = V.chart()
    sage: M.declare_union(U,V)
    sage: stereoN_to_S = stereoN.transition_map(stereoS,
    ....:                [x/(x^2+y^2), y/(x^2+y^2)], intersection_name='W',
    ....:                restrictions1= x^2+y^2!=0, restrictions2= u^2+v^2!=0)
    sage: W = U.intersection(V)
    sage: stereoN_to_S
    Change of coordinates from Chart (W, (x, y)) to Chart (W, (u, v))
    sage: stereoN_to_S.display()
    u = x/(x^2 + y^2)
    v = y/(x^2 + y^2)
    sage: stereoN_to_S.inverse().display()
    x = u/(u^2 + v^2)
    y = v/(u^2 + v^2)

We get the metric defining the Riemannian structure by::

    sage: g = M.metric()
    sage: g
    Riemannian metric g on the 2-dimensional Riemannian manifold S^2

At this stage, the metric `g` is defined as a Python object but there remains to
initialize it by setting its components with respect to the vector frames
associated with the stereographic coordinates. Let us begin with the frame
of chart ``stereoN``::

    sage: eU = stereoN.frame()
    sage: g[eU, 0, 0] = 4/(1 + x^2 + y^2)^2
    sage: g[eU, 1, 1] = 4/(1 + x^2 + y^2)^2

The metric components in the frame of chart ``stereoS`` are obtained by
continuation of the expressions found in `W = U\cap V` from the known
change-of-coordinate formulas::

    sage: eV = stereoS.frame()
    sage: g.add_comp_by_continuation(eV, W)

At this stage, the metric `g` is well defined in all `S^2`::

    sage: g.display(eU)
    g = 4/(x^2 + y^2 + 1)^2 dx⊗dx + 4/(x^2 + y^2 + 1)^2 dy⊗dy
    sage: g.display(eV)
    g = 4/(u^4 + v^4 + 2*(u^2 + 1)*v^2 + 2*u^2 + 1) du⊗du
     + 4/(u^4 + v^4 + 2*(u^2 + 1)*v^2 + 2*u^2 + 1) dv⊗dv

The expression in frame ``eV`` can be given a shape similar to that in frame
``eU``, by factorizing the components::

    sage: g[eV, 0, 0].factor()
    4/(u^2 + v^2 + 1)^2
    sage: g[eV, 1, 1].factor()
    4/(u^2 + v^2 + 1)^2
    sage: g.display(eV)
    g = 4/(u^2 + v^2 + 1)^2 du⊗du + 4/(u^2 + v^2 + 1)^2 dv⊗dv

Let us consider a scalar field `f` on `S^2`::

    sage: f = M.scalar_field({stereoN: 1/(1+x^2+y^2)}, name='f')
    sage: f.add_expr_by_continuation(stereoS, W)
    sage: f.display()
    f: S^2 → ℝ
    on U: (x, y) ↦ 1/(x^2 + y^2 + 1)
    on V: (u, v) ↦ (u^2 + v^2)/(u^2 + v^2 + 1)

The gradient of `f` (with respect to the metric `g`) is::

    sage: gradf = f.gradient()
    sage: gradf
    Vector field grad(f) on the 2-dimensional Riemannian manifold S^2
    sage: gradf.display(eU)
    grad(f) = -1/2*x ∂/∂x - 1/2*y ∂/∂y
    sage: gradf.display(eV)
    grad(f) = 1/2*u ∂/∂u + 1/2*v ∂/∂v

It is possible to write ``grad(f)`` instead of ``f.gradient()``, by importing
the standard differential operators of vector calculus::

    sage: from sage.manifolds.operators import *
    sage: grad(f) == gradf
    True

The Laplacian of `f`  (with respect to the metric `g`) is obtained either
as ``f.laplacian()`` or, thanks to the above import, as ``laplacian(f)``::

    sage: Df = laplacian(f)
    sage: Df
    Scalar field Delta(f) on the 2-dimensional Riemannian manifold S^2
    sage: Df.display()
    Delta(f): S^2 → ℝ
    on U: (x, y) ↦ (x^2 + y^2 - 1)/(x^2 + y^2 + 1)
    on V: (u, v) ↦ -(u^2 + v^2 - 1)/(u^2 + v^2 + 1)

Let us check the standard formula
`\Delta f = \mathrm{div}( \mathrm{grad}\,  f )`::

    sage: Df == div(gradf)
    True

Since each open subset of `S^2` inherits the structure of a Riemannian
manifold, we can get the metric on it via the method ``metric()``::

    sage: gU = U.metric()
    sage: gU
    Riemannian metric g on the Open subset U of the 2-dimensional Riemannian
     manifold S^2
    sage: gU.display()
    g = 4/(x^2 + y^2 + 1)^2 dx⊗dx + 4/(x^2 + y^2 + 1)^2 dy⊗dy

Of course, ``gU`` is nothing but the restriction of `g` to `U`::

    sage: gU is g.restrict(U)
    True

.. RUBRIC:: Example 2: Minkowski spacetime as a Lorentzian manifold of
  dimension 4

We start by declaring a 4-dimensional Lorentzian manifold `M`::

    sage: M = Manifold(4, 'M', structure='Lorentzian')
    sage: M
    4-dimensional Lorentzian manifold M

We define Minkowskian coordinates on `M`::

    sage: X.<t,x,y,z> = M.chart()

We construct the metric tensor by::

    sage: g = M.metric()
    sage: g
    Lorentzian metric g on the 4-dimensional Lorentzian manifold M

and initialize it to the Minkowskian value::

    sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1, 1, 1, 1
    sage: g.display()
    g = -dt⊗dt + dx⊗dx + dy⊗dy + dz⊗dz
    sage: g[:]
    [-1  0  0  0]
    [ 0  1  0  0]
    [ 0  0  1  0]
    [ 0  0  0  1]

We may check that the metric is flat, i.e. has a vanishing Riemann curvature
tensor::

    sage: g.riemann().display()
    Riem(g) = 0

A vector field on `M`::

    sage: u = M.vector_field(name='u')
    sage: u[0] = cosh(t)
    sage: u[1] = sinh(t)
    sage: u.display()
    u = cosh(t) ∂/∂t + sinh(t) ∂/∂x

The scalar square of `u` is::

    sage: s = u.dot(u); s
    Scalar field u.u on the 4-dimensional Lorentzian manifold M

Scalar products are taken with respect to the metric tensor::

    sage: u.dot(u) == g(u,u)
    True

`u` is a unit timelike vector, i.e. its scalar square is identically `-1`::

    sage: s.display()
    u.u: M → ℝ
       (t, x, y, z) ↦ -1
    sage: s.expr()
    -1

Let us consider a unit spacelike vector::

    sage: v = M.vector_field(name='v')
    sage: v[0] = sinh(t)
    sage: v[1] = cosh(t)
    sage: v.display()
    v = sinh(t) ∂/∂t + cosh(t) ∂/∂x
    sage: v.dot(v).display()
    v.v: M → ℝ
       (t, x, y, z) ↦ 1
    sage: v.dot(v).expr()
    1

`u` and `v` are orthogonal vectors with respect to Minkowski metric::

    sage: u.dot(v).display()
        u.v: M → ℝ
       (t, x, y, z) ↦ 0
    sage: u.dot(v).expr()
    0

The divergence of `u` is::

    sage: s = u.div(); s
    Scalar field div(u) on the 4-dimensional Lorentzian manifold M
    sage: s.display()
    div(u): M → ℝ
       (t, x, y, z) ↦ sinh(t)

while its d'Alembertian is::

    sage: Du = u.dalembertian(); Du
    Vector field Box(u) on the 4-dimensional Lorentzian manifold M
    sage: Du.display()
    Box(u) = -cosh(t) ∂/∂t - sinh(t) ∂/∂x

AUTHORS:

- Eric Gourgoulhon (2018): initial version

REFERENCES:

- \B. O'Neill : *Semi-Riemannian Geometry* [ONe1983]_
- \J. M. Lee : *Riemannian Manifolds* [Lee1997]_

"""

#*****************************************************************************
#       Copyright (C) 2018 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.infinity import infinity
from sage.manifolds.structure import (PseudoRiemannianStructure,
                                      RiemannianStructure, LorentzianStructure)
from sage.manifolds.differentiable.manifold import DifferentiableManifold

###############################################################################

class PseudoRiemannianManifold(DifferentiableManifold):
    r"""
    PseudoRiemannian manifold.

    A *pseudo-Riemannian manifold* is a pair `(M,g)` where `M` is a real
    differentiable manifold `M` (see
    :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`)
    and `g` is a field of non-degenerate symmetric bilinear forms on `M`, which
    is called the *metric tensor*, or simply the *metric* (see
    :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`).

    Two important subcases are

    - *Riemannian manifold*: the metric `g` is positive definite, i.e. its
      signature is `n = \dim M`;
    - *Lorentzian manifold*: the metric `g` has signature `n-2` (positive
      convention) or `2-n` (negative convention).

    INPUT:

    - ``n`` -- positive integer; dimension of the manifold
    - ``name`` -- string; name (symbol) given to the manifold
    - ``metric_name`` -- (default: ``None``) string; name (symbol) given to the
      metric; if ``None``, ``'g'`` is used
    - ``signature`` -- (default: ``None``) signature `S` of the metric as a
      single integer: `S = n_+ - n_-`, where `n_+` (resp. `n_-`) is the
      number of positive terms (resp. number of negative terms) in any
      diagonal writing of the metric components; if ``signature`` is not
      provided, `S` is set to the manifold's dimension (Riemannian
      signature)
    - ``base_manifold`` -- (default: ``None``) if not ``None``, must be a
      differentiable manifold; the created object is then an open subset of
      ``base_manifold``
    - ``diff_degree`` -- (default: ``infinity``) degree `k` of
      differentiability
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the manifold; if none is provided, it is set to ``name``
    - ``metric_latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the metric; if none is provided, it is set to ``metric_name``
    - ``start_index`` -- (default: 0) integer; lower value of the range of
      indices used for "indexed objects" on the manifold, e.g. coordinates
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

    Pseudo-Riemannian manifolds are constructed via the generic function
    :func:`~sage.manifolds.manifold.Manifold`, using the keyword
    ``structure``::

        sage: M = Manifold(4, 'M', structure='pseudo-Riemannian', signature=0)
        sage: M
        4-dimensional pseudo-Riemannian manifold M
        sage: M.category()
        Category of smooth manifolds over Real Field with 53 bits of precision

    The metric associated with ``M`` is::

        sage: M.metric()
        Pseudo-Riemannian metric g on the 4-dimensional pseudo-Riemannian
         manifold M
        sage: M.metric().signature()
        0
        sage: M.metric().tensor_type()
        (0, 2)

    Its value has to be initialized either by setting its components in various
    vector frames (see the above examples regarding the 2-sphere and Minkowski
    spacetime) or by making it equal to a given field of symmetric bilinear
    forms (see the method
    :meth:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric.set`
    of the metric class). Both methods are also covered in the
    documentation of method :meth:`metric` below.

    The metric object belongs to the class
    :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`::

        sage: isinstance(M.metric(), sage.manifolds.differentiable.metric.
        ....:                        PseudoRiemannianMetric)
        True

    See the documentation of this class for all operations available on
    metrics.

    The default name of the metric is ``g``; it can be customized::

        sage: M = Manifold(4, 'M', structure='pseudo-Riemannian',
        ....:              metric_name='gam', metric_latex_name=r'\gamma')
        sage: M.metric()
        Riemannian metric gam on the 4-dimensional Riemannian manifold M
        sage: latex(M.metric())
        \gamma

    A Riemannian manifold is constructed by the proper setting of the keyword
    ``structure``::

        sage: M = Manifold(4, 'M', structure='Riemannian'); M
        4-dimensional Riemannian manifold M
        sage: M.metric()
        Riemannian metric g on the 4-dimensional Riemannian manifold M
        sage: M.metric().signature()
        4

    Similarly, a Lorentzian manifold is obtained by::

        sage: M = Manifold(4, 'M', structure='Lorentzian'); M
        4-dimensional Lorentzian manifold M
        sage: M.metric()
        Lorentzian metric g on the 4-dimensional Lorentzian manifold M

    The default Lorentzian signature is taken to be positive::

        sage: M.metric().signature()
        2

    but one can opt for the negative convention via the keyword ``signature``::

        sage: M = Manifold(4, 'M', structure='Lorentzian', signature='negative')
        sage: M.metric()
        Lorentzian metric g on the 4-dimensional Lorentzian manifold M
        sage: M.metric().signature()
        -2

    """
    def __init__(self, n, name, metric_name=None, signature=None,
                 base_manifold=None, diff_degree=infinity, latex_name=None,
                 metric_latex_name=None, start_index=0, category=None,
                 unique_tag=None):
        r"""
        Construct a pseudo-Riemannian manifold.

        TESTS::

            sage: M = Manifold(4, 'M', structure='pseudo-Riemannian',
            ....:              signature=0)
            sage: M
            4-dimensional pseudo-Riemannian manifold M
            sage: type(M)
            <class 'sage.manifolds.differentiable.pseudo_riemannian.PseudoRiemannianManifold_with_category'>
            sage: X.<w,x,y,z> = M.chart()
            sage: M.metric()
            Pseudo-Riemannian metric g on the 4-dimensional pseudo-Riemannian manifold M
            sage: TestSuite(M).run()

        """
        if base_manifold and not isinstance(base_manifold, PseudoRiemannianManifold):
            raise TypeError("the argument 'base_manifold' must be a " +
                            "pseudo-Riemannian manifold")
        if signature is None or signature == n:
            structure = RiemannianStructure()
        elif signature == n-2 or signature == 2-n:
            structure = LorentzianStructure()
        else:
            structure = PseudoRiemannianStructure()
        DifferentiableManifold.__init__(self, n, name, 'real', structure,
                                        base_manifold=base_manifold,
                                        diff_degree=diff_degree,
                                        latex_name=latex_name,
                                        start_index=start_index,
                                        category=category)
        self._metric = None # to be initialized by metric()
        self._metric_signature = signature
        if metric_name is None:
            metric_name = 'g'
        elif not isinstance(metric_name, str):
            raise TypeError("{} is not a string".format(metric_name))
        self._metric_name = metric_name
        if metric_latex_name is None:
            self._metric_latex_name = self._metric_name
        else:
            if not isinstance(metric_latex_name, str):
                raise TypeError("{} is not a string".format(metric_latex_name))
            self._metric_latex_name = metric_latex_name

    def metric(self, name=None, signature=None, latex_name=None,
               dest_map=None):
        r"""
        Return the metric giving the pseudo-Riemannian structure to the
        manifold, or define a new metric tensor on the manifold.

        INPUT:

        - ``name`` -- (default: ``None``) name given to the metric; if
          ``name`` is ``None`` or matches the name of the metric defining the
          pseudo-Riemannian structure of ``self``, the latter metric is
          returned
        - ``signature`` -- (default: ``None``; ignored if ``name`` is ``None``)
          signature `S` of the metric as a single integer: `S = n_+ - n_-`,
          where `n_+` (resp. `n_-`) is the number of positive terms (resp.
          number of negative terms) in any diagonal writing of the metric
          components; if ``signature`` is not provided, `S` is set to the
          manifold's dimension (Riemannian signature)
        - ``latex_name`` -- (default: ``None``; ignored if ``name`` is ``None``)
          LaTeX symbol to denote the metric; if ``None``, it is formed from
          ``name``
        - ``dest_map`` -- (default: ``None``; ignored if ``name`` is ``None``)
          instance of
          class :class:`~sage.manifolds.differentiable.diff_map.DiffMap`
          representing the destination map `\Phi:\ U \rightarrow M`, where `U`
          is the current manifold; if ``None``, the identity map is assumed
          (case of a metric tensor field *on* `U`)

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`

        EXAMPLES:

        Metric of a 3-dimensional Riemannian manifold::

            sage: M = Manifold(3, 'M', structure='Riemannian', start_index=1)
            sage: X.<x,y,z> = M.chart()
            sage: g = M.metric(); g
            Riemannian metric g on the 3-dimensional Riemannian manifold M

        The metric remains to be initialized, for instance by setting its
        components in the coordinate frame associated to the chart ``X``::

            sage: g[1,1], g[2,2], g[3,3] = 1, 1, 1
            sage: g.display()
            g = dx⊗dx + dy⊗dy + dz⊗dz

        Alternatively, the metric can be initialized from a given field of
        nondegenerate symmetric bilinear forms; we may create the former
        object by::

            sage: X.coframe()
            Coordinate coframe (M, (dx,dy,dz))
            sage: dx, dy, dz = X.coframe()[1], X.coframe()[2], X.coframe()[3]
            sage: b = dx*dx + dy*dy + dz*dz
            sage: b
            Field of symmetric bilinear forms dx⊗dx+dy⊗dy+dz⊗dz on the
             3-dimensional Riemannian manifold M

        We then use the metric method
        :meth:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric.set`
        to make ``g`` being equal to ``b`` as a symmetric tensor field of
        type ``(0,2)``::

            sage: g.set(b)
            sage: g.display()
            g = dx⊗dx + dy⊗dy + dz⊗dz

        Another metric can be defined on ``M`` by specifying a metric name
        distinct from that chosen at the creation of the manifold (which
        is ``g`` by default, but can be changed thanks to the keyword
        ``metric_name`` in :func:`~sage.manifolds.manifold.Manifold`)::

            sage: h = M.metric('h'); h
            Riemannian metric h on the 3-dimensional Riemannian manifold M
            sage: h[1,1], h[2,2], h[3,3] = 1+y^2, 1+z^2, 1+x^2
            sage: h.display()
            h = (y^2 + 1) dx⊗dx + (z^2 + 1) dy⊗dy + (x^2 + 1) dz⊗dz

        The metric tensor ``h`` is distinct from the metric entering in the
        definition of the Riemannian manifold ``M``::

            sage: h is M.metric()
            False

        while we have of course::

            sage: g is M.metric()
            True

        Providing the same name as the manifold's default metric returns the
        latter::

            sage: M.metric('g') is M.metric()
            True

        In the present case (``M`` is diffeomorphic to `\RR^3`), we can even
        create a Lorentzian metric on ``M``::

            sage: h = M.metric('h', signature=1); h
            Lorentzian metric h on the 3-dimensional Riemannian manifold M

        """
        if name is None or name == self._metric_name:
            # Default metric associated with the manifold
            if self._metric is None:
                if self._manifold is not self and self._manifold._metric is not None:
                    # case of an open subset with a metric already defined on
                    # the ambient manifold:
                    self._metric = self._manifold._metric.restrict(self)
                else:
                    # creation from scratch:
                    self._metric = DifferentiableManifold.metric(self,
                                           self._metric_name,
                                           signature=self._metric_signature,
                                           latex_name=self._metric_latex_name)
            return self._metric
        # Metric distinct from the default one: it is created by the method
        # metric of the superclass for generic differentiable manifolds:
        return DifferentiableManifold.metric(self, name, signature=signature,
                                             latex_name=latex_name,
                                             dest_map=dest_map)

    def volume_form(self, contra=0):
        r"""
        Volume form (Levi-Civita tensor) `\epsilon` associated with ``self``.

        This assumes that ``self`` is an orientable manifold, with a
        preferred orientation; see
        :meth:`~sage.manifolds.differentiable.manifold.DifferentiableManifold.orientation`
        for details.

        The volume form `\epsilon` is a `n`-form (`n` being the manifold's
        dimension) such that, for any vector frame `(e_i)` that is orthonormal
        with respect to the metric of the pseudo-Riemannian manifold ``self``,

        .. MATH::

            \epsilon(e_1,\ldots,e_n) = \pm 1

        There are only two such `n`-forms, which are opposite of each other.
        The volume form `\epsilon` is selected as the one that returns `+1` for
        any right-handed vector frame with respect to the chosen orientation of
        ``self``.

        INPUT:

        - ``contra`` -- (default: 0) number of contravariant indices of the
          returned tensor

        OUTPUT:

        - if ``contra = 0`` (default value): the volume `n`-form `\epsilon`, as
          an instance of
          :class:`~sage.manifolds.differentiable.diff_form.DiffForm`
        - if ``contra = k``, with `1\leq k \leq n`, the tensor field of type
          (k,n-k) formed from `\epsilon` by raising the first k indices with
          the metric (see method
          :meth:`~sage.manifolds.differentiable.tensorfield.TensorField.up`);
          the output is then an instance of
          :class:`~sage.manifolds.differentiable.tensorfield.TensorField`, with
          the appropriate antisymmetries, or of the subclass
          :class:`~sage.manifolds.differentiable.multivectorfield.MultivectorField`
          if `k=n`

        EXAMPLES:

        Volume form of the Euclidean 3-space::

            sage: M = Manifold(3, 'M', structure='Riemannian', start_index=1)
            sage: X.<x,y,z> = M.chart()
            sage: g = M.metric()
            sage: g[1,1], g[2,2], g[3,3] = 1, 1, 1
            sage: eps = M.volume_form(); eps
            3-form eps_g on the 3-dimensional Riemannian manifold M
            sage: eps.display()
            eps_g = dx∧dy∧dz

        Raising the first index::

            sage: eps1 = M.volume_form(1); eps1
            Tensor field of type (1,2) on the 3-dimensional Riemannian
             manifold M
            sage: eps1.display()
            ∂/∂x⊗dy⊗dz - ∂/∂x⊗dz⊗dy - ∂/∂y⊗dx⊗dz + ∂/∂y⊗dz⊗dx + ∂/∂z⊗dx⊗dy
             - ∂/∂z⊗dy⊗dx
            sage: eps1.symmetries()
            no symmetry; antisymmetry: (1, 2)

        Raising the first and second indices::

            sage: eps2 = M.volume_form(2); eps2
            Tensor field of type (2,1) on the 3-dimensional Riemannian
             manifold M
            sage: eps2.display()
            ∂/∂x⊗∂/∂y⊗dz - ∂/∂x⊗∂/∂z⊗dy - ∂/∂y⊗∂/∂x⊗dz + ∂/∂y⊗∂/∂z⊗dx
             + ∂/∂z⊗∂/∂x⊗dy - ∂/∂z⊗∂/∂y⊗dx
            sage: eps2.symmetries()
            no symmetry; antisymmetry: (0, 1)

        Fully contravariant version::

            sage: eps3 = M.volume_form(3); eps3
            3-vector field on the 3-dimensional Riemannian manifold M
            sage: eps3.display()
            ∂/∂x∧∂/∂y∧∂/∂z

        """
        return self.metric().volume_form(contra=contra)

    def open_subset(self, name, latex_name=None, coord_def={}, supersets=None):
        r"""
        Create an open subset of ``self``.

        An open subset is a set that is (i) included in the manifold and (ii)
        open with respect to the manifold's topology. It is a differentiable
        manifold by itself. Moreover, equipped with the restriction of the
        manifold metric to itself, it is a pseudo-Riemannian manifold. Hence
        the returned object is an instance of
        :class:`PseudoRiemannianManifold`.

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

        - instance of :class:`PseudoRiemannianManifold` representing the
          created open subset

        EXAMPLES:

        Open subset of a 2-dimensional Riemannian manifold::

            sage: M = Manifold(2, 'M', structure='Riemannian')
            sage: X.<x,y> = M.chart()
            sage: U = M.open_subset('U', coord_def={X: x>0}); U
            Open subset U of the 2-dimensional Riemannian manifold M
            sage: type(U)
            <class 'sage.manifolds.differentiable.pseudo_riemannian.PseudoRiemannianManifold_with_category'>

        We initialize the metric of ``M``::

            sage: g = M.metric()
            sage: g[0,0], g[1,1] = 1, 1

        Then the metric on ``U`` is determined as the restriction of ``g`` to
        ``U``::

            sage: gU = U.metric(); gU
            Riemannian metric g on the Open subset U of the 2-dimensional Riemannian manifold M
            sage: gU.display()
            g = dx⊗dx + dy⊗dy
            sage: gU is g.restrict(U)
            True

        TESTS:

        Open subset created after the initialization of the metric::

            sage: V = M.open_subset('V', coord_def={X: x<0}); V
            Open subset V of the 2-dimensional Riemannian manifold M
            sage: gV = V.metric()
            sage: gV.display()
            g = dx⊗dx + dy⊗dy
            sage: gV is g.restrict(V)
            True

        """
        resu = PseudoRiemannianManifold(self._dim, name,
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
