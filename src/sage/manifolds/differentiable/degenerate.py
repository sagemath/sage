r"""
Degenerate manifolds
"""
# *****************************************************************************
#  Copyright (C) 2019 Hans Fotsing Tetsing <hans.fotsing@aims-cameroon.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# *****************************************************************************

from __future__ import annotations
from typing import TYPE_CHECKING

from sage.rings.infinity import infinity
from sage.manifolds.structure import DegenerateStructure
from sage.manifolds.differentiable.manifold import DifferentiableManifold

if TYPE_CHECKING:
    from sage.manifolds.differentiable.metric import DegenerateMetric

###############################################################################


class DegenerateManifold(DifferentiableManifold):
    r"""

    Degenerate Manifolds

    A *degenerate manifold* (or a *null manifold*) is a pair `(M,g)`
    where `M` is a real differentiable manifold  (see
    :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`)
    and `g` is a field of degenerate symmetric bilinear forms on `M` (see
    :class:`~sage.manifolds.differentiable.metric.DegenerateMetric`).

    INPUT:

    - ``n`` -- positive integer; dimension of the manifold
    - ``name`` -- string; name (symbol) given to the manifold
    - ``metric_name`` -- (default: ``None``) string; name (symbol) given to the
      metric; if ``None``, ``'g'`` is used
    - ``signature`` -- (default: ``None``) signature `S` of the metric as a
      tuple: `S = (n_+, n_-, n_0)`, where `n_+` (resp. `n_-`, resp. `n_0`) is the
      number of positive terms (resp. negative terms, resp. zero tems) in any
      diagonal writing of the metric components; if ``signature`` is not
      provided, `S` is set to `(ndim-1, 0, 1)`, being `ndim` the manifold's dimension
    - ``ambient`` -- (default: ``None``) if not ``None``, must be a
      differentiable manifold; the created object is then an open subset of
      ``ambient``
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

    A degenerate manifold is constructed via the generic function
    :func:`~sage.manifolds.manifold.Manifold`, using the keyword
    ``structure``::

        sage: M = Manifold(3, 'M', structure='degenerate_metric')
        sage: M
        3-dimensional degenerate_metric manifold M
        sage: M.parent()
        <class 'sage.manifolds.differentiable.degenerate.DegenerateManifold_with_category'>

    The metric associated with ``M`` is::

        sage: g = M.metric()
        sage: g
        degenerate metric g on the 3-dimensional degenerate_metric manifold M
        sage: g.signature()
        (0, 2, 1)

    Its value has to be initialized either by setting its components in various
    vector frames (see the above examples regarding the 2-sphere and Minkowski
    spacetime) or by making it equal to a given field of symmetric bilinear
    forms (see the method
    :meth:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric.set`
    of the metric class). Both methods are also covered in the
    documentation of method :meth:`metric` below.

    REFERENCES:

    - [DB1996]_
    - [DS2010]_
    """
    def __init__(self, n, name, metric_name=None, signature=None,
                 base_manifold=None, diff_degree=infinity, latex_name=None,
                 metric_latex_name=None, start_index=0, category=None,
                 unique_tag=None):
        r"""
        Construct a degenerate manifold.

        TESTS::

            sage: M = Manifold(3, 'M', structure='degenerate_metric')
            sage: M
            3-dimensional degenerate_metric manifold M
            sage: type(M)
            <class 'sage.manifolds.differentiable.degenerate.DegenerateManifold_with_category'>
            sage: X.<x,y,z> = M.chart()
            sage: M.metric()
            degenerate metric g on the 3-dimensional degenerate_metric manifold M
            sage: TestSuite(M).run()
        """
        if base_manifold and not isinstance(base_manifold, DegenerateManifold):
            raise TypeError("the argument 'base_manifold' must be a " +
                            "Degenerate manifold")
        structure = DegenerateStructure()
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
               dest_map=None) -> DegenerateMetric:
        r"""
        Return the metric giving the null manifold structure to the
        manifold, or define a new metric tensor on the manifold.

        INPUT:

        - ``name`` -- (default: ``None``) name given to the metric; if
          ``name`` is ``None`` or matches the name of the metric defining the
          null manifold structure of ``self``, the latter metric is
          returned
        - ``signature`` -- (default: ``None``; ignored if ``name`` is ``None``)
          signature `S` of the metric as a tuple: `S = (n_+, n_-, n_0)`,
          where `n_+` (resp. `n_-`, resp. `n_0`) is the number of positive
          terms (resp. negative terms, resp. zero tems) in any diagonal writing
          of the metric components; if ``signature`` is not provided, `S` is set
          to `(ndim-1, 0, 1)`, being `ndim` the manifold's dimension
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
          :class:`~sage.manifolds.differentiable.metric.DegenerateMetric`

        EXAMPLES:

        Metric of a 3-dimensional degenerate manifold::

            sage: M = Manifold(3, 'M', structure='degenerate_metric', start_index=1)
            sage: X.<x,y,z> = M.chart()
            sage: g = M.metric(); g
            degenerate metric g on the 3-dimensional degenerate_metric manifold M

        The metric remains to be initialized, for instance by setting its
        components in the coordinate frame associated to the chart ``X``::

            sage: g[1,1], g[2,2] = -1, 1
            sage: g.display()
            g = -dx⊗dx + dy⊗dy
            sage: g[:]
            [-1  0  0]
            [ 0  1  0]
            [ 0  0  0]

        Alternatively, the metric can be initialized from a given field of
        degenerate symmetric bilinear forms; we may create the former
        object by::

            sage: X.coframe()
            Coordinate coframe (M, (dx,dy,dz))
            sage: dx, dy = X.coframe()[1], X.coframe()[2]
            sage: b = dx*dx + dy*dy
            sage: b
            Field of symmetric bilinear forms dx⊗dx+dy⊗dy on the 3-dimensional
            degenerate_metric manifold M

        We then use the metric method
        :meth:`~sage.manifolds.differentiable.metric.DegenerateMetric.set`
        to make ``g`` being equal to ``b`` as a symmetric tensor field of
        type ``(0,2)``::

            sage: g.set(b)
            sage: g.display()
            g = dx⊗dx + dy⊗dy

        Another metric can be defined on ``M`` by specifying a metric name
        distinct from that chosen at the creation of the manifold (which
        is ``g`` by default, but can be changed thanks to the keyword
        ``metric_name`` in :func:`~sage.manifolds.manifold.Manifold`)::

            sage: h = M.metric('h'); h
            degenerate metric h on the 3-dimensional degenerate_metric manifold M
            sage: h[1,1], h[2,2], h[3,3] = 1+y^2, 1+z^2, 1+x^2
            sage: h.display()
            h = (y^2 + 1) dx⊗dx + (z^2 + 1) dy⊗dy + (x^2 + 1) dz⊗dz

        The metric tensor ``h`` is distinct from the metric entering in the
        definition of the degenerate manifold ``M``::

            sage: h is M.metric()
            False

        while we have of course::

            sage: g is M.metric()
            True

        Providing the same name as the manifold's default metric returns the
        latter::

            sage: M.metric('g') is M.metric()
            True
        """
        if signature is None:
            signature = self._metric_signature
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

    def open_subset(self, name, latex_name=None, coord_def={}):
        r"""
        Create an open subset of ``self``.

        An open subset is a set that is (i) included in the manifold and (ii)
        open with respect to the manifold's topology. It is a differentiable
        manifold by itself. Moreover, equipped with the restriction of the
        manifold metric to itself, it is a null manifold. Hence
        the returned object is an instance of
        :class:`DegenerateManifold`.

        INPUT:

        - ``name`` -- name given to the open subset
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          subset; if none is provided, it is set to ``name``
        - ``coord_def`` -- (default: {}) definition of the subset in
          terms of coordinates; ``coord_def`` must a be dictionary with keys
          charts in the manifold's atlas and values the symbolic expressions
          formed by the coordinates to define the subset.

        OUTPUT:

        - instance of :class:`DegenerateManifold` representing the
          created open subset

        EXAMPLES:

        Open subset of a 3-dimensional degenerate manifold::

            sage: M = Manifold(3, 'M', structure='degenerate_metric', start_index=1)
            sage: X.<x,y,z> = M.chart()
            sage: U = M.open_subset('U', coord_def={X: [x>0, y>0]}); U
            Open subset U of the 3-dimensional degenerate_metric manifold M
            sage: type(U)
            <class 'sage.manifolds.differentiable.degenerate.DegenerateManifold_with_category'>

        We initialize the metric of ``M``::

            sage: g = M.metric()
            sage: g[1,1], g[2,2] = -1, 1

        Then the metric on ``U`` is determined as the restriction of ``g`` to
        ``U``::

            sage: gU = U.metric(); gU
            degenerate metric g on the Open subset U of the 3-dimensional
            degenerate_metric manifold M
            sage: gU.display()
            g = -dx⊗dx + dy⊗dy
            sage: gU is g.restrict(U)
            True

        TESTS:

        Open subset created after the initialization of the metric::

            sage: V = M.open_subset('V', coord_def={X: x<0}); V
            Open subset V of the 3-dimensional degenerate_metric manifold M
            sage: gV = V.metric()
            sage: gV.display()
            g = -dx⊗dx + dy⊗dy
            sage: gV is g.restrict(V)
            True
        """
        resu = DegenerateManifold(self._dim, name,
                                        metric_name=self._metric_name,
                                        signature=self._metric_signature,
                                        base_manifold=self._manifold,
                                        diff_degree=self._diff_degree,
                                        latex_name=latex_name,
                                        metric_latex_name=self._metric_latex_name,
                                        start_index=self._sindex)
        resu._calculus_method = self._calculus_method
        resu._supersets.update(self._supersets)
        for sd in self._supersets:
            sd._subsets.add(resu)
        self._top_subsets.add(resu)
        # Charts on the result from the coordinate definition:
        for chart, restrictions in coord_def.items():
            if chart not in self._atlas:
                raise ValueError("the {} does not belong to ".format(chart) +
                                 "the atlas of {}".format(self))
            chart.restrict(resu, restrictions)
        # Transition maps on the result inferred from those of self:
        for chart1 in coord_def:
            for chart2 in coord_def:
                if chart2 != chart1 and (chart1, chart2) in self._coord_changes:
                    self._coord_changes[(chart1, chart2)].restrict(resu)
        #!# update non-coordinate vector frames and change of frames
        #
        return resu


#*******************************************************************************************

from sage.manifolds.differentiable.tensorfield_paral import TensorFieldParal
from sage.manifolds.differentiable.tensorfield import TensorField


class TangentTensor(TensorFieldParal):
    r"""
    Let ``S`` be a lightlike submanifold embedded in a pseudo-Riemannian
    manifold ``(M,g)`` with ``Phi`` the embedding map. Let ``T1`` be a tensor
    on ``M`` along ``S`` or not. ``TangentTensor(T1, Phi)`` returns the
    restriction ``T2`` of ``T1`` along ``S`` that in addition can be applied
    only on vector fields tangent to ``S``, when ``T1`` has a covariant part.

    INPUT:

    - ``tensor`` -- a tensor field on the ambient manifold
    - ``embedding`` -- the embedding map ``Phi``

    EXAMPLES:

        Section of the lightcone of the Minkowski space with a hyperplane
        passing through the origin::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [sqrt(u^2+v^2), u, v, 0]},
            ....:                  name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x, y]}, name='Phi_inv',
            ....:                      latex_name=r'\Phi^{-1}')
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1, 1, 1, 1
            sage: V = M.vector_field(0,0,0,1)
            sage: S.set_transverse(rigging=t, normal=V)
            sage: xi = M.vector_field(sqrt(x^2+y^2+z^2), x, y, 0)
            sage: U = M.vector_field(0, -y, x, 0)
            sage: Sc = S.screen('Sc', U, xi);
            sage: T1 = M.tensor_field(1,1).along(Phi); T1[0,0] = 1
            sage: V1 = M.vector_field().along(Phi); V1[0] = 1; V1[1]=1
            sage: T1(V1).display()
            ∂/∂t
            sage: from sage.manifolds.differentiable.degenerate_submanifold import TangentTensor
            sage: T2 = TangentTensor(T1, Phi)
            sage: T2
            Tensor field of type (1,1) along the 2-dimensional degenerate
             submanifold S embedded in 4-dimensional differentiable manifold M
             with values on the 4-dimensional Lorentzian manifold M
            sage: V2 = S.projection(V1)
            sage: T2(V2).display()
            u/sqrt(u^2 + v^2) ∂/∂t

        Of course `T1` and `T2` give the same output on vector fields tangent to S::

            sage: T1(xi.along(Phi)).display()
            sqrt(u^2 + v^2) ∂/∂t
            sage: T2(xi.along(Phi)).display()
            sqrt(u^2 + v^2) ∂/∂t
    """
    def __init__(self, tensor, embedding, screen=None):
        r"""

        TESTS::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(3, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v,w> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [u, u, v, w]},
            ....:                  name='Phi', latex_name=r'\Phi');
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x,y, z]}, name='Phi_inv',
            ....:                      latex_name=r'\Phi^{-1}');
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: v = M.vector_field(); v[1] = 1
            sage: S.set_transverse(rigging=v)
            sage: xi = M.vector_field(); xi[0] = 1; xi[1] = 1
            sage: U = M.vector_field(); U[2] = 1; V = M.vector_field(); V[3] = 1
            sage: Sc = S.screen('Sc', (U,V), xi);
            sage: T1 = M.tensor_field(1,1).along(Phi); T1[0,0] = 1
            sage: from sage.manifolds.differentiable.degenerate_submanifold import TangentTensor
            sage: T2 = TangentTensor(T1, Phi); T2
            Tensor field of type (1,1) along the degenerate hypersurface S embedded in
            4-dimensional differentiable manifold M with values on the 4-dimensional
            Lorentzian manifold M
        """
        if not isinstance(tensor, TensorField):
            raise TypeError("the second argument must be a tensor field")
        self._tensor = tensor
        self._tensor.set_name(tensor._name, latex_name=tensor._latex_name)
        self._embedding = embedding
        try:
            tensor = tensor.along(embedding)
        except ValueError:
            pass
        if isinstance(tensor, TensorFieldParal):
            TensorFieldParal.__init__(self, tensor._vmodule, tensor._tensor_type, name=tensor._name,
                   latex_name=tensor._latex_name, sym=tensor._sym, antisym=tensor._antisym)
        else:
            TensorField.__init__(self, tensor._vmodule, tensor._tensor_type, name=tensor._name,
                   latex_name=tensor._latex_name, sym=tensor._sym, antisym=tensor._antisym)
        f = tensor._domain._ambient.default_frame().along(embedding)
        self[f, :] = tensor[f, :]
        frame = self._domain.adapted_frame(screen)
        self.display(frame)
        for i in self._domain._ambient.index_generator(tensor.tensor_rank()):
            for j in range(len(i)):
                if i[j] == self._domain._ambient._dim-self._domain._sindex-1:
                    self[frame, i] = 0

    def __call__(self, *args):
        r"""

        TESTS::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(3, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v,w> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [u, u, v, w]},
            ....:         name='Phi', latex_name=r'\Phi');
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x,y, z]}, name='Phi_inv',
            ....:           latex_name=r'\Phi^{-1}');
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: v = M.vector_field(); v[1] = 1
            sage: S.set_transverse(rigging=v)
            sage: xi = M.vector_field(); xi[0] = 1; xi[1] = 1
            sage: U = M.vector_field(); U[2] = 1; V = M.vector_field(); V[3] = 1
            sage: Sc = S.screen('Sc', (U,V), xi);
            sage: T1 = M.tensor_field(1,1).along(Phi); T1[0,0] = 1
            sage: from sage.manifolds.differentiable.degenerate_submanifold import TangentTensor
            sage: T2 = TangentTensor(T1, Phi); T2(xi.along(Phi))
            Vector field along the degenerate hypersurface S embedded in
            4-dimensional differentiable manifold M with values on the 4-dimensional
            Lorentzian manifold M
        """
        for vector in args:
            try:
                vector = vector.along(self._embedding)
            except ValueError:
                pass
            if not self._domain.is_tangent(vector):
                raise ValueError("The provided vector field is not " +
                        "tangent to {}".format(self._domain._name))
        try:
            return TensorField.__call__(self._tensor.along(self._embedding), *args)
        except ValueError:
            return TensorField.__call__(self._tensor, *args)

    def extension(self):
        r"""

        Return initial tensor

        EXAMPLES:

        Section of the lightcone of the Minkowski space with a hyperplane
        passing through the origin::

            sage: M = Manifold(4, 'M', structure='Lorentzian')
            sage: X.<t,x,y,z> = M.chart()
            sage: S = Manifold(2, 'S', ambient=M, structure='degenerate_metric')
            sage: X_S.<u,v> = S.chart()
            sage: Phi = S.diff_map(M, {(X_S, X): [sqrt(u^2+v^2), u, v, 0]},
            ....:               name='Phi', latex_name=r'\Phi')
            sage: Phi_inv = M.diff_map(S, {(X, X_S): [x, y]}, name='Phi_inv',
            ....:                       latex_name=r'\Phi^{-1}')
            sage: S.set_immersion(Phi, inverse=Phi_inv); S.declare_embedding()
            sage: g = M.metric()
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1,1,1,1
            sage: V = M.vector_field(); V[3] = 1
            sage: S.set_transverse(rigging=t, normal=V)
            sage: xi = M.vector_field(); xi[0] = sqrt(x^2+y^2+z^2); xi[1] = x; xi[2] = y
            sage: U = M.vector_field(); U[1] = -y; U[2] = x
            sage: Sc = S.screen('Sc', U, xi);
            sage: T1 = M.tensor_field(1,1).along(Phi); T1[0,0] = 1
            sage: from sage.manifolds.differentiable.degenerate_submanifold import TangentTensor
            sage: T2 = TangentTensor(T1, Phi); T3 = T2.extension()
            sage: T3 is T2
            False
            sage: T3 is T1
            True
        """
        return self._tensor
