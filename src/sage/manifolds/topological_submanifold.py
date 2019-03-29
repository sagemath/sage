r"""
Submanifolds of topological manifolds

Given a topological manifold `M` over a topological field `K`, a *topological
submanifold of* `M` is defined by a topological manifold `N` over the same
field `K` of dimension lower than the dimension of `M` and a topological
embedding `\phi` from `N` to `M` (i.e. `\phi` is a homeomorphism onto its
image).

In the case where the map `\phi` is only an embedding locally, it is called an
*topological immersion*, and defines an *immersed submanifold*.

The global embedding property cannot be checked in sage, so the immersed or
embedded aspect of the manifold must be declared by the user, by calling either
:meth:`~sage.manifolds.topological_submanifold.TopologicalSubmanifold.set_embedding`
or
:meth:`~sage.manifolds.topological_submanifold.TopologicalSubmanifold.set_immersion`
while declaring the map `\phi`.

The map `\phi: N\to M` can also depend on one or multiple parameters. As long
as `\phi` remains injective in these parameters, it represents a *foliation*.
The *dimension* of the foliation is defined as the number of parameters.

AUTHORS:

- Florentin Jaffredo (2018): initial version

REFERENCES:

- [Lee2013]_

"""


# *****************************************************************************
#  Copyright (C) 2018 Florentin Jaffredo <florentin.jaffredo@polytechnique.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************

from sage.manifolds.manifold import TopologicalManifold
from sage.manifolds.continuous_map import ContinuousMap
from sage.symbolic.expression import Expression
from sage.symbolic.assumptions import assumptions, assume
from sage.plot.plot3d.parametric_surface import ParametricSurface

#############################################################################
# Global options

#############################################################################
# Class


class TopologicalSubmanifold(TopologicalManifold):
    r"""
    Submanifold of a topological manifold.

    Given a topological manifold `M` over a topological field `K`, a
    *topological submanifold of* `M` is defined by a topological manifold `N`
    over the same field `K` of dimension lower than the dimension of `M` and
    a topological embedding `\phi` from `N` to `M` (i.e. `\phi` is an
    homeomorphism onto its image).

    In the case where `\phi` is only an topological immersion (i.e. is only
    locally an embedding), one says that `N` is an *immersed submanifold*.

    The map `\phi` can also depend on one or multiple parameters.
    As long as `\phi` remains injective in these parameters, it represents
    a *foliation*. The *dimension* of the foliation is defined as the number of
    parameters.

    INPUT:

    - ``n`` -- positive integer; dimension of the manifold
    - ``name`` -- string; name (symbol) given to the manifold
    - ``field`` -- field `K` on which the manifold is defined; allowed values
      are

      - ``'real'`` or an object of type ``RealField`` (e.g., ``RR``) for
        a manifold over `\RR`
      - ``'complex'`` or an object of type ``ComplexField`` (e.g., ``CC``)
        for a manifold over `\CC`
      - an object in the category of topological fields (see
        :class:`~sage.categories.fields.Fields` and
        :class:`~sage.categories.topological_spaces.TopologicalSpaces`)
        for other types of manifolds

    - ``structure`` -- manifold structure (see
      :class:`~sage.manifolds.structure.TopologicalStructure` or
      :class:`~sage.manifolds.structure.RealTopologicalStructure`)
    - ``ambient`` -- (default: ``None``) manifold of destination
      of the immersion. If ``None``, set to ``self``
    - ``base_manifold`` -- (default: ``None``) if not ``None``, must be a
      topological manifold; the created object is then an open subset of
      ``base_manifold``
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the manifold; if none are provided, it is set to ``name``
    - ``start_index`` -- (default: 0) integer; lower value of the range of
      indices used for "indexed objects" on the manifold, e.g., coordinates
      in a chart
    - ``category`` -- (default: ``None``) to specify the category; if
      ``None``, ``Manifolds(field)`` is assumed (see the category
      :class:`~sage.categories.manifolds.Manifolds`)
    - ``unique_tag`` -- (default: ``None``) tag used to force the construction
      of a new object when all the other arguments have been used previously
      (without ``unique_tag``, the
      :class:`~sage.structure.unique_representation.UniqueRepresentation`
      behavior inherited from
      :class:`~sage.manifolds.subset.ManifoldSubset`
      would return the previously constructed object corresponding to these
      arguments)

    EXAMPLES:

    Let `N` be a 2-dimensional submanifold of a 3-dimensional manifold `M`::

        sage: M = Manifold(3, 'M', structure="topological")
        sage: N = Manifold(2, 'N', ambient=M, structure="topological")
        sage: N
        2-dimensional submanifold N embedded in 3-dimensional manifold M
        sage: CM.<x,y,z> = M.chart()
        sage: CN.<u,v> = N.chart()

    Let us define a 1-dimensional foliation indexed by `t`. The inverse map is
    needed in order to compute the adapted chart in the ambient manifold::

        sage: t = var('t')
        sage: phi = N.continuous_map(M, {(CN,CM):[u, v, t+u**2+v**2]}); phi
        Continuous map from the 2-dimensional submanifold N embedded in
         3-dimensional manifold M to the 3-dimensional topological manifold M
        sage: phi_inv = M.continuous_map(N, {(CM, CN):[x, y]})
        sage: phi_inv_t = M.scalar_field({CM: z-x**2-y**2})

    `\phi` can then be declared as an embedding `N\to M`::

        sage: N.set_embedding(phi, inverse=phi_inv, var=t,
        ....:                 t_inverse={t: phi_inv_t})

    The foliation can also be used to find new charts on the ambient manifold
    that are adapted to the foliation, i.e. in which the expression of the
    immersion is trivial. At the same time, the appropriate coordinate changes
    are computed::

        sage: N.adapted_chart()
        [Chart (M, (u_M, v_M, t_M))]
        sage: len(M.coord_changes())
        2

    The foliations parameters are always added as the last coordinates.

    .. SEEALSO::

        :mod:`~sage.manifolds.manifold`

    """
    def __init__(self, n, name, field, structure, ambient=None,
                 base_manifold=None, latex_name=None, start_index=0,
                 category=None, unique_tag=None):
        r"""
        Construct a submanifold of a topological manifold.

        EXAMPLES::
            sage: M = Manifold(3, 'M', structure="topological")
            sage: N = Manifold(2, 'N', ambient=M, structure="topological")
            sage: N
            2-dimensional submanifold N embedded in 3-dimensional manifold M

        """
        TopologicalManifold.__init__(self, n, name, field, structure,
                                     base_manifold=base_manifold,
                                     latex_name=latex_name,
                                     start_index=start_index,
                                     category=category)
        self._immersion = None
        self._immersion_inv = None
        self._var = None
        self._dim_foliation = 0
        self._t_inverse = {}
        if ambient is None:
            self._ambient = self
        else:
            self._ambient = ambient
        self._immersed = False
        self._embedded = False
        self._adapted_charts = None
        self._subs = None

    def _repr_(self):
        r"""
        Return a string representation of the submanifold.

        If no ambient manifold is specified, the submanifold is considered as
        a manifold.

        TESTS::

            sage: M = Manifold(3, 'M', structure="topological")
            sage: N = Manifold(2, 'N', ambient=M, structure="topological")
            sage: N._repr_()
            '2-dimensional submanifold N embedded in 3-dimensional manifold M'

        """
        if self._ambient == self:
            return super(TopologicalManifold, self).__repr__()
        return "{}-dimensional submanifold {} embedded in {}-" \
               "dimensional manifold {}".format(self._dim, self._name,
                                                self._ambient._dim,
                                                self._ambient._name)

    def set_immersion(self, phi, inverse=None, var=None,
                      t_inverse=None):
        r"""
        Register the immersion of the immersed submanifold.

        A *topological immersion* is a continuous map that is locally a
        topological embedding (i.e. a homeomorphism onto its image).
        A *differentiable immersion* is a differentiable map whose differential
        is injective at each point.

        If an inverse of the immersion onto its image exists, it can be
        registered at the same time. If the immersion depends on parameters,
        they must also be declared here.

        INPUTS:

        - ``phi`` -- continuous map `\phi` from self to self._ambient
        - ``inverse`` -- (default: ``None``) inverse of `\phi` onto its image,
          used for computing changes of chart from or to adapted charts. No
          verification is made
        - ``var`` -- (default: ``None``) list of parameters appearing in `\phi`
        - ``t_inverse`` -- (default: ``None``) dictionary of scalar field on
          self._ambient indexed by elements of ``var`` representing the missing
          information in ``inverse``

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure="topological")
            sage: N = Manifold(2, 'N', ambient=M, structure="topological")
            sage: N
            2-dimensional submanifold N embedded in 3-dimensional manifold M
            sage: CM.<x,y,z> = M.chart()
            sage: CN.<u,v> = N.chart()
            sage: t = var('t')
            sage: phi = N.continuous_map(M,{(CN,CM):[u,v,t+u**2+v**2]}); phi
            Continuous map from the 2-dimensional submanifold N embedded in
             3-dimensional manifold M to the 3-dimensional topological
             manifold M
            sage: phi_inv = M.continuous_map(N,{(CM,CN):[x,y]})
            sage: phi_inv_t = M.scalar_field({CM:z-x**2-y**2})
            sage: N.set_immersion(phi, inverse=phi_inv, var=t,
            ....:                 t_inverse={t: phi_inv_t})

        """
        if not isinstance(phi, ContinuousMap):
            raise TypeError("phi must be a or differentiable (or at least"
                            " continuous) map")
        if phi._domain is not self or phi._codomain is not self._ambient:
            raise ValueError("{} is not a homeomorphism "
                             "from {} to {}".format(phi._name, self._name,
                                                    self._ambient.name()))
        self._immersion = phi

        if inverse is not None:
            self._immersion._inverse = inverse
            self._immersion_inv = inverse

        if var is not None:
            try:
                iter(var)
                for v in var:
                    if not isinstance(v, Expression):
                        raise TypeError()
            except TypeError:
                if not isinstance(var, Expression):
                    raise TypeError("var must be a variable "
                                    "or list of variables")

            if isinstance(var, Expression):
                self._var = [var]
                self._dim_foliation = 1
            else:
                self._var = var
                self._dim_foliation = len(var)
        if t_inverse is None:
            t_inverse = {}

        self._t_inverse = t_inverse
        self._immersed = True

    def declare_embedding(self):
        r"""
        Declare that the immersion provided by :meth:`set_immersion` is in
        fact an embedding.

        A *topological embedding* is a continuous map that is a homeomorphism
        onto its image. A *differentiable embedding* is a topological embedding
        that is also a differentiable immersion.

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure="topological")
            sage: N = Manifold(2, 'N', ambient=M, structure="topological")
            sage: N
            2-dimensional submanifold N embedded in 3-dimensional manifold M
            sage: CM.<x,y,z> = M.chart()
            sage: CN.<u,v> = N.chart()
            sage: t = var('t')
            sage: phi = N.continuous_map(M,{(CN,CM):[u,v,t+u**2+v**2]}); phi
            Continuous map from the 2-dimensional submanifold N embedded in
             3-dimensional manifold M to the 3-dimensional topological
             manifold M
            sage: phi_inv = M.continuous_map(N,{(CM,CN):[x,y]})
            sage: phi_inv_t = M.scalar_field({CM:z-x**2-y**2})
            sage: N.set_immersion(phi, inverse=phi_inv, var=t,
            ....:                 t_inverse={t: phi_inv_t})
            sage: N._immersed
            True
            sage: N._embedded
            False
            sage: N.declare_embedding()
            sage: N._immersed
            True
            sage: N._embedded
            True

        """
        if not self._immersed:
            raise ValueError("please declare an embedding using set_immersion "
                             "before calling declare_embedding()")
        self._immersion._is_isomorphism = True
        self._embedded = True

    def set_embedding(self, phi, inverse=None, var=None,
                      t_inverse=None):
        r"""
        Register the embedding of an embedded submanifold.

        A *topological embedding* is a continuous map that is a homeomorphism
        onto its image. A *differentiable embedding* is a topological embedding
        that is also a differentiable immersion.

        INPUTS:

        - ``phi`` -- continuous map `\phi` from self to self._ambient
        - ``inverse`` -- (default: ``None``) inverse of `\phi` onto its image,
          used for computing changes of chart from or to adapted charts. No
          verification is made
        - ``var`` -- (default: ``None``) list of parameters appearing in `\phi`
        - ``t_inverse`` -- (default: ``None``) dictionary of scalar field on
          self._ambient indexed by elements of ``var`` representing the missing
          information in ``inverse``

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure="topological")
            sage: N = Manifold(2, 'N', ambient=M, structure="topological")
            sage: N
            2-dimensional submanifold N embedded in 3-dimensional manifold M
            sage: CM.<x,y,z> = M.chart()
            sage: CN.<u,v> = N.chart()
            sage: t = var('t')
            sage: phi = N.continuous_map(M,{(CN,CM):[u,v,t+u**2+v**2]}); phi
            Continuous map from the 2-dimensional submanifold N embedded in
             3-dimensional manifold M to the 3-dimensional topological
             manifold M
            sage: phi_inv = M.continuous_map(N,{(CM,CN):[x,y]})
            sage: phi_inv_t = M.scalar_field({CM:z-x**2-y**2})
            sage: N.set_embedding(phi, inverse=phi_inv, var=t,
            ....:                 t_inverse={t: phi_inv_t})

        """
        self.set_immersion(phi, inverse, var, t_inverse)
        self.declare_embedding()


    def adapted_chart(self, index="", latex_index=""):
        r"""
        Create charts and changes of charts in the ambient manifold adapted
        to the foliation.

        A manifold `M` of dimension `m` can be foliated by submanifolds `N` of
        dimension `n`. The corresponding embedding needs `m-n` free parameters
        to describe the whole manifold.

        A set of coordinates adapted to a foliation is a set of coordinates
        `(x_1,...,x_n,t_1,...t_{m-n})` such that `(x_1,...x_n)` are coordinates
        of `N` and `(t_1,...t_{m-n})` are the `m-n` free parameters of the
        foliation.

        Provided that an embedding with free variables is already defined, this
        function constructs such charts and coordinates changes whenever
        it is possible.

        If there are restrictions of the coordinates on the starting chart,
        these restrictions are also propagated.

        INPUT:

        - ``index`` -- (default: ``""``) string defining the name of the
          coordinates in the new chart. This string will be added at the end of
          the names of the old coordinates. By default, it is replaced by
          ``"_"+self._ambient._name``
        - ``latex_index`` -- (default: ``""``) string defining the latex name
          of the coordinates in the new chart. This string will be added at the
          end of the latex names of the old coordinates. By default, it is
          replaced by ``"_"+self._ambient._latex_()``

        OUTPUT:

        - list of charts created from the charts of ``self``

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure="topological")
            sage: N = Manifold(2, 'N', ambient=M, structure="topological")
            sage: N
            2-dimensional submanifold N embedded in 3-dimensional manifold M
            sage: CM.<x,y,z> = M.chart()
            sage: CN.<u,v> = N.chart()
            sage: t = var('t')
            sage: phi = N.continuous_map(M,{(CN,CM):[u,v,t+u**2+v**2]})
            sage: phi_inv = M.continuous_map(N,{(CM,CN):[x,y]})
            sage: phi_inv_t = M.scalar_field({CM:z-x**2-y**2})
            sage: N.set_immersion(phi, inverse=phi_inv, var=t,
            ....:                 t_inverse={t:phi_inv_t})
            sage: N.declare_embedding()
            sage: N.adapted_chart()
            [Chart (M, (u_M, v_M, t_M))]
        """
        if not self._embedded:
            raise ValueError("an embedding is required")

        if self._dim_foliation+self._dim != self._ambient._dim:
            raise ValueError("a foliation of dimension Dim(M)-Dim(N) is "
                             "needed to find an adapted chart")
        if not isinstance(index, str):
            raise TypeError("index must be a string")

        res = []
        self._subs = []

        # All possible expressions for the immersion
        domains = self._immersion._coord_expression.keys()
        postscript = index
        if index == "":
            postscript = "_" + self._ambient._name

        latex_postscript = latex_index
        if latex_postscript == "":
            latex_postscript = "_" + self._ambient._latex_()

        for domain in list(domains):
            name = " ".join([domain[0][i]._repr_()+postscript+":{"
                             + domain[0][i]._latex_()+"}"+latex_postscript
                             for i in self.irange()]) + " "\
                   + " ".join(v._repr_()+postscript for v in self._var)
            chart = domain[1]._domain.chart(name)
            if chart not in res:

                # Construct restrictions on coordinates:
                subs = {domain[0][i]: chart[:][i] for i in range(self._dim)}
                for i in range(len(self._var)):
                    subs[self._var[i]] = chart[:][self._dim+i]
                for rest in domain[0]._restrictions:
                    chart.add_restrictions(rest.subs(subs))
                for _a in assumptions(*(domain[0][:]+tuple(self._var))):
                    if isinstance(_a, Expression):
                        assume(_a.subs(subs))

                self._subs.append(subs)
                res.append(chart)
                self._immersion.add_expr(domain[0], chart,
                                         list(domain[0][:]) + self._var)
                self._immersion_inv.add_expr(chart, domain[0],
                                             chart[:][0:self._dim])
                for i in range(len(self._var)):
                    self._t_inverse[self._var[i]].add_expr(
                        chart[:][self._dim:][i], chart=chart)

        for (chartNV, chartMV) in self._immersion._coord_expression:
            for (chartNU, chartMU) in self._immersion._coord_expression:
                if chartMU is not chartMV and\
                        (chartMU, chartMV) not in self._ambient._coord_changes:
                    if (chartNU, chartNV) in self._coord_changes or \
                            chartNU is chartNV:
                        _f = self._immersion.coord_functions(chartNV, chartMV)
                        _g = self._coord_changes[(chartNU, chartNV)]._transf \
                            if chartNU is not chartNV else lambda *x: x
                        _h = self._immersion_inv.coord_functions(chartMU,
                                                                 chartNU)
                        expr = list(_f(*_g(*_h(*chartMU[:]))))
                        substitutions = {v: self._t_inverse[v].expr(chartMU)
                                         for v in self._var}
                        for i in range(len(expr)):
                            expr[i] = expr[i].subs(substitutions)

                        chartMU.transition_map(chartMV, expr)
        self._adapted_charts = res
        return res

    def plot(self, param, u, v, chart1=None, chart2=None, **kwargs):
        r"""
        Plot an embedding.

        Plot the embedding defined by the foliation and a set of values for the
        free parameters. This function can only plot 2-dimensional surfaces
        embedded in 3-dimensional manifolds. It ultimately calls
        :class:`~sage.plot.plot3d.parametric_surface.ParametricSurface`.

        INPUT:

        - ``param`` -- dictionary of values indexed by the free variables
          appearing in the foliation.
        - ``u`` -- iterable of the values taken by the first coordinate of the
          surface to plot
        - ``v`` -- iterable of the values taken by the second coordinate of the
          surface to plot
        - ``chart1`` -- (default: ``None``) chart in which ``u`` and ``v`` are
          considered. By default, the default chart of the submanifold is used
        - ``chart1`` -- (default: ``None``) destination chart. By default, the
          default chart of the manifold is used
        - ``**kwargs`` -- other arguments as used in
          :class:`~sage.plot.plot3d.parametric_surface.ParametricSurface`

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure="topological")
            sage: N = Manifold(2, 'N', ambient = M, structure="topological")
            sage: CM.<x,y,z> = M.chart()
            sage: CN.<u,v> = N.chart()
            sage: t = var('t')
            sage: phi = N.continuous_map(M,{(CN,CM):[u,v,t+u**2+v**2]})
            sage: phi_inv = M.continuous_map(N,{(CM,CN):[x,y]})
            sage: phi_inv_t = M.scalar_field({CM:z-x**2-y**2})
            sage: N.set_immersion(phi, inverse=phi_inv, var=t,
            ....:                 t_inverse = {t:phi_inv_t})
            sage: N.declare_embedding()
            sage: N.adapted_chart()
            [Chart (M, (u_M, v_M, t_M))]
            sage: P0 = N.plot({t:0}, srange(-1, 1, 0.1), srange(-1, 1, 0.1),
            ....:             CN, CM, opacity=0.3, mesh=True)
            sage: P1 = N.plot({t:1}, srange(-1, 1, 0.1), srange(-1, 1, 0.1),
            ....:             CN, CM, opacity=0.3, mesh=True)
            sage: P2 = N.plot({t:2}, srange(-1, 1, 0.1), srange(-1, 1, 0.1),
            ....:             CN, CM, opacity=0.3, mesh=True)
            sage: P3 = N.plot({t:3}, srange(-1, 1, 0.1), srange(-1, 1, 0.1),
            ....:             CN, CM, opacity=0.3, mesh=True)
            sage: show(P0+P1+P2+P3)

        .. SEEALSO::

            :class:`~sage.plot.plot3d.parametric_surface.ParametricSurface`

        """

        if self._dim != 2 or self._ambient._dim != 3:
            raise ValueError("plot only for 2-dimensional hypersurfaces")
        if chart1 is None:
            chart1 = self.default_chart()
        if chart2 is None:
            chart2 = self._ambient.default_chart()
        expr = list(self._immersion.coord_functions(chart1, chart2))
        for i in range(len(expr)):
            expr[i] = expr[i].expr().subs(param)
        fx = expr[0].function(*chart1[:])
        fy = expr[1].function(*chart1[:])
        fz = expr[2].function(*chart1[:])

        return ParametricSurface((fx, fy, fz), (u, v), **kwargs)

    def ambient(self):
        r"""
        Return the ambient manifold in which ``self`` is immersed or embedded.

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure="topological")
            sage: N = Manifold(2, 'N', ambient=M, structure="topological")
            sage: N.ambient()
            3-dimensional topological manifold M
        """
        return self._ambient

    def immersion(self):
        r"""
        Return the immersion of the submanifold.

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure="topological")
            sage: N = Manifold(2, 'N', ambient=M, structure="topological")
            sage: CM.<x,y,z> = M.chart()
            sage: CN.<u,v> = N.chart()
            sage: t = var('t')
            sage: phi = N.continuous_map(M,{(CN,CM):[u,v,t+u**2+v**2]})
            sage: phi_inv = M.continuous_map(N,{(CM,CN):[x,y]})
            sage: phi_inv_t = M.scalar_field({CM:z-x**2-y**2})
            sage: N.set_immersion(phi, inverse=phi_inv, var=t,
            ....:                 t_inverse={t: phi_inv_t})
            sage: N.immersion()
            Continuous map from the 2-dimensional submanifold N embedded in
             3-dimensional manifold M to the 3-dimensional topological
             manifold M
        """
        if not self._immersed:
            raise ValueError("the submanifold is not immersed")
        return self._immersion

    def embedding(self):
        r"""
        Return the embedding of the submanifold.

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure="topological")
            sage: N = Manifold(2, 'N', ambient=M, structure="topological")
            sage: CM.<x,y,z> = M.chart()
            sage: CN.<u,v> = N.chart()
            sage: t = var('t')
            sage: phi = N.continuous_map(M,{(CN,CM):[u,v,t+u**2+v**2]})
            sage: phi_inv = M.continuous_map(N,{(CM,CN):[x,y]})
            sage: phi_inv_t = M.scalar_field({CM:z-x**2-y**2})
            sage: N.set_embedding(phi, inverse=phi_inv, var=t,
            ....:                 t_inverse={t: phi_inv_t})
            sage: N.embedding()
            Homeomorphism from the 2-dimensional submanifold N embedded in
             3-dimensional manifold M to the 3-dimensional topological manifold
             M
        """
        if not self._embedded:
            raise ValueError("the submanifold is not embedded")
        return self._immersion
