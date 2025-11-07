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
- Eric Gourgoulhon (2018-2019): add documentation
- Matthias Koeppe (2021): open subsets of submanifolds

REFERENCES:

- \J. M. Lee:  *Introduction to Smooth Manifolds* [Lee2013]_
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

from sage.manifolds.continuous_map import ContinuousMap
from sage.manifolds.manifold import TopologicalManifold
from sage.misc.lazy_import import lazy_import
from sage.symbolic.assumptions import assume, assumptions
from sage.symbolic.expression import Expression

lazy_import("sage.plot.plot3d.parametric_surface", "ParametricSurface")

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

    In the case where `\phi` is only a topological immersion (i.e. is only
    locally an embedding), one says that `N` is an *immersed submanifold*.

    The map `\phi` can also depend on one or multiple parameters.
    As long as `\phi` remains injective in these parameters, it represents
    a *foliation*. The *dimension* of the foliation is defined as the number of
    parameters.

    INPUT:

    - ``n`` -- positive integer; dimension of the submanifold
    - ``name`` -- string; name (symbol) given to the submanifold
    - ``field`` -- field `K` on which the submanifold is defined; allowed
      values are

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
    - ``ambient`` -- (default: ``None``) codomain `M` of the immersion `\phi`;
      must be a topological manifold. If ``None``, it is set to ``self``
    - ``base_manifold`` -- (default: ``None``) if not ``None``, must be a
      topological manifold; the created object is then an open subset of
      ``base_manifold``
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the submanifold; if none are provided, it is set to ``name``
    - ``start_index`` -- (default: 0) integer; lower value of the range of
      indices used for "indexed objects" on the submanifold, e.g., coordinates
      in a chart
    - ``category`` -- (default: ``None``) to specify the category; if
      ``None``, ``Manifolds(field)`` is assumed (see the category
      :class:`~sage.categories.manifolds.Manifolds`)
    - ``unique_tag`` -- (default: ``None``) tag used to force the construction
      of a new object when all the other arguments have been used previously
      (without ``unique_tag``, the
      :class:`~sage.structure.unique_representation.UniqueRepresentation`
      behavior inherited from
      :class:`~sage.manifolds.subset.ManifoldSubset` via
      :class:`~sage.manifolds.manifold.TopologicalManifold`
      would return the previously constructed object corresponding to these
      arguments)

    EXAMPLES:

    Let `N` be a 2-dimensional submanifold of a 3-dimensional manifold `M`::

        sage: M = Manifold(3, 'M', structure='topological')
        sage: N = Manifold(2, 'N', ambient=M, structure='topological')
        sage: N
        2-dimensional topological submanifold N immersed in the 3-dimensional
         topological manifold M
        sage: CM.<x,y,z> = M.chart()
        sage: CN.<u,v> = N.chart()

    Let us define a 1-dimensional foliation indexed by `t`::

        sage: t = var('t')
        sage: phi = N.continuous_map(M, {(CN,CM): [u, v, t+u^2+v^2]})
        sage: phi.display()
        N → M
           (u, v) ↦ (x, y, z) = (u, v, u^2 + v^2 + t)

    The foliation inverse maps are needed for computing the adapted chart on
    the ambient manifold::

        sage: phi_inv = M.continuous_map(N, {(CM, CN): [x, y]})
        sage: phi_inv.display()
        M → N
           (x, y, z) ↦ (u, v) = (x, y)
        sage: phi_inv_t = M.scalar_field({CM: z-x^2-y^2})
        sage: phi_inv_t.display()
        M → ℝ
        (x, y, z) ↦ -x^2 - y^2 + z

    `\phi` can then be declared as an embedding `N\to M`::

        sage: N.set_embedding(phi, inverse=phi_inv, var=t,
        ....:                 t_inverse={t: phi_inv_t})

    The foliation can also be used to find new charts on the ambient manifold
    that are adapted to the foliation, i.e. in which the expression of the
    immersion is trivial. At the same time, the appropriate coordinate changes
    are computed::

        sage: N.adapted_chart()
        [Chart (M, (u_M, v_M, t_M))]
        sage: M.atlas()
        [Chart (M, (x, y, z)), Chart (M, (u_M, v_M, t_M))]
        sage: len(M.coord_changes())
        2

    The foliation parameters are always added as the last coordinates.

    .. SEEALSO::

        :mod:`~sage.manifolds.manifold`
    """
    def __init__(self, n, name, field, structure, ambient=None,
                 base_manifold=None, latex_name=None, start_index=0,
                 category=None, unique_tag=None):
        r"""
        Construct a submanifold of a topological manifold.

        TESTS::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: N = Manifold(2, 'N', ambient=M, structure='topological')
            sage: N
            2-dimensional topological submanifold N immersed in the
             3-dimensional topological manifold M
        """
        TopologicalManifold.__init__(self, n, name, field, structure,
                                     base_manifold=base_manifold,
                                     latex_name=latex_name,
                                     start_index=start_index,
                                     category=category)
        if not (ambient is None
                or isinstance(ambient, TopologicalManifold)):
            raise TypeError("ambient must be a manifold")
        self._init_immersion(ambient=ambient)

    def _init_immersion(self, ambient=None):
        r"""
        Initialize the attributes relative to the immersion of ``self`` in
        the ambient manifold.

        INPUT:

        - ``ambient`` -- (default: ``None``) codomain of the immersion;
          must be a topological manifold. If ``None``, it is set to ``self``

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: N = Manifold(1, 'N', ambient=M, structure='topological')
            sage: N._init_immersion(ambient=M)
        """
        self._immersion = None
        self._immersion_inv = None
        self._var = None
        self._dim_foliation = 0
        self._t_inverse = {}
        if ambient is None:
            self._ambient = self
        else:
            self._ambient = ambient
            self._codim = ambient._dim - self._dim
            if self._codim < 0:
                raise ValueError("the submanifold must be of smaller "
                                 + "dimension than its ambient manifold")
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

            sage: M = Manifold(3, 'M', structure='topological')
            sage: N = Manifold(2, 'N', ambient=M, structure='topological')
            sage: N
            2-dimensional topological submanifold N immersed in the
             3-dimensional topological manifold M
            sage: phi = N.continuous_map(M)
            sage: N.set_embedding(phi)
            sage: N
            2-dimensional topological submanifold N embedded in the
             3-dimensional topological manifold M
        """
        if self is not self._manifold:
            return "Open subset {} of the {}".format(self._name, self._manifold)
        if self._ambient is self:
            return super(TopologicalManifold, self).__repr__()
        if self._embedded:
            return "{}-dimensional {} submanifold {} embedded in the {}".format(
                self._dim, self._structure.name, self._name, self._ambient)
        return "{}-dimensional {} submanifold {} immersed in the {}".format(
                self._dim, self._structure.name, self._name, self._ambient)

    def open_subset(self, name, latex_name=None, coord_def={}, supersets=None):
        r"""
        Create an open subset of the manifold.

        An open subset is a set that is (i) included in the manifold and (ii)
        open with respect to the manifold's topology. It is a topological
        manifold by itself.

        As ``self`` is a submanifold of its ambient manifold,
        the new open subset is also considered a submanifold of that.
        Hence the returned object is an instance of
        :class:`TopologicalSubmanifold`.

        INPUT:

        - ``name`` -- name given to the open subset
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote
          the subset; if none are provided, it is set to ``name``
        - ``coord_def`` -- (default: {}) definition of the subset in
          terms of coordinates; ``coord_def`` must a be dictionary with keys
          charts on the manifold and values the symbolic expressions formed
          by the coordinates to define the subset
        - ``supersets`` -- (default: only ``self``) list of sets that the
          new open subset is a subset of

        OUTPUT: the open subset, as an instance of :class:`TopologicalSubmanifold`

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: N = Manifold(2, 'N', ambient=M, structure='topological'); N
            2-dimensional topological submanifold N immersed in the
             3-dimensional topological manifold M
            sage: S = N.subset('S'); S
            Subset S of the
             2-dimensional topological submanifold N immersed in the
              3-dimensional topological manifold M
            sage: O = N.subset('O', is_open=True); O  # indirect doctest
            Open subset O of the
             2-dimensional topological submanifold N immersed in the
              3-dimensional topological manifold M

            sage: phi = N.continuous_map(M)
            sage: N.set_embedding(phi)
            sage: N
            2-dimensional topological submanifold N embedded in the
             3-dimensional topological manifold M
            sage: S = N.subset('S'); S
            Subset S of the
             2-dimensional topological submanifold N embedded in the
              3-dimensional topological manifold M
            sage: O = N.subset('O', is_open=True); O  # indirect doctest
            Open subset O of the
             2-dimensional topological submanifold N embedded in the
              3-dimensional topological manifold M
        """
        resu = TopologicalSubmanifold(self._dim, name, self._field,
                                      self._structure, self._ambient,
                                      base_manifold=self._manifold,
                                      latex_name=latex_name,
                                      start_index=self._sindex)
        if supersets is None:
            supersets = [self]
        for superset in supersets:
            superset._init_open_subset(resu, coord_def=coord_def)
        return resu

    def _init_open_subset(self, resu, coord_def):
        r"""
        Initialize ``resu`` as an open subset of ``self``.

        INPUT:

        - ``resu`` -- an instance of :class:`TopologicalManifold` or
          a subclass

        - ``coord_def`` -- (default: ``{}``) definition of the subset in
          terms of coordinates; ``coord_def`` must a be dictionary with keys
          charts on the manifold and values the symbolic expressions formed
          by the coordinates to define the subset

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: N = Manifold(2, 'N', ambient=M, structure='topological')
            sage: phi = N.continuous_map(M)
            sage: N.set_embedding(phi)
            sage: N
            2-dimensional topological submanifold N embedded in the
             3-dimensional topological manifold M
            sage: from sage.manifolds.topological_submanifold import TopologicalSubmanifold
            sage: O = TopologicalSubmanifold(3, 'O', field=M._field, structure=M._structure,
            ....:                            ambient=M, base_manifold=N)
            sage: N._init_open_subset(O, {})
            sage: O
            Open subset O of the
             2-dimensional topological submanifold N embedded in the
              3-dimensional topological manifold M
            sage: O.embedding()
            Continuous map
             from the Open subset O of the 2-dimensional topological submanifold N
              embedded in the 3-dimensional topological manifold M
             to the 3-dimensional topological manifold M
        """
        super()._init_open_subset(resu, coord_def=coord_def)
        ## Extras for Submanifold
        if self._immersed:
            resu.set_immersion(self._immersion.restrict(resu),
                               var=self._var, t_inverse=self._t_inverse)
        if self._embedded:
            resu.declare_embedding()

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

        INPUT:

        - ``phi`` -- continuous map `\phi` from ``self`` to ``self.ambient()``
        - ``inverse`` -- (default: ``None``) continuous map from
          ``self.ambient()`` to ``self``, which once restricted to the image
          of `\phi` is the inverse of `\phi` onto its image if the latter
          exists (NB: no check of this is performed)
        - ``var`` -- (default: ``None``) list of parameters involved in the
          definition of `\phi` (case of foliation); if `\phi` depends on a
          single parameter ``t``, one can write ``var=t`` as a shortcut for
          ``var=[t]``
        - ``t_inverse`` -- (default: ``None``) dictionary of scalar fields on
          ``self.ambient()`` providing the values of the parameters involved
          in the definition of `\phi` (case of foliation), the keys being
          the parameters

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: N = Manifold(2, 'N', ambient=M, structure='topological')
            sage: N
            2-dimensional topological submanifold N immersed in the
             3-dimensional topological manifold M
            sage: CM.<x,y,z> = M.chart()
            sage: CN.<u,v> = N.chart()
            sage: t = var('t')
            sage: phi = N.continuous_map(M, {(CN,CM): [u,v,t+u^2+v^2]})
            sage: phi.display()
            N → M
               (u, v) ↦ (x, y, z) = (u, v, u^2 + v^2 + t)
            sage: phi_inv = M.continuous_map(N, {(CM,CN): [x,y]})
            sage: phi_inv.display()
            M → N
                (x, y, z) ↦ (u, v) = (x, y)
            sage: phi_inv_t = M.scalar_field({CM: z-x^2-y^2})
            sage: phi_inv_t.display()
            M → ℝ
            (x, y, z) ↦ -x^2 - y^2 + z
            sage: N.set_immersion(phi, inverse=phi_inv, var=t,
            ....:                 t_inverse={t: phi_inv_t})
        """
        if not isinstance(phi, ContinuousMap):
            raise TypeError("the argument phi must be a continuous map")
        if phi.domain() is not self or phi.codomain() is not self._ambient:
            raise ValueError("{} is not a map from {} to {}".format(phi, self,
                                                                self._ambient))
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

            sage: M = Manifold(3, 'M', structure='topological')
            sage: N = Manifold(2, 'N', ambient=M, structure='topological')
            sage: N
            2-dimensional topological submanifold N immersed in the
             3-dimensional topological manifold M
            sage: CM.<x,y,z> = M.chart()
            sage: CN.<u,v> = N.chart()
            sage: t = var('t')
            sage: phi = N.continuous_map(M, {(CN,CM): [u,v,t+u^2+v^2]})
            sage: phi_inv = M.continuous_map(N, {(CM,CN): [x,y]})
            sage: phi_inv_t = M.scalar_field({CM: z-x^2-y^2})
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
        self._embedded = True

    def set_embedding(
        self, phi: ContinuousMap, inverse=None, var=None, t_inverse=None
    ):
        r"""
        Register the embedding of an embedded submanifold.

        A *topological embedding* is a continuous map that is a homeomorphism
        onto its image. A *differentiable embedding* is a topological embedding
        that is also a differentiable immersion.

        INPUT:

        - ``phi`` -- continuous map `\phi` from ``self`` to ``self.ambient()``
        - ``inverse`` -- (default: ``None``) continuous map from
          ``self.ambient()`` to ``self``, which once restricted to the image
          of `\phi` is the inverse of `\phi` onto its image (NB: no check of
          this is performed)
        - ``var`` -- (default: ``None``) list of parameters involved in the
          definition of `\phi` (case of foliation); if `\phi` depends on a
          single parameter ``t``, one can write ``var=t`` as a shortcut for
          ``var=[t]``
        - ``t_inverse`` -- (default: ``None``) dictionary of scalar fields on
          ``self.ambient()`` providing the values of the parameters involved
          in the definition of `\phi` (case of foliation), the keys being
          the parameters

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: N = Manifold(2, 'N', ambient=M, structure='topological')
            sage: N
            2-dimensional topological submanifold N immersed in the
             3-dimensional topological manifold M
            sage: CM.<x,y,z> = M.chart()
            sage: CN.<u,v> = N.chart()
            sage: t = var('t')
            sage: phi = N.continuous_map(M, {(CN,CM): [u,v,t+u^2+v^2]})
            sage: phi.display()
            N → M
               (u, v) ↦ (x, y, z) = (u, v, u^2 + v^2 + t)
            sage: phi_inv = M.continuous_map(N, {(CM,CN): [x,y]})
            sage: phi_inv.display()
            M → N
                (x, y, z) ↦ (u, v) = (x, y)
            sage: phi_inv_t = M.scalar_field({CM: z-x^2-y^2})
            sage: phi_inv_t.display()
            M → ℝ
            (x, y, z) ↦ -x^2 - y^2 + z
            sage: N.set_embedding(phi, inverse=phi_inv, var=t,
            ....:                 t_inverse={t: phi_inv_t})

        Now ``N`` appears as an embedded submanifold::

            sage: N
            2-dimensional topological submanifold N embedded in the
             3-dimensional topological manifold M
        """
        self.set_immersion(phi, inverse, var, t_inverse)
        self.declare_embedding()

    def adapted_chart(self, postscript=None, latex_postscript=None):
        r"""
        Create charts and changes of charts in the ambient manifold adapted
        to the foliation.

        A manifold `M` of dimension `m` can be foliated by submanifolds `N` of
        dimension `n`. The corresponding embedding needs `m-n` free parameters
        to describe the whole manifold.

        A chart adapted to the foliation is a set of coordinates
        `(x_1,\ldots,x_n,t_1,\ldots,t_{m-n})` on `M` such that
        `(x_1,\ldots,x_n)` are coordinates on `N` and `(t_1,\ldots,t_{m-n})`
        are the `m-n` free parameters of the foliation.

        Provided that an embedding with free variables is already defined, this
        function constructs such charts and coordinates changes whenever
        it is possible.

        If there are restrictions of the coordinates on the starting chart,
        these restrictions are also propagated.

        INPUT:

        - ``postscript`` -- (default: ``None``) string defining the name of the
          coordinates of the adapted chart. This string will be appended to
          the names of the coordinates `(x_1,\ldots,x_n)` and of the parameters
          `(t_1,\ldots,t_{m-n})`. If ``None``, ``"_" + self.ambient()._name``
          is used
        - ``latex_postscript`` -- (default: ``None``) string defining the LaTeX
          name of the coordinates of the adapted chart. This string will be
          appended to the LaTeX names of the coordinates `(x_1,\ldots,x_n)` and
          of the parameters `(t_1,\ldots,t_{m-n})`, If ``None``,
          ``"_" + self.ambient()._latex_()`` is used

        OUTPUT: list of adapted charts on `M` created from the charts of ``self``

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='topological',
            ....:              latex_name=r"\mathcal{M}")
            sage: N = Manifold(2, 'N', ambient=M, structure='topological')
            sage: N
            2-dimensional topological submanifold N immersed in the
             3-dimensional topological manifold M
            sage: CM.<x,y,z> = M.chart()
            sage: CN.<u,v> = N.chart()
            sage: t = var('t')
            sage: phi = N.continuous_map(M, {(CN,CM): [u,v,t+u^2+v^2]})
            sage: phi_inv = M.continuous_map(N, {(CM,CN): [x,y]})
            sage: phi_inv_t = M.scalar_field({CM: z-x^2-y^2})
            sage: N.set_embedding(phi, inverse=phi_inv, var=t,
            ....:                 t_inverse={t:phi_inv_t})
            sage: N.adapted_chart()
            [Chart (M, (u_M, v_M, t_M))]
            sage: latex(_)
            \left[\left(\mathcal{M},({{u}_{\mathcal{M}}}, {{v}_{\mathcal{M}}},
             {{t}_{\mathcal{M}}})\right)\right]

        The adapted chart has been added to the atlas of ``M``::

            sage: M.atlas()
            [Chart (M, (x, y, z)), Chart (M, (u_M, v_M, t_M))]
            sage: N.atlas()
            [Chart (N, (u, v))]

        The names of the adapted coordinates can be customized::

            sage: N.adapted_chart(postscript='1', latex_postscript='_1')
            [Chart (M, (u1, v1, t1))]
            sage: latex(_)
            \left[\left(\mathcal{M},({{u}_1}, {{v}_1}, {{t}_1})\right)\right]
        """
        if not self._embedded:
            raise ValueError("an embedding is required")

        if self._dim_foliation + self._dim != self._ambient._dim:
            raise ValueError("a foliation of dimension dim(M) - dim(N) is "
                             "needed to find an adapted chart")
        res = []
        self._subs = []

        if postscript is None:
            postscript = "_" + self._ambient._name.replace("^", "")
            # NB: "^" is deleted from the name of ambient to get valid
            # Python identifiers for the symbolic variables representing the
            # coordinates
        if latex_postscript is None:
            latex_postscript = "_{" + self._ambient._latex_() + "}"

        # All possible expressions for the immersion
        chart_pairs = list(self._immersion._coord_expression.keys())
        for (chart1, chart2) in chart_pairs:
            name = " ".join(chart1[i]._repr_() + postscript + ":{"
                             + chart1[i]._latex_() + "}" + latex_postscript
                             for i in self.irange()) + " " \
                   + " ".join(v._repr_() + postscript + ":{" + v._latex_()
                              + "}" + latex_postscript for v in self._var)
            chart = chart2.domain().chart(name)
            if chart not in res:

                # Construct restrictions on coordinates:
                subs = {chart1[:][i]: chart[:][i] for i in range(self._dim)}
                # NB: chart1[:][i] is used instead of chart1[i] to allow for
                #     start_index != 0
                for i in range(len(self._var)):
                    subs[self._var[i]] = chart[:][self._dim + i]
                for rest in chart1._restrictions:
                    chart.add_restrictions(rest.subs(subs))
                for _a in assumptions(*(chart1[:] + tuple(self._var))):
                    if isinstance(_a, Expression):
                        assume(_a.subs(subs))

                self._subs.append(subs)
                res.append(chart)
                self._immersion.add_expr(chart1, chart,
                                         list(chart1[:]) + self._var)
                self._immersion_inv.add_expr(chart, chart1,
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
          appearing in the foliation
        - ``u`` -- iterable of the values taken by the first coordinate of the
          surface to plot
        - ``v`` -- iterable of the values taken by the second coordinate of the
          surface to plot
        - ``chart1`` -- (default: ``None``) chart in which ``u`` and ``v`` are
          considered. By default, the default chart of the submanifold is used
        - ``chart2`` -- (default: ``None``) chart in the codomain of the
          embedding. By default, the default chart of the codomain is used
        - ``**kwargs`` -- other arguments as used in
          :class:`~sage.plot.plot3d.parametric_surface.ParametricSurface`

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: N = Manifold(2, 'N', ambient = M, structure='topological')
            sage: CM.<x,y,z> = M.chart()
            sage: CN.<u,v> = N.chart()
            sage: t = var('t')
            sage: phi = N.continuous_map(M, {(CN,CM): [u,v,t+u^2+v^2]})
            sage: phi_inv = M.continuous_map(N, {(CM,CN): [x,y]})
            sage: phi_inv_t = M.scalar_field({CM: z-x^2-y^2})
            sage: N.set_embedding(phi, inverse=phi_inv, var=t,
            ....:                 t_inverse = {t:phi_inv_t})
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
            sage: P0 + P1 + P2 + P3
            Graphics3d Object

        .. PLOT::

            M = Manifold(3, 'M', structure='topological')
            N = Manifold(2, 'N', ambient = M, structure='topological')
            CM = M.chart('x y z'); x, y, z = CM[:]
            CN = N.chart('u v'); u, v = CN[:]
            t = var('t')
            phi = N.continuous_map(M, {(CN,CM): [u,v,t+u**2+v**2]})
            phi_inv = M.continuous_map(N, {(CM,CN): [x,y]})
            phi_inv_t = M.scalar_field({CM: z-x**2-y**2})
            N.set_embedding(phi, inverse=phi_inv, var=t,
                            t_inverse = {t:phi_inv_t})
            N.adapted_chart()
            P0 = N.plot({t:0}, srange(-1, 1, 0.1), srange(-1, 1, 0.1),
                        CN, CM, opacity=0.3, mesh=True)
            P1 = N.plot({t:1}, srange(-1, 1, 0.1), srange(-1, 1, 0.1),
                        CN, CM, opacity=0.3, mesh=True)
            P2 = N.plot({t:2}, srange(-1, 1, 0.1), srange(-1, 1, 0.1),
                        CN, CM, opacity=0.3, mesh=True)
            P3 = N.plot({t:3}, srange(-1, 1, 0.1), srange(-1, 1, 0.1),
                        CN, CM, opacity=0.3, mesh=True)
            sphinx_plot(P0 + P1 + P2 + P3)

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

    def ambient(self) -> TopologicalManifold:
        r"""
        Return the manifold in which ``self`` is immersed or embedded.

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: N = Manifold(2, 'N', ambient=M, structure='topological')
            sage: N.ambient()
            3-dimensional topological manifold M
        """
        return self._ambient

    def immersion(self) -> ContinuousMap:
        r"""
        Return the immersion of ``self`` into the ambient manifold.

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: N = Manifold(2, 'N', ambient=M, structure='topological')
            sage: CM.<x,y,z> = M.chart()
            sage: CN.<u,v> = N.chart()
            sage: t = var('t')
            sage: phi = N.continuous_map(M, {(CN,CM): [u,v,t+u^2+v^2]})
            sage: phi_inv = M.continuous_map(N, {(CM,CN): [x,y]})
            sage: phi_inv_t = M.scalar_field({CM: z-x^2-y^2})
            sage: N.set_immersion(phi, inverse=phi_inv, var=t,
            ....:                 t_inverse={t: phi_inv_t})
            sage: N.immersion()
            Continuous map from the 2-dimensional topological submanifold N
             immersed in the 3-dimensional topological manifold M to the
             3-dimensional topological manifold M
        """
        if not self._immersed:
            raise ValueError("the submanifold is not immersed")
        assert self._immersion
        return self._immersion

    def embedding(self) -> ContinuousMap:
        r"""
        Return the embedding of ``self`` into the ambient manifold.

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: N = Manifold(2, 'N', ambient=M, structure='topological')
            sage: CM.<x,y,z> = M.chart()
            sage: CN.<u,v> = N.chart()
            sage: t = var('t')
            sage: phi = N.continuous_map(M, {(CN,CM): [u,v,t+u^2+v^2]})
            sage: phi_inv = M.continuous_map(N, {(CM,CN): [x,y]})
            sage: phi_inv_t = M.scalar_field({CM: z-x^2-y^2})
            sage: N.set_embedding(phi, inverse=phi_inv, var=t,
            ....:                 t_inverse={t: phi_inv_t})
            sage: N.embedding()
            Continuous map from the 2-dimensional topological submanifold N
             embedded in the 3-dimensional topological manifold M to the
             3-dimensional topological manifold M
        """
        if not self._embedded:
            raise ValueError("the submanifold is not embedded")
        assert self._immersion
        return self._immersion

    def as_subset(self):
        r"""
        Return ``self`` as a subset of the ambient manifold.

        ``self`` must be an embedded submanifold.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: N = Manifold(1, 'N', ambient=M, structure='topological')
            sage: CM.<x,y> = M.chart()
            sage: CN.<u> = N.chart(coord_restrictions=lambda u: [u > -1, u < 1])
            sage: phi = N.continuous_map(M, {(CN,CM): [u, u^2]})
            sage: N.set_embedding(phi)
            sage: N
            1-dimensional topological submanifold N
              embedded in the 2-dimensional topological manifold M
            sage: N.as_subset()
            Image of the Continuous map
              from the 1-dimensional topological submanifold N
                embedded in the 2-dimensional topological manifold M
              to the 2-dimensional topological manifold M
        """
        return self.embedding().image()
