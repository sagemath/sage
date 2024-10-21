# sage_setup: distribution = sagemath-symbolics
r"""
Submanifolds of differentiable manifolds

Given two differentiable manifolds `N` and `M`, an *immersion* `\phi` is a
differentiable map `N\to M` whose differential is everywhere
injective. One then says that `N` is an *immersed submanifold* of `M`, via
`\phi`.

If in addition, `\phi` is a differentiable embedding (i.e. `\phi` is an
immersion that is a homeomorphism onto its image), then `N` is called an
*embedded submanifold* of `M` (or simply a *submanifold*).

`\phi` can also depend on one or multiple parameters. As long as the
differential of `\phi` remains injective in these parameters, it represents a
*foliation*. The *dimension* of the foliation is defined as the number of
parameters.

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

from sage.manifolds.differentiable.manifold import DifferentiableManifold
from sage.manifolds.topological_submanifold import TopologicalSubmanifold
from sage.rings.infinity import infinity


class DifferentiableSubmanifold(DifferentiableManifold, TopologicalSubmanifold):
    r"""
    Submanifold of a differentiable manifold.

    Given two differentiable manifolds `N` and `M`, an *immersion* `\phi` is a
    differentiable map `N\to M` whose differential is everywhere
    injective. One then says that `N` is an *immersed submanifold* of `M`, via
    `\phi`.

    If in addition, `\phi` is a differentiable embedding (i.e. `\phi` is an
    immersion that is a homeomorphism onto its image), then `N` is called an
    *embedded submanifold* of `M` (or simply a *submanifold*).

    `\phi` can also depend on one or multiple parameters. As long as the
    differential of `\phi` remains injective in these parameters, it represents
    a *foliation*. The *dimension* of the foliation is defined as the number of
    parameters.

    INPUT:

    - ``n`` -- positive integer; dimension of the submanifold
    - ``name`` -- string; name (symbol) given to the submanifold
    - ``field`` -- field `K` on which the sub manifold is defined; allowed
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
      must be a differentiable manifold. If ``None``, it is set to ``self``
    - ``base_manifold`` -- (default: ``None``) if not ``None``, must be a
      differentiable manifold; the created object is then an open subset of
      ``base_manifold``
    - ``diff_degree`` -- (default: ``infinity``) degree of differentiability
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the submanifold; if none are provided, it is set to ``name``
    - ``start_index`` -- (default: 0) integer; lower value of the range of
      indices used for "indexed objects" on the submanifold, e.g., coordinates
      in a chart
    - ``category`` -- (default: ``None``) to specify the category; if ``None``,
      ``Manifolds(field).Differentiable()`` (or ``Manifolds(field).Smooth()``
      if ``diff_degree`` = ``infinity``) is assumed (see the category
      :class:`~sage.categories.manifolds.Manifolds`)
    - ``unique_tag`` -- (default: ``None``) tag used to force the construction
      of a new object when all the other arguments have been used previously
      (without ``unique_tag``, the
      :class:`~sage.structure.unique_representation.UniqueRepresentation`
      behavior inherited from
      :class:`~sage.manifolds.subset.ManifoldSubset` via
      :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`
      would return the previously constructed object corresponding to these
      arguments)

    EXAMPLES:

    Let `N` be a 2-dimensional submanifold of a 3-dimensional manifold `M`::

        sage: M = Manifold(3, 'M')
        sage: N = Manifold(2, 'N', ambient=M)
        sage: N
        2-dimensional differentiable submanifold N immersed in the
         3-dimensional differentiable manifold M
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
    that are adapted to the foliation, ie in which the expression of the
    immersion is trivial. At the same time, the appropriate coordinate changes
    are computed::

        sage: N.adapted_chart()
        [Chart (M, (u_M, v_M, t_M))]
        sage: M.atlas()
        [Chart (M, (x, y, z)), Chart (M, (u_M, v_M, t_M))]
        sage: len(M.coord_changes())
        2

    .. SEEALSO::

        :mod:`~sage.manifolds.manifold` and
        :mod:`~sage.manifolds.topological_submanifold`
    """
    def __init__(self, n, name, field, structure, ambient=None,
                 base_manifold=None, diff_degree=infinity,
                 latex_name=None, start_index=0, category=None,
                 unique_tag=None):
        r"""
        Construct a submanifold of a differentiable manifold.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: N = Manifold(2, 'N', ambient=M)
            sage: N
            2-dimensional differentiable submanifold N immersed in the
             3-dimensional differentiable manifold M
            sage: S = Manifold(2, 'S', latex_name=r'\Sigma', ambient=M,
            ....:              start_index=1)
            sage: latex(S)
            \Sigma
            sage: S.start_index()
            1
        """
        DifferentiableManifold.__init__(self, n, name, field, structure,
                                        base_manifold=base_manifold,
                                        diff_degree=diff_degree,
                                        latex_name=latex_name,
                                        start_index=start_index,
                                        category=category)
        if not (ambient is None
                or isinstance(ambient, DifferentiableManifold)):
            raise TypeError("ambient must be a differentiable manifold")
        self._init_immersion(ambient=ambient)

    def _repr_(self):
        r"""
        Return a string representation of the submanifold.

        If no ambient manifold is specified, the submanifold is considered as
        a manifold.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: N = Manifold(2, 'N', ambient=M)
            sage: N
            2-dimensional differentiable submanifold N immersed in the
             3-dimensional differentiable manifold M
            sage: phi = N.diff_map(M)
            sage: N.set_embedding(phi)
            sage: N
            2-dimensional differentiable submanifold N embedded in the
             3-dimensional differentiable manifold M
        """
        if self is not self._manifold:
            return "Open subset {} of the {}".format(self._name, self._manifold)
        if self._ambient is None:
            return super(DifferentiableManifold, self).__repr__()
        if self._embedded:
            return "{}-dimensional {} submanifold {} embedded in the {}".format(
                self._dim, self._structure.name, self._name, self._ambient)
        return "{}-dimensional {} submanifold {} immersed in the {}".format(
                self._dim, self._structure.name, self._name, self._ambient)

    def open_subset(self, name, latex_name=None, coord_def={}, supersets=None):
        r"""
        Create an open subset of the manifold.

        An open subset is a set that is (i) included in the manifold and (ii)
        open with respect to the manifold's topology. It is a differentiable
        manifold by itself.

        As ``self`` is a submanifold of its ambient manifold,
        the new open subset is also considered a submanifold of that.
        Hence the returned object is an instance of
        :class:`DifferentiableSubmanifold`.

        INPUT:

        - ``name`` -- name given to the open subset
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          subset; if none is provided, it is set to ``name``
        - ``coord_def`` -- (default: {}) definition of the subset in
          terms of coordinates; ``coord_def`` must a be dictionary with keys
          charts in the manifold's atlas and values the symbolic expressions
          formed by the coordinates to define the subset.
        - ``supersets`` -- (default: only ``self``) list of sets that the
          new open subset is a subset of

        OUTPUT: the open subset, as an instance of :class:`DifferentiableSubmanifold`

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='differentiable')
            sage: N = Manifold(2, 'N', ambient=M, structure='differentiable'); N
            2-dimensional differentiable submanifold N immersed in the
             3-dimensional differentiable manifold M
            sage: S = N.subset('S'); S
            Subset S of the
             2-dimensional differentiable submanifold N immersed in the
              3-dimensional differentiable manifold M
            sage: O = N.subset('O', is_open=True); O  # indirect doctest
            Open subset O of the
             2-dimensional differentiable submanifold N immersed in the
              3-dimensional differentiable manifold M

            sage: phi = N.diff_map(M)
            sage: N.set_embedding(phi)
            sage: N
            2-dimensional differentiable submanifold N embedded in the
             3-dimensional differentiable manifold M
            sage: S = N.subset('S'); S
            Subset S of the
             2-dimensional differentiable submanifold N embedded in the
              3-dimensional differentiable manifold M
            sage: O = N.subset('O', is_open=True); O  # indirect doctest
            Open subset O of the
             2-dimensional differentiable submanifold N embedded in the
              3-dimensional differentiable manifold M
        """
        resu = DifferentiableSubmanifold(self._dim, name, self._field,
                                         self._structure, ambient=self._ambient,
                                         base_manifold=self._manifold,
                                         diff_degree=self._diff_degree,
                                         latex_name=latex_name,
                                         start_index=self._sindex)
        if supersets is None:
            supersets = [self]
        for superset in supersets:
            superset._init_open_subset(resu, coord_def=coord_def)
        return resu
