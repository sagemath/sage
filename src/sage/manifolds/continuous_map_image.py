r"""
Images of Manifold Subsets under Continuous Maps as Subsets of the Codomain

:class:`ImageManifoldSubset` implements the image of a continuous map `\Phi`
from a manifold `M` to some manifold `N` as a subset `\Phi(M)` of `N`,
or more generally, the image `\Phi(S)` of a subset `S \subseteq M` as a
subset of `N`.
"""

# ****************************************************************************
#       Copyright (C) 2021 Matthias Koeppe <mkoeppe@math.ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.manifolds.subset import ManifoldSubset

class ImageManifoldSubset(ManifoldSubset):
    r"""
    Subset of a topological manifold that is a continuous image of a manifold subset.

    INPUT:

    - ``map`` -- continuous map `\Phi`
    - ``inverse`` -- (default: ``None``) continuous map from
      ``map.codomain()`` to ``map.domain()``, which once restricted to the image
      of `\Phi` is the inverse of `\Phi` onto its image if the latter
      exists (NB: no check of this is performed)
    - ``name`` -- (default: computed from the names of the map and the subset)
      string; name (symbol) given to the subset
    - ``latex_name`` -- string (default: ``None``); LaTeX symbol to
      denote the subset; if none is provided, it is set to ``name``
    - ``domain_subset`` -- (default: the domain of ``map``) a subset of the domain of
      ``map``
    """

    def __init__(self, map, inverse=None, name=None, latex_name=None, domain_subset=None):
        r"""
        Construct a manifold subset that is the image of a continuous map.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: N = Manifold(1, 'N', ambient=M, structure='topological')
            sage: CM.<x,y> = M.chart()
            sage: CN.<u> = N.chart(coord_restrictions=lambda u: [u > -1, u < 1])
            sage: Phi = N.continuous_map(M, {(CN,CM): [u, 1 + u^2]}, name='Phi')
            sage: Phi_inv = M.continuous_map(N, {(CM, CN): [x]}, name='Phi_inv')
            sage: Phi_N = Phi.image(inverse=Phi_inv)
            sage: TestSuite(Phi_N).run()
        """
        self._map = map
        self._inverse = inverse
        if domain_subset is None:
            domain_subset = map.domain()
        self._domain_subset = domain_subset
        base_manifold = map.codomain()
        map_name = map._name or 'f'
        map_latex_name = map._latex_name or map_name
        if latex_name is None:
            if name is None:
                latex_name = map_latex_name + r'(' + domain_subset._latex_name + ')'
            else:
                latex_name = name
        if name is None:
            name = map_name + '_' + domain_subset._name
        ManifoldSubset.__init__(self, base_manifold, name, latex_name=latex_name)

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: N = Manifold(1, 'N', ambient=M, structure='topological')
            sage: CM.<x,y> = M.chart()
            sage: CN.<u> = N.chart(coord_restrictions=lambda u: [u > -1, u < 1])
            sage: Phi = N.continuous_map(M, {(CN,CM): [u, 1 + u^2]}, name='Phi')
            sage: Phi.image()                                 # indirect doctest
            Image of the Continuous map Phi
             from the 1-dimensional topological submanifold N immersed in the
              2-dimensional topological manifold M
             to the 2-dimensional topological manifold M
            sage: S = N.subset('S')
            sage: Phi.image(S)                                # indirect doctest
            Image of the
             Subset S of the
              1-dimensional topological submanifold N immersed in the
               2-dimensional topological manifold M
             under the Continuous map Phi
             from the 1-dimensional topological submanifold N immersed in the
              2-dimensional topological manifold M
             to the 2-dimensional topological manifold M
        """
        if self._domain_subset is self._map.domain():
            return f"Image of the {self._map}"
        else:
            return f"Image of the {self._domain_subset} under the {self._map}"

    def _an_element_(self):
        r"""
        Construct some point in the subset.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: N = Manifold(1, 'N', ambient=M, structure='topological')
            sage: CM.<x,y> = M.chart()
            sage: CN.<u> = N.chart(coord_restrictions=lambda u: [u > -1, u < 1])
            sage: Phi = N.continuous_map(M, {(CN,CM): [u, 1 + u^2]}, name='Phi')
            sage: Phi_N = Phi.image()
            sage: p = Phi_N.an_element(); p   # indirect doctest
            Point on the 2-dimensional topological manifold M
            sage: p.coordinates()
            (0, 1)
        """
        return self._map(self._domain_subset.an_element())

    def __contains__(self, point):
        r"""
        Check whether ``point`` is contained in ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: N = Manifold(1, 'N', ambient=M, structure='topological')
            sage: CM.<x,y> = M.chart()
            sage: CN.<u> = N.chart(coord_restrictions=lambda u: [u > -1, u < 1])
            sage: Phi = N.continuous_map(M, {(CN,CM): [u, 1 + u^2]}, name='Phi')
            sage: Phi_inv = M.continuous_map(N, {(CM, CN): [x]}, name='Phi_inv')
            sage: Phi_N = Phi.image(inverse=Phi_inv)
            sage: M((0, 0)) in Phi_N
            False
            sage: M((0, 1)) in Phi_N
            True
        """
        if super().__contains__(point):
            return True
        if point not in self._map.codomain():
            return False
        if self._inverse is not None:
            preimage = self._inverse(point)
            if preimage not in self._domain_subset:
                return False
            return self._map(preimage) == point
        raise NotImplementedError
