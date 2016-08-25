r"""
Sets of Morphisms between Differentiable Manifolds

The class :class:`DifferentiableManifoldHomset` implements sets of morphisms between
two differentiable manifolds over the same topological field `K` (in most
applications, `K = \RR` or `K = \CC`), a morphism being a *differentiable map*
for the category of differentiable manifolds.

AUTHORS:

- Eric Gourgoulhon (2015): initial version

REFERENCES:

- [Lee13]_ \J.M. Lee : *Introduction to Smooth Manifolds*,
  2nd ed., Springer (New York) (2013)
- [KN63]_ \S. Kobayashi & K. Nomizu : *Foundations of Differential
  Geometry*, vol. 1, Interscience Publishers (New York) (1963)

"""
#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.manifolds.manifold_homset import TopologicalManifoldHomset
from sage.manifolds.differentiable.diff_map import DiffMap

class DifferentiableManifoldHomset(TopologicalManifoldHomset):
    r"""
    Set of differentiable maps between two differentiable manifolds.

    Given two differentiable manifolds `M` and `N` over a topological field `K`,
    the class :class:`DifferentiableManifoldHomset` implements the set
    `\mathrm{Hom}(M,N)` of morphisms (i.e. differentiable maps)
    `M\rightarrow N`.

    This is a Sage *parent* class, whose *element* class is
    :class:`~sage.manifolds.differentiable.diff_map.DiffMap`.

    INPUT:

    - ``domain`` -- differentiable manifold `M` (domain of the morphisms),
      as an instance of
      :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`
    - ``codomain`` -- differentiable manifold `N` (codomain of the morphisms),
      as an instance of
      :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`
    - ``name`` -- (default: ``None``) string; name given to the homset; if
      ``None``, Hom(M,N) will be used
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
      homset; if ``None``, `\mathrm{Hom}(M,N)` will be used

    EXAMPLES:

    Set of differentiable maps between a 2-dimensional differentiable manifold
    and a 3-dimensional one::

        sage: M = Manifold(2, 'M')
        sage: X.<x,y> = M.chart()
        sage: N = Manifold(3, 'N')
        sage: Y.<u,v,w> = N.chart()
        sage: H = Hom(M, N) ; H
        Set of Morphisms from 2-dimensional differentiable manifold M to
         3-dimensional differentiable manifold N in Category of smooth
         manifolds over Real Field with 53 bits of precision
        sage: type(H)
        <class 'sage.manifolds.differentiable.manifold_homset.DifferentiableManifoldHomset_with_category'>
        sage: H.category()
        Category of homsets of topological spaces
        sage: latex(H)
        \mathrm{Hom}\left(M,N\right)
        sage: H.domain()
        2-dimensional differentiable manifold M
        sage: H.codomain()
        3-dimensional differentiable manifold N

    An element of ``H`` is a differentiable map from ``M`` to ``N``::

        sage: H.Element
        <class 'sage.manifolds.differentiable.diff_map.DiffMap'>
        sage: f = H.an_element() ; f
        Differentiable map from the 2-dimensional differentiable manifold M to the
         3-dimensional differentiable manifold N
        sage: f.display()
        M --> N
           (x, y) |--> (u, v, w) = (0, 0, 0)

    The test suite is passed::

        sage: TestSuite(H).run()

    When the codomain coincides with the domain, the homset is a set of
    *endomorphisms* in the category of differentiable manifolds::

        sage: E = Hom(M, M) ; E
        Set of Morphisms from 2-dimensional differentiable manifold M to
         2-dimensional differentiable manifold M in Category of smooth
         manifolds over Real Field with 53 bits of precision
        sage: E.category()
        Category of endsets of topological spaces
        sage: E.is_endomorphism_set()
        True
        sage: E is End(M)
        True

    In this case, the homset is a monoid for the law of morphism composition::

        sage: E in Monoids()
        True

    This was of course not the case for ``H = Hom(M, N)``::

        sage: H in Monoids()
        False

    The identity element of the monoid is of course the identity map of ``M``::

        sage: E.one()
        Identity map Id_M of the 2-dimensional differentiable manifold M
        sage: E.one() is M.identity_map()
        True
        sage: E.one().display()
        Id_M: M --> M
           (x, y) |--> (x, y)

    The test suite is passed by ``E``::

        sage: TestSuite(E).run()

    This test suite includes more tests than in the case of ``H``, since ``E``
    has some extra structure (monoid).

    """

    Element = DiffMap

    def __init__(self, domain, codomain, name=None, latex_name=None):
        r"""
        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: N = Manifold(3, 'N')
            sage: Y.<u,v,w> = N.chart()
            sage: H = Hom(M, N) ; H
            Set of Morphisms from 2-dimensional differentiable manifold M to
             3-dimensional differentiable manifold N in Category of smooth
             manifolds over Real Field with 53 bits of precision
            sage: TestSuite(H).run()

        Test for an endomorphism set::

            sage: E = Hom(M, M) ; E
            Set of Morphisms from 2-dimensional differentiable manifold M to
             2-dimensional differentiable manifold M in Category of smooth
             manifolds over Real Field with 53 bits of precision
            sage: TestSuite(E).run()

        """
        from sage.manifolds.differentiable.manifold import \
                                                         DifferentiableManifold
        if not isinstance(domain, DifferentiableManifold):
            raise TypeError("domain = {} is not an ".format(domain) +
                            "instance of DifferentiableManifold")
        if not isinstance(codomain, DifferentiableManifold):
            raise TypeError("codomain = {} is not an ".format(codomain) +
                            "instance of DifferentiableManifold")
        TopologicalManifoldHomset.__init__(self, domain, codomain, name=name,
                                   latex_name=latex_name)

    #### Parent methods ####

    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to ``self`` exists from other parent.

        EXAMPLE::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: N = Manifold(3, 'N')
            sage: Y.<u,v,w> = N.chart()
            sage: H = Hom(M,N)
            sage: H._coerce_map_from_(ZZ)
            False
            sage: H._coerce_map_from_(M)
            False
            sage: H._coerce_map_from_(N)
            False

        """
        #!# for the time being:
        return False

    #### End of parent methods ####
