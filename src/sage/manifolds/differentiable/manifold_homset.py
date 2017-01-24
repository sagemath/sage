r"""
Sets of Morphisms between Differentiable Manifolds

The class :class:`DifferentiableManifoldHomset` implements sets of morphisms between
two differentiable manifolds over the same topological field `K` (in most
applications, `K = \RR` or `K = \CC`), a morphism being a *differentiable map*
for the category of differentiable manifolds.

The subclass :class:`DifferentiableCurveSet` is devoted to the specific case of
differential curves, i.e. morphisms whose domain is an open interval of
`\RR`.

AUTHORS:

- Eric Gourgoulhon (2015): initial version
- Travis Scrimshaw (2016): review tweaks

REFERENCES:

- [Lee2013]_
- [KN1963]_

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
from sage.manifolds.differentiable.curve import DifferentiableCurve

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


#******************************************************************************

class DifferentiableCurveSet(DifferentiableManifoldHomset):
    r"""
    Set of differentiable curves in a differentiable manifold.

    Given an open interval `I` of `\RR` (possibly `I = \RR`) and
    a differentiable manifold `M` over `\RR`, this is the set
    `\mathrm{Hom}(I,M)` of morphisms (i.e. differentiable curves) `I \to M`.

    INPUT:

    - ``domain`` --
      :class:`~sage.manifolds.differentiable.real_line.OpenInterval`
      if an open interval `I \subset \RR` (domain of the morphisms),
      or :class:`~sage.manifolds.differentiable.real_line.RealLine`
      if `I = \RR`
    - ``codomain`` --
      :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`;
      differentiable manifold `M` (codomain of the morphisms)
    - ``name`` -- (default: ``None``) string; name given to the set of
      curves; if ``None``, ``Hom(I, M)`` will be used
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote
      the set of curves; if ``None``, `\mathrm{Hom}(I,M)` will be used

    EXAMPLES:

    Set of curves `\RR \longrightarrow M`, where `M` is a 2-dimensional
    manifold::

        sage: M = Manifold(2, 'M')
        sage: X.<x,y> = M.chart()
        sage: R.<t> = RealLine() ; R
        Real number line R
        sage: H = Hom(R, M) ; H
        Set of Morphisms from Real number line R to 2-dimensional
         differentiable manifold M in Category of smooth manifolds over Real
         Field with 53 bits of precision
        sage: H.category()
        Category of homsets of topological spaces
        sage: latex(H)
        \mathrm{Hom}\left(\Bold{R},M\right)
        sage: H.domain()
        Real number line R
        sage: H.codomain()
        2-dimensional differentiable manifold M

    An element of ``H`` is a curve in ``M``::

        sage: c = H.an_element(); c
        Curve in the 2-dimensional differentiable manifold M
        sage: c.display()
        R --> M
           t |--> (x, y) = (1/(t^2 + 1) - 1/2, 0)

    The test suite is passed::

        sage: TestSuite(H).run()

    The set of curves `(0,1) \longrightarrow U`, where `U` is an open
    subset of `M`::

        sage: I = R.open_interval(0, 1) ; I
        Real interval (0, 1)
        sage: U = M.open_subset('U', coord_def={X: x^2+y^2<1}) ; U
        Open subset U of the 2-dimensional differentiable manifold M
        sage: H = Hom(I, U) ; H
        Set of Morphisms from Real interval (0, 1) to Open subset U of the
         2-dimensional differentiable manifold M in Join of Category of
         subobjects of sets and Category of smooth manifolds over Real Field
         with 53 bits of precision

    An element of ``H`` is a curve in ``U``::

        sage: c = H.an_element() ; c
        Curve in the Open subset U of the 2-dimensional differentiable
         manifold M
        sage: c.display()
        (0, 1) --> U
           t |--> (x, y) = (1/(t^2 + 1) - 1/2, 0)

    The set of curves `\RR \longrightarrow \RR` is a set of (manifold)
    endomorphisms::

        sage: E = Hom(R, R) ; E
        Set of Morphisms from Real number line R to Real number line R in
         Category of smooth manifolds over Real Field with 53 bits of precision
        sage: E.category()
        Category of endsets of topological spaces
        sage: E.is_endomorphism_set()
        True
        sage: E is End(R)
        True

    It is a monoid for the law of morphism composition::

        sage: E in Monoids()
        True

    The identity element of the monoid is the identity map of `\RR`::

        sage: E.one()
        Identity map Id_R of the Real number line R
        sage: E.one() is R.identity_map()
        True
        sage: E.one().display()
        Id_R: R --> R
           t |--> t

    A "typical" element of the monoid::

        sage: E.an_element().display()
        R --> R
           t |--> 1/(t^2 + 1) - 1/2

    The test suite is passed by ``E``::

        sage: TestSuite(E).run()

    Similarly, the set of curves `I \longrightarrow I` is a monoid, whose
    elements are (manifold) endomorphisms::

        sage: EI = Hom(I, I) ; EI
        Set of Morphisms from Real interval (0, 1) to Real interval (0, 1) in
         Join of Category of subobjects of sets and
             Category of smooth manifolds over Real Field with 53 bits of precision
        sage: EI.category()
        Category of endsets of subobjects of sets and topological spaces
        sage: EI is End(I)
        True
        sage: EI in Monoids()
        True

    The identity element and a "typical" element of this monoid::

        sage: EI.one()
        Identity map Id_(0, 1) of the Real interval (0, 1)
        sage: EI.one().display()
        Id_(0, 1): (0, 1) --> (0, 1)
           t |--> t
        sage: EI.an_element().display()
        (0, 1) --> (0, 1)
           t |--> 1/2/(t^2 + 1) + 1/4

    The test suite is passed by ``EI``::

        sage: TestSuite(EI).run()

    """
    Element = DifferentiableCurve

    def __init__(self, domain, codomain, name=None, latex_name=None):
        r"""
        Initialize ``self``.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: R.<t> = RealLine()
            sage: H = Hom(R, M); H
            Set of Morphisms from Real number line R to 3-dimensional
             differentiable manifold M in Category of smooth manifolds over
             Real Field with 53 bits of precision
            sage: TestSuite(H).run()
            sage: Hom(R, M) is Hom(R, M)
            True
            sage: H = Hom(R, R); H
            Set of Morphisms from Real number line R to Real number line R in
             Category of smooth manifolds over Real Field with 53 bits of
             precision
            sage: TestSuite(H).run()
            sage: I = R.open_interval(-1, 2)
            sage: H = Hom(I, M); H
            Set of Morphisms from Real interval (-1, 2) to 3-dimensional
             differentiable manifold M in Join of Category of subobjects of
             sets and Category of smooth manifolds over Real Field with 53 bits
             of precision
            sage: TestSuite(H).run()
            sage: H = Hom(I, I); H
            Set of Morphisms from Real interval (-1, 2) to Real interval (-1, 2)
             in Join of Category of subobjects of sets and Category of smooth
             manifolds over Real Field with 53 bits of precision
            sage: TestSuite(H).run()

        """
        from sage.manifolds.differentiable.real_line import OpenInterval
        if not isinstance(domain, OpenInterval):
            raise TypeError("{} is not an open real interval".format(domain))
        DifferentiableManifoldHomset.__init__(self, domain, codomain, name=name,
                                              latex_name=latex_name)

    #### Parent methods ####

    def _element_constructor_(self, coord_expression, name=None,
                              latex_name=None, is_isomorphism=False,
                              is_identity=False):
        r"""
        Construct an element of ``self``, i.e. a differentiable curve
        `I \to M`, where `I` is a real interval and `M` some
        differentiable manifold.

        OUTPUT:

        - :class:`~sage.manifolds.differentiable.curve.DifferentiableCurve`

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: R.<t> = RealLine() ; R
            Real number line R
            sage: H = Hom(R, M)
            sage: c = H({X: [sin(t), sin(2*t)/2]}, name='c') ; c
            Curve c in the 2-dimensional differentiable manifold M
            sage: c = Hom(R, R)({}, is_identity=True) ; c
            Identity map Id_R of the Real number line R

        """
        # Standard construction
        return self.element_class(self, coord_expression=coord_expression,
                                  name=name, latex_name=latex_name,
                                  is_isomorphism=is_isomorphism,
                                  is_identity=is_identity)

    def _an_element_(self):
        r"""
        Construct some element of ``self``.

        OUTPUT:

        - :class:`~sage.manifolds.differentiable.curve.DifferentiableCurve`

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: R.<t> = RealLine()
            sage: c = Hom(R,M)._an_element_() ; c
            Curve in the 3-dimensional differentiable manifold M
            sage: c.display()
            R --> M
               t |--> (x, y, z) = (1/(t^2 + 1) - 1/2, 0, 0)

        ::

            sage: I = R.open_interval(0, pi)
            sage: c = Hom(I,M)._an_element_() ; c
            Curve in the 3-dimensional differentiable manifold M
            sage: c.display()
            (0, pi) --> M
               t |--> (x, y, z) = (1/(t^2 + 1) - 1/2, 0, 0)

        ::

            sage: c = Hom(I,I)._an_element_() ; c
            Differentiable map from the Real interval (0, pi) to itself
            sage: c.display()
            (0, pi) --> (0, pi)
               t |--> 1/4*pi + 1/2*pi/(t^2 + 1)

        """
        from sage.rings.infinity import Infinity
        from sage.rings.rational_field import QQ
        dom = self.domain()
        codom = self.codomain()
        # A simple curve is constructed around a point of the codomain:
        chart2 = codom.default_chart()
        target_point = chart2.domain().an_element()
        target_coord = list(target_point.coord(chart2))
        bounds = chart2._bounds[0]  # bounds of first coordinate
        # Determination of an interval (x1, x2) arround target_point:
        xmin = bounds[0][0]
        xmax = bounds[1][0]
        one_half = QQ(1) / QQ(2)
        if xmin == -Infinity:
            if xmax == Infinity:
                x1 = - one_half
                x2 = one_half
            else:
                x1 = xmax - 3*one_half
                x2 = xmax - one_half
        else:
            if xmax == Infinity:
                x1 = xmin + one_half
                x2 = xmin + 3*one_half
            else:
                dx = (xmax - xmin) / 4
                x1 = xmin + dx
                x2 = xmax - dx
        # The coordinate function defining the curve:
        t = dom.canonical_coordinate()
        target_coord[0] = x1 + (x2-x1) / (1+t*t)
        coord_expression = {chart2: target_coord}
        return self.element_class(self, coord_expression)

