# sage_setup: distribution = sagemath-symbolics
r"""
Subsets of Topological Manifolds

The class :class:`ManifoldSubset` implements generic subsets of a
topological manifold. Open subsets are implemented by the class
:class:`~sage.manifolds.manifold.TopologicalManifold` (since an
open subset of a manifold is a manifold by itself), which inherits
from :class:`ManifoldSubset`.  Besides, subsets that are images of
a manifold subset under a continuous map are implemented by the
subclass :class:`~sage.manifolds.continuous_map_image.ImageManifoldSubset`.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015): initial version
- Travis Scrimshaw (2015): review tweaks; removal of facade parents
- Matthias Koeppe (2021): Families and posets of subsets

REFERENCES:

- [Lee2011]_

EXAMPLES:

Two subsets on a manifold::

    sage: M = Manifold(2, 'M', structure='topological')
    sage: a = M.subset('A'); a
    Subset A of the 2-dimensional topological manifold M
    sage: b = M.subset('B'); b
    Subset B of the 2-dimensional topological manifold M
    sage: M.subset_family()
    Set {A, B, M} of subsets of the 2-dimensional topological manifold M

The intersection of the two subsets::

    sage: c = a.intersection(b); c
    Subset A_inter_B of the 2-dimensional topological manifold M

Their union::

    sage: d = a.union(b); d
    Subset A_union_B of the 2-dimensional topological manifold M

Families of subsets after the above operations::

    sage: M.subset_family()
    Set {A, A_inter_B, A_union_B, B, M} of subsets of the 2-dimensional topological manifold M
    sage: a.subset_family()
    Set {A, A_inter_B} of subsets of the 2-dimensional topological manifold M
    sage: c.subset_family()
    Set {A_inter_B} of subsets of the 2-dimensional topological manifold M
    sage: d.subset_family()
    Set {A, A_inter_B, A_union_B, B} of subsets of the 2-dimensional topological manifold M
"""
#*****************************************************************************
#       Copyright (C) 2015-2020 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015      Michal Bejger <bejger@camk.edu.pl>
#       Copyright (C) 2015-2016 Travis Scrimshaw <tscrimsh@umn.edu>
#       Copyright (C) 2021      Matthias Koeppe <mkoeppe@math.ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import annotations
from typing import Optional
from collections import defaultdict
import itertools
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.superseded import deprecation
from sage.categories.sets_cat import Sets
from sage.manifolds.family import ManifoldObjectFiniteFamily, ManifoldSubsetFiniteFamily
from sage.manifolds.point import ManifoldPoint

class ManifoldSubset(UniqueRepresentation, Parent):
    r"""
    Subset of a topological manifold.

    The class :class:`ManifoldSubset` inherits from the generic
    class :class:`~sage.structure.parent.Parent`.
    The corresponding element class is
    :class:`~sage.manifolds.point.ManifoldPoint`.

    Note that open subsets are not implemented directly by this class, but
    by the derived class :class:`~sage.manifolds.manifold.TopologicalManifold`
    (an open subset of a topological manifold being itself a topological
    manifold).

    INPUT:

    - ``manifold`` -- topological manifold on which the subset is defined
    - ``name`` -- string; name (symbol) given to the subset
    - ``latex_name`` -- string (default: ``None``); LaTeX symbol to
      denote the subset; if none are provided, it is set to ``name``
    - ``category`` -- (default: ``None``) to specify the category;
      if ``None``, the category for generic subsets is used

    EXAMPLES:

    A subset of a manifold::

        sage: M = Manifold(2, 'M', structure='topological')
        sage: from sage.manifolds.subset import ManifoldSubset
        sage: A = ManifoldSubset(M, 'A', latex_name=r'\mathcal{A}')
        sage: A
        Subset A of the 2-dimensional topological manifold M
        sage: latex(A)
        \mathcal{A}
        sage: A.is_subset(M)
        True

    Instead of importing :class:`ManifoldSubset` in the global
    namespace, it is recommended to use the method
    :meth:`~sage.manifolds.subset.ManifoldSubset.subset` to create a new
    subset::

        sage: B = M.subset('B', latex_name=r'\mathcal{B}'); B
        Subset B of the 2-dimensional topological manifold M
        sage: M.subset_family()
        Set {A, B, M} of subsets of the 2-dimensional topological manifold M

    The manifold is itself a subset::

        sage: isinstance(M, ManifoldSubset)
        True
        sage: M in M.subsets()
        True

    Instances of :class:`ManifoldSubset` are parents::

        sage: isinstance(A, Parent)
        True
        sage: A.category()
        Category of subobjects of sets
        sage: p = A.an_element(); p
        Point on the 2-dimensional topological manifold M
        sage: p.parent()
        Subset A of the 2-dimensional topological manifold M
        sage: p in A
        True
        sage: p in M
        True
    """

    Element = ManifoldPoint

    _name: str

    def __init__(self, manifold, name: str, latex_name=None, category=None):
        r"""
        Construct a manifold subset.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A'); A
            Subset A of the 2-dimensional topological manifold M
            sage: type(A)
            <class 'sage.manifolds.subset.ManifoldSubset_with_category'>
            sage: A.category()
            Category of subobjects of sets
            sage: TestSuite(A).run(skip='_test_elements')

        .. NOTE::

            ``_test_elements`` cannot be passed without a proper
            coordinate definition of the subset.
        """
        if not isinstance(name, str):
            raise TypeError("{} is not a string".format(name))
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            if not isinstance(latex_name, str):
                raise TypeError("{} is not a string".format(latex_name))
            self._latex_name = latex_name
        if category is None:
            category = Sets().Subobjects()
            base = None
        else:
            base = manifold._field
        Parent.__init__(self, base=base, category=category)
        if self is not manifold:
            for dom in manifold._subsets:
                if name == dom._name:
                    raise ValueError("the name '" + name +
                                     "' is already used for another " +
                                     "subset of the {}".format(manifold))
            manifold._subsets.add(self)
        self._supersets = set([manifold, self])  # subsets containing self
        self._subsets = set([self])  # subsets of self
        self._top_subsets = set([self])  # subsets contained in self but not
                                         # in another strict subset of self
        self._intersections = {}  # dict. of intersections with other subsets
                                  # (key: subset name)
        self._unions = {}  # dict. of unions with other subsets (key: subset
                           # name)
        self._open_covers = []  # list of open covers of self
        self._is_open = False   # a priori (may be redefined by subclasses)
        self._manifold = manifold  # the ambient manifold
        self._has_defined_points = False

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: A._repr_()
            'Subset A of the 2-dimensional topological manifold M'
            sage: repr(A)  # indirect doctest
            'Subset A of the 2-dimensional topological manifold M'
        """
        return "Subset {} of the {}".format(self._name, self._manifold)

    def _latex_(self):
        r"""
        LaTeX representation of ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: A._latex_()
            'A'
            sage: B = A.subset('B', latex_name=r'\mathcal{B}')
            sage: B._latex_()
            '\\mathcal{B}'
            sage: latex(B)  # indirect doctest
            \mathcal{B}

            sage: M = Manifold(3, 'M', structure='topological')
            sage: M._latex_()
            'M'
            sage: latex(M)
            M
            sage: M = Manifold(3, 'M', latex_name=r'\mathcal{M}',
            ....:              structure='topological')
            sage: M._latex_()
            '\\mathcal{M}'
            sage: latex(M)  # indirect doctest
            \mathcal{M}
        """
        return self._latex_name

    #### Methods required for any Parent in the category of sets:

    def _element_constructor_(self, coords=None, chart=None, name=None,
                              latex_name=None, check_coords=True):
        r"""
        Construct a point in the subset from its coordinates in some chart.

        INPUT:

        - ``coords`` -- (default: ``None``) either (i) the point coordinates
          (as a tuple or a list) in the chart ``chart`` or (ii) another point
          in the subset
        - ``chart`` -- (default: ``None``) chart in which the coordinates
          are given; if none are provided, the coordinates are assumed to
          refer to the subset's default chart
        - ``name`` -- (default: ``None``) name given to the point
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          point; if none are provided, the LaTeX symbol is set to ``name``
        - ``check_coords`` -- boolean (default: ``True``); determines whether
          ``coords`` are valid coordinates for the chart ``chart``.
          For symbolic coordinates, it is recommended to set ``check_coords``
          to ``False``.

        OUTPUT:

        An instance of :class:`~sage.manifolds.point.ManifoldPoint`
        representing a point in the current subset.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: p = M((-2,3)); p  # coord in the default chart
            Point on the 2-dimensional topological manifold M
            sage: X(p)
            (-2, 3)

        A generic subset has no default chart, so the chart must
        be explicitly given::

            sage: A = M.subset('A')
            sage: p = A((-2,3), chart=X); p
            Point on the 2-dimensional topological manifold M
            sage: X(p)
            (-2, 3)
            sage: p.parent()
            Subset A of the 2-dimensional topological manifold M
            sage: p in A
            True

        Coordinates in a chart with some coordinate restrictions::

            sage: Y.<u,v> = M.chart('u:(-1,1) v:(-1,1)')
            sage: p = A((0,1/2), chart=Y); p
            Point on the 2-dimensional topological manifold M
            sage: Y(p)
            (0, 1/2)
            sage: p = A((0,1/2), chart=Y, check_coords=False); p
            Point on the 2-dimensional topological manifold M
            sage: Y(p)
            (0, 1/2)
            sage: p = A((3,1/2), chart=Y)
            Traceback (most recent call last):
            ...
            ValueError: the coordinates (3, 1/2) are not valid on the Chart (M, (u, v))

        Specifying the name of the point::

            sage: p = A((-2,3), chart=X, name='p'); p
            Point p on the 2-dimensional topological manifold M

        A point as entry::

            sage: q = A(p); q
            Point p on the 2-dimensional topological manifold M
            sage: X(q)
            (-2, 3)
        """
        if isinstance(coords, ManifoldPoint):
            point = coords # for readability
            # This should actually never happen by the coercion framework...
            if point.parent() is self:
                return point
            if point in self:
                resu = self.element_class(self, name=point._name,
                                          latex_name=point._latex_name)
                for chart, coords in point._coordinates.items():
                    resu._coordinates[chart] = coords
                return resu
            else:
                raise ValueError("the {}".format(point) +
                                 " is not in {}".format(self))
        return self.element_class(self, coords=coords, chart=chart,
                                  name=name, latex_name=latex_name,
                                  check_coords=check_coords)

    def _an_element_(self):
        r"""
        Construct some point in the subset.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: A = M.subset('A')
            sage: p = A._an_element_(); p
            Point on the 2-dimensional topological manifold M
            sage: p in A
            True
        """
        #!# should be improved...
        return self.element_class(self)

    #### End of methods required for any Parent in the category of sets

    def __contains__(self, point):
        r"""
        Check whether ``point`` is contained in ``self``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: A = M.subset('A')
            sage: p = A((-2,3), chart=X); p
            Point on the 2-dimensional topological manifold M
            sage: A.__contains__(p)
            True
            sage: p in A  # indirect doctest
            True
            sage: A.__contains__(A.an_element())
            True
            sage: q = M((0,0), chart=X); q
            Point on the 2-dimensional topological manifold M
            sage: A.__contains__(q)
            False
        """
        # for efficiency, a quick test first:
        if point.parent() is self:
            return True
        if point.parent().is_subset(self):
            return True
        #!# should be improved once coordinate definition have been introduced
        # in ManifoldSubset
        return False

    def lift(self, p):
        r"""
        Return the lift of ``p`` to the ambient manifold of ``self``.

        INPUT:

        - ``p`` -- point of the subset

        OUTPUT: the same point, considered as a point of the ambient manifold

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: A = M.open_subset('A', coord_def={X: x>0})
            sage: p = A((1, -2)); p
            Point on the 2-dimensional topological manifold M
            sage: p.parent()
            Open subset A of the 2-dimensional topological manifold M
            sage: q = A.lift(p); q
            Point on the 2-dimensional topological manifold M
            sage: q.parent()
            2-dimensional topological manifold M
            sage: q.coord()
            (1, -2)
            sage: (p == q) and (q == p)
            True
        """
        return self._manifold(p)

    def retract(self, p):
        r"""
        Return the retract of ``p`` to ``self``.

        INPUT:

        - ``p`` -- point of the ambient manifold

        OUTPUT: the same point, considered as a point of the subset

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: A = M.open_subset('A', coord_def={X: x>0})
            sage: p = M((1, -2)); p
            Point on the 2-dimensional topological manifold M
            sage: p.parent()
            2-dimensional topological manifold M
            sage: q = A.retract(p); q
            Point on the 2-dimensional topological manifold M
            sage: q.parent()
            Open subset A of the 2-dimensional topological manifold M
            sage: q.coord()
            (1, -2)
            sage: (q == p) and (p == q)
            True

        Of course, if the point does not belong to ``A``, the ``retract``
        method fails::

            sage: p = M((-1, 3))  #  x < 0, so that p is not in A
            sage: q = A.retract(p)
            Traceback (most recent call last):
            ...
            ValueError: the Point on the 2-dimensional topological manifold M
             is not in Open subset A of the 2-dimensional topological manifold M
        """
        return self(p)

    #### Accessors

    def manifold(self):
        r"""
        Return the ambient manifold of ``self``.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: A.manifold()
            2-dimensional topological manifold M
            sage: A.manifold() is M
            True
            sage: B = A.subset('B')
            sage: B.manifold() is M
            True

        An alias is ``ambient``::

            sage: A.ambient() is A.manifold()
            True
        """
        return self._manifold

    ambient = manifold

    def is_open(self):
        r"""
        Return if ``self`` is an open set.

        This method always returns ``False``, since open subsets must be
        constructed as instances of the subclass
        :class:`~sage.manifolds.manifold.TopologicalManifold`
        (which redefines ``is_open``)

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: A.is_open()
            False
        """
        return False

    def is_closed(self):
        r"""
        Return if ``self`` is a closed set.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: M.is_closed()
            True
            sage: also_M = M.subset('also_M')
            sage: M.declare_subset(also_M)
            sage: also_M.is_closed()
            True

            sage: A = M.subset('A')
            sage: A.is_closed()
            False
            sage: A.declare_empty()
            sage: A.is_closed()
            True

            sage: N = M.open_subset('N')
            sage: N.is_closed()
            False
            sage: complement_N = M.subset('complement_N')
            sage: M.declare_union(N, complement_N, disjoint=True)
            sage: complement_N.is_closed()
            True
        """
        if self.manifold().is_subset(self):
            return True
        if self.is_empty():
            return True
        for other_name, intersection in self._intersections.items():
            if intersection.is_empty():
                other = self.manifold().subset_family()[other_name]
                if other.is_open():
                    try:
                        union = self._unions[other_name]
                    except KeyError:
                        pass
                    else:
                        if union.is_open():
                            # self is complement of open other in open union
                            return True
        return False

    def open_covers(self, trivial=True, supersets=False):
        r"""
        Generate the open covers of the current subset.

        If the current subset, `A` say, is a subset of the manifold `M`, an
        *open cover* of `A` is a :class:`ManifoldSubsetFiniteFamily` `F`
        of open subsets `U \in F` of `M` such that

        .. MATH::

            A \subset \bigcup_{U \in F} U.

        If `A` is open, we ask that the above inclusion is actually an
        identity:

        .. MATH::

            A = \bigcup_{U \in F} U.

        .. NOTE::

            To get the open covers as a family, sorted lexicographically by the
            names of the subsets forming the open covers, use the method
            :meth:`open_cover_family` instead.

        INPUT:

        - ``trivial`` -- boolean (default: ``True``); if ``self`` is open,
          include the trivial open cover of ``self`` by itself
        - ``supersets`` -- boolean (default: ``False``); if ``True``, include
          open covers of all the supersets. It can also be an iterable of
          supersets to include.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: M.open_covers()
            <generator ...>
            sage: list(M.open_covers())
            [Set {M} of open subsets of the 2-dimensional topological manifold M]
            sage: U = M.open_subset('U')
            sage: list(U.open_covers())
            [Set {U} of open subsets of the 2-dimensional topological manifold M]
            sage: A = U.open_subset('A')
            sage: B = U.open_subset('B')
            sage: U.declare_union(A,B)
            sage: list(U.open_covers())
            [Set {U} of open subsets of the 2-dimensional topological manifold M,
             Set {A, B} of open subsets of the 2-dimensional topological manifold M]
            sage: list(U.open_covers(trivial=False))
            [Set {A, B} of open subsets of the 2-dimensional topological manifold M]
            sage: V = M.open_subset('V')
            sage: M.declare_union(U,V)
            sage: list(M.open_covers())
            [Set {M} of open subsets of the 2-dimensional topological manifold M,
             Set {U, V} of open subsets of the 2-dimensional topological manifold M,
             Set {A, B, V} of open subsets of the 2-dimensional topological manifold M]
        """
        if supersets is False:
            supersets = [self]
        elif supersets is True:
            supersets = self._supersets
        for superset in supersets:
            for oc in superset._open_covers:
                if not trivial:
                    if any(x in supersets for x in oc):
                        continue
                yield ManifoldSubsetFiniteFamily(oc)

    def open_cover_family(self, trivial=True, supersets=False):
        r"""
        Return the family of open covers of the current subset.

        If the current subset, `A` say, is a subset of the manifold `M`, an
        *open cover* of `A` is a :class:`ManifoldSubsetFiniteFamily` `F`
        of open subsets `U \in F` of `M` such that

        .. MATH::

            A \subset \bigcup_{U \in F} U.

        If `A` is open, we ask that the above inclusion is actually an
        identity:

        .. MATH::

            A = \bigcup_{U \in F} U.

        The family is sorted lexicographically by the names of the subsets
        forming the open covers.

        .. NOTE::

            If you only need to iterate over the open covers in arbitrary
            order, you can use the generator method :meth:`open_covers`
            instead.

        INPUT:

        - ``trivial`` -- boolean (default: ``True``); if ``self`` is open,
          include the trivial open cover of ``self`` by itself
        - ``supersets`` -- boolean (default: ``False``); if ``True``, include
          open covers of all the supersets. It can also be an iterable of
          supersets to include

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: M.open_cover_family()
            Set {{M}} of objects of the 2-dimensional topological manifold M
            sage: U = M.open_subset('U')
            sage: U.open_cover_family()
            Set {{U}} of objects of the 2-dimensional topological manifold M
            sage: A = U.open_subset('A')
            sage: B = U.open_subset('B')
            sage: U.declare_union(A,B)
            sage: U.open_cover_family()
            Set {{A, B}, {U}} of objects of the 2-dimensional topological manifold M
            sage: U.open_cover_family(trivial=False)
            Set {{A, B}} of objects of the 2-dimensional topological manifold M
            sage: V = M.open_subset('V')
            sage: M.declare_union(U,V)
            sage: M.open_cover_family()
            Set {{A, B, V}, {M}, {U, V}} of objects of the 2-dimensional topological manifold M
        """
        return ManifoldObjectFiniteFamily(self.open_covers(
            trivial=trivial, supersets=supersets))

    def open_supersets(self):
        r"""
        Generate the open supersets of ``self``.

        .. NOTE::

            To get the open supersets as a family, sorted by name, use the method
            :meth:`open_superset_family` instead.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: U = M.open_subset('U')
            sage: V = U.subset('V')
            sage: W = V.subset('W')
            sage: sorted(W.open_supersets(), key=lambda S: S._name)
            [2-dimensional topological manifold M,
             Open subset U of the 2-dimensional topological manifold M]
        """
        for superset in self._supersets:
            if superset.is_open():
                yield superset

    def open_superset_family(self):
        r"""
        Return the family of open supersets of ``self``.

        The family is sorted by the alphabetical names of the subsets.

        OUTPUT:

        - a :class:`ManifoldSubsetFiniteFamily` instance containing all the
          open supersets that have been defined on the current subset

        .. NOTE::

            If you only need to iterate over the open supersets in arbitrary
            order, you can use the generator method :meth:`open_supersets`
            instead.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: U = M.open_subset('U')
            sage: V = U.subset('V')
            sage: W = V.subset('W')
            sage: W.open_superset_family()
            Set {M, U} of open subsets of the 2-dimensional topological manifold M
        """
        return ManifoldSubsetFiniteFamily(self.open_supersets())

    def subsets(self):
        r"""
        Generate the subsets that have been defined on the current subset.

        .. NOTE::

            To get the subsets as a family, sorted by name, use the method
            :meth:`subset_family` instead.

        EXAMPLES:

        Subsets of a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: U = M.open_subset('U')
            sage: V = M.subset('V')
            sage: frozenset(M.subsets())  # random (set output)
            {Subset V of the 2-dimensional topological manifold M,
             2-dimensional topological manifold M,
             Open subset U of the 2-dimensional topological manifold M}
            sage: U in M.subsets()
            True

        The method :meth:`subset_family` returns a family (sorted
        alphabetically by the subset names)::

            sage: M.subset_family()
            Set {M, U, V} of subsets of the 2-dimensional topological manifold M
        """
        yield from self._subsets

    def list_of_subsets(self):
        r"""
        Return the list of subsets that have been defined on the current
        subset.

        The list is sorted by the alphabetical names of the subsets.

        OUTPUT:

        - a list containing all the subsets that have been defined on
          the current subset

        .. NOTE::

            This method is deprecated.

            To get the subsets as a :class:`ManifoldSubsetFiniteFamily`
            instance (which sorts its elements alphabetically by name),
            use :meth:`subset_family` instead.

            To loop over the subsets in an arbitrary order, use the
            generator method :meth:`subsets` instead.

        EXAMPLES:

        List of subsets of a 2-dimensional manifold (deprecated)::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: U = M.open_subset('U')
            sage: V = M.subset('V')
            sage: M.list_of_subsets()
            doctest:...: DeprecationWarning: the method list_of_subsets of ManifoldSubset
             is deprecated; use subset_family or subsets instead...
            [2-dimensional topological manifold M,
             Open subset U of the 2-dimensional topological manifold M,
             Subset V of the 2-dimensional topological manifold M]

        Using :meth:`subset_family` instead (recommended when order matters)::

            sage: M.subset_family()
            Set {M, U, V} of subsets of the 2-dimensional topological manifold M

        The method :meth:`subsets` generates the subsets in an unspecified order.
        To create a set::

            sage: frozenset(M.subsets())  # random (set output)
            {Subset V of the 2-dimensional topological manifold M,
             2-dimensional topological manifold M,
             Open subset U of the 2-dimensional topological manifold M}
        """
        deprecation(31727, "the method list_of_subsets of ManifoldSubset is deprecated; use subset_family or subsets instead")
        return sorted(self._subsets, key=lambda x: x._name)

    def subset_family(self):
        r"""
        Return the family of subsets that have been defined on the current subset.

        The family is sorted by the alphabetical names of the subsets.

        OUTPUT:

        - a :class:`ManifoldSubsetFiniteFamily` instance containing all the
          subsets that have been defined on the current subset

        .. NOTE::

            If you only need to iterate over the subsets in arbitrary order,
            you can use the generator method :meth:`subsets` instead.

        EXAMPLES:

        Subsets of a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: U = M.open_subset('U')
            sage: V = M.subset('V')
            sage: M.subset_family()
            Set {M, U, V} of subsets of the 2-dimensional topological manifold M
        """
        return ManifoldSubsetFiniteFamily(self.subsets())

    def subset_digraph(self, loops=False, quotient=False, open_covers=False, points=False, lower_bound=None):
        r"""
        Return the digraph whose arcs represent subset relations among the subsets of ``self``.

        INPUT:

        - ``loops`` -- boolean (default: ``False``); whether to include the trivial containment
          of each subset in itself as loops of the digraph
        - ``quotient`` -- boolean (default: ``False``); whether to contract directed
          cycles in the graph, replacing equivalence classes of equal subsets by a
          single vertex. In this case, each vertex of the digraph is a set of
          :class:`ManifoldSubset` instances.
        - ``open_covers`` -- boolean (default: ``False``); whether to include vertices for open covers
        - ``points`` -- boolean (default: ``False``); whether to include vertices for declared points;
          this can also be an iterable for the points to include
        - ``lower_bound`` -- (default: ``None``) only include supersets of this

        OUTPUT: a digraph; each vertex of the digraph is either:

        - a :class:`ManifoldSubsetFiniteFamily` containing one instance of :class:`ManifoldSubset`.
        - (if ``open_covers`` is ``True``) a tuple of :class:`ManifoldSubsetFiniteFamily` instances,
          representing an open cover.

        EXAMPLES::

            sage: # needs sage.graphs
            sage: M = Manifold(3, 'M')
            sage: U = M.open_subset('U'); V = M.open_subset('V'); W = M.open_subset('W')
            sage: D = M.subset_digraph(); D
            Digraph on 4 vertices
            sage: D.edges(sort=True, key=lambda e: (e[0]._name, e[1]._name))            # needs sage.graphs
            [(Set {U} of open subsets of the 3-dimensional differentiable manifold M,
              Set {M} of open subsets of the 3-dimensional differentiable manifold M,
              None),
             (Set {V} of open subsets of the 3-dimensional differentiable manifold M,
              Set {M} of open subsets of the 3-dimensional differentiable manifold M,
              None),
             (Set {W} of open subsets of the 3-dimensional differentiable manifold M,
              Set {M} of open subsets of the 3-dimensional differentiable manifold M,
              None)]
            sage: D.plot(layout='acyclic')                                              # needs sage.plot
            Graphics object consisting of 8 graphics primitives
            sage: def label(element):
            ....:     try:
            ....:         return element._name
            ....:     except AttributeError:
            ....:         return '[' + ', '.join(sorted(x._name for x in element)) + ']'
            sage: D.relabel(label, inplace=False).plot(layout='acyclic')                # needs sage.plot
            Graphics object consisting of 8 graphics primitives
            sage: VW = V.union(W)
            sage: D = M.subset_digraph(); D
            Digraph on 5 vertices
            sage: D.relabel(label, inplace=False).plot(layout='acyclic')                # needs sage.plot
            Graphics object consisting of 12 graphics primitives

        If ``open_covers`` is ``True``, the digraph includes a special vertex for
        each nontrivial open cover of a subset::

            sage: D = M.subset_digraph(open_covers=True)                                # needs sage.graphs
            sage: D.relabel(label, inplace=False).plot(layout='acyclic')                # needs sage.graphs sage.plot
            Graphics object consisting of 14 graphics primitives

        .. PLOT::

            def label(element):
                try:
                    return element._name
                except AttributeError:
                    return '[' + ', '.join(sorted(x._name for x in element)) + ']'
            M = Manifold(3, 'M')
            U = M.open_subset('U'); V = M.open_subset('V'); W = M.open_subset('W')
            D = M.subset_digraph()
            g1 = D.relabel(label, inplace=False).plot(layout='acyclic')
            VW = V.union(W)
            D = M.subset_digraph()
            g2 = D.relabel(label, inplace=False).plot(layout='acyclic')
            D = M.subset_digraph(open_covers=True)
            g3 = D.relabel(label, inplace=False).plot(layout='acyclic')
            sphinx_plot(graphics_array([g1, g2, g3]), figsize=(8, 3))
        """
        from sage.graphs.digraph import DiGraph
        D = DiGraph(multiedges=False, loops=loops)

        if loops:
            add_edges = D.add_edges
        else:
            def add_edges(edges):
                for u, v in edges:
                    if u != v:
                        D.add_edge((u, v))

        if quotient:
            def vertex_family(subset):
                return ManifoldSubsetFiniteFamily(subset.equal_subsets())
        else:
            def vertex_family(subset):
                return ManifoldSubsetFiniteFamily([subset])
        subset_to_vertex = {}

        def vertex(subset):
            try:
                return subset_to_vertex[subset]
            except KeyError:
                family = vertex_family(subset)
                for S in family:
                    subset_to_vertex[S] = family
                return family

        if lower_bound is not None:
            if not lower_bound.is_subset(self):
                return D
        visited = set()
        to_visit = [self]
        while to_visit:
            S = to_visit.pop()
            if S in visited:
                continue
            visited.add(S)

            if lower_bound is None:
                subsets = S._subsets
            else:
                subsets = [subset for subset in S._subsets
                           if lower_bound.is_subset(subset)]

            add_edges((vertex(subset), vertex(S)) for subset in subsets)

            subsets_without_S = [subset for subset in subsets
                                 if subset is not S]
            to_visit.extend(subsets_without_S)

        # Make sure to include isolated vertices in the graph
        D.add_vertices(subset_to_vertex.values())

        if open_covers:

            def open_cover_vertex(open_cover):
                return tuple(sorted(ManifoldSubsetFiniteFamily([subset]) for subset in open_cover))

            for S in visited:
                add_edges((vertex(S), open_cover_vertex(open_cover))
                          for open_cover in S.open_covers(trivial=False))

        if points is not False:
            subset_to_points = defaultdict(list)
            if points is not True:
                # Manifolds do not keep track of the points defined on them.
                # Use the provided iterator.
                def point_vertex(point):
                    return point

                for point in points:
                    S = point.parent()
                    subset_to_points[S].append(point)
                    D.add_edge((point_vertex(point), vertex(S)))

            # Add a placeholder vertex under each subset that has a defined
            # point that we do not know about.
            def anonymous_point_vertex(S):
                return f"p{S._name}"

            add_edges((anonymous_point_vertex(S), vertex(S))
                      for S in visited
                      if S.has_defined_points(subsets=False)
                      and S not in subset_to_points)

        return D

    def subset_poset(self, open_covers=False, points=False, lower_bound=None):
        r"""
        Return the poset of equivalence classes of the subsets of ``self``.

        Each element of the poset is a set of :class:`ManifoldSubset` instances,
        which are known to be equal.

        INPUT:

        - ``open_covers`` -- boolean (default: ``False``); whether to include vertices for open covers
        - ``points`` -- boolean (default: ``False``); whether to include vertices for declared points;
          this can also be an iterable for the points to include
        - ``lower_bound`` -- (default: ``None``) only include supersets of this

        EXAMPLES::

            sage: # needs sage.graphs
            sage: M = Manifold(3, 'M')
            sage: U = M.open_subset('U'); V = M.open_subset('V'); W = M.open_subset('W')
            sage: P = M.subset_poset(); P
            Finite poset containing 4 elements
            sage: P.plot(element_labels={element: element._name for element in P})      # needs sage.plot
            Graphics object consisting of 8 graphics primitives
            sage: VW = V.union(W)
            sage: P = M.subset_poset(); P
            Finite poset containing 5 elements
            sage: P.maximal_elements()
            [Set {M} of open subsets of the 3-dimensional differentiable manifold M]
            sage: sorted(P.minimal_elements(), key=lambda v: v._name)
             [Set {U} of open subsets of the 3-dimensional differentiable manifold M,
              Set {V} of open subsets of the 3-dimensional differentiable manifold M,
              Set {W} of open subsets of the 3-dimensional differentiable manifold M]
            sage: from sage.manifolds.subset import ManifoldSubsetFiniteFamily
            sage: sorted(P.lower_covers(ManifoldSubsetFiniteFamily([M])), key=str)
             [Set {U} of open subsets of the 3-dimensional differentiable manifold M,
              Set {V_union_W} of open subsets of the 3-dimensional differentiable manifold M]
            sage: P.plot(element_labels={element: element._name for element in P})      # needs sage.plot
            Graphics object consisting of 10 graphics primitives

        If ``open_covers`` is ``True``, the poset includes a special vertex for
        each nontrivial open cover of a subset::

            sage: # needs sage.graphs
            sage: P = M.subset_poset(open_covers=True); P
            Finite poset containing 6 elements
            sage: from sage.manifolds.subset import ManifoldSubsetFiniteFamily
            sage: sorted(P.upper_covers(ManifoldSubsetFiniteFamily([VW])), key=str)
            [(Set {V} of open subsets of the 3-dimensional differentiable manifold M,
              Set {W} of open subsets of the 3-dimensional differentiable manifold M),
             Set {M} of open subsets of the 3-dimensional differentiable manifold M]
            sage: def label(element):
            ....:     try:
            ....:         return element._name
            ....:     except AttributeError:
            ....:         return '[' + ', '.join(sorted(x._name for x in element)) + ']'
            sage: P.plot(element_labels={element: label(element) for element in P})     # needs sage.plot
            Graphics object consisting of 12 graphics primitives

        .. PLOT::

            def label(element):
                try:
                    return element._name
                except AttributeError:
                    return '[' + ', '.join(sorted(x._name for x in element)) + ']'
            M = Manifold(3, 'M')
            U = M.open_subset('U'); V = M.open_subset('V'); W = M.open_subset('W')
            P = M.subset_poset()
            g1 = P.plot(element_labels={element: label(element) for element in P})
            VW = V.union(W)
            P = M.subset_poset()
            g2 = P.plot(element_labels={element: label(element) for element in P})
            P = M.subset_poset(open_covers=True)
            g3 = P.plot(element_labels={element: label(element) for element in P})
            sphinx_plot(graphics_array([g1, g2, g3]), figsize=(8, 3))
        """
        from sage.combinat.posets.posets import Poset
        return Poset(self.subset_digraph(open_covers=open_covers, points=points,
                                         quotient=True, lower_bound=lower_bound))

    def equal_subsets(self):
        r"""
        Generate the declared manifold subsets that are equal to ``self``.

        .. NOTE::

            To get the equal subsets as a family, sorted by name, use the method
            :meth:`equal_subset_family` instead.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: U = M.open_subset('U')
            sage: V = U.subset('V')
            sage: V.declare_equal(M)
            sage: sorted(V.equal_subsets(), key=lambda v: v._name)
            [2-dimensional topological manifold M,
             Open subset U of the 2-dimensional topological manifold M,
             Subset V of the 2-dimensional topological manifold M]
        """
        for S in self.supersets():
            if S in self._subsets:
                yield S

    def equal_subset_family(self):
        r"""
        Generate the declared manifold subsets that are equal to ``self``.

        .. NOTE::

            If you only need to iterate over the equal sets in arbitrary order,
            you can use the generator method :meth:`equal_subsets` instead.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: U = M.open_subset('U')
            sage: V = U.subset('V')
            sage: V.declare_equal(M)
            sage: V.equal_subset_family()
            Set {M, U, V} of subsets of the 2-dimensional topological manifold M
        """
        return ManifoldSubsetFiniteFamily(self.equal_subsets())

    def supersets(self):
        r"""
        Generate the declared supersets of the current subset.

        .. NOTE::

            To get the supersets as a family, sorted by name, use the method
            :meth:`superset_family` instead.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: U = M.open_subset('U')
            sage: V = M.subset('V')
            sage: sorted(V.supersets(), key=lambda v: v._name)
            [2-dimensional topological manifold M,
             Subset V of the 2-dimensional topological manifold M]
        """
        yield from self._supersets

    def superset_family(self):
        r"""
        Return the family of declared supersets of the current subset.

        The family is sorted by the alphabetical names of the supersets.

        OUTPUT:

        - a :class:`ManifoldSubsetFiniteFamily` instance containing all the
          supersets

        .. NOTE::

            If you only need to iterate over the supersets in arbitrary order,
            you can use the generator method :meth:`supersets` instead.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: U = M.open_subset('U')
            sage: V = M.subset('V')
            sage: V.superset_family()
            Set {M, V} of subsets of the 2-dimensional topological manifold M
        """
        return ManifoldSubsetFiniteFamily(self.supersets())

    def superset_digraph(self, loops=False, quotient=False, open_covers=False, points=False, upper_bound=None):
        r"""
        Return the digraph whose arcs represent subset relations among the supersets of ``self``.

        INPUT:

        - ``loops`` -- boolean (default: ``False``); whether to include the trivial containment
          of each subset in itself as loops of the digraph
        - ``quotient`` -- boolean (default: ``False``); whether to contract
          directed cycles in the graph, replacing equivalence classes of equal
          subsets by a single vertex. In this case, each vertex of the digraph
          is a set of :class:`ManifoldSubset` instances.
        - ``open_covers`` -- boolean (default: ``False``); whether to include vertices for open covers
        - ``points`` -- boolean (default: ``False``); whether to include vertices for declared points;
          this can also be an iterable for the points to include
        - ``upper_bound`` -- (default: ``None``) only include subsets of this

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: U = M.open_subset('U'); V = M.open_subset('V'); W = M.open_subset('W')
            sage: VW = V.union(W)
            sage: P = V.superset_digraph(loops=False, upper_bound=VW); P                # needs sage.graphs
            Digraph on 2 vertices
        """
        if upper_bound is None:
            upper_bound = self._manifold
        return upper_bound.subset_digraph(loops=loops, open_covers=open_covers, points=points,
                                          quotient=quotient, lower_bound=self)

    def superset_poset(self, open_covers=False, points=False, upper_bound=None):
        r"""
        Return the poset of the supersets of ``self``.

        INPUT:

        - ``open_covers`` -- boolean (default: ``False``); whether to include vertices for open covers
        - ``points`` -- boolean (default: ``False``); whether to include vertices for declared points;
          this can also be an iterable for the points to include
        - ``upper_bound`` -- (default: ``None``) only include subsets of this

        EXAMPLES::

            sage: M = Manifold(3, 'M')
            sage: U = M.open_subset('U'); V = M.open_subset('V'); W = M.open_subset('W')
            sage: VW = V.union(W)
            sage: P = V.superset_poset(); P                                             # needs sage.graphs
            Finite poset containing 3 elements
            sage: P.plot(element_labels={element: element._name for element in P})      # needs sage.graphs sage.plot
            Graphics object consisting of 6 graphics primitives
        """
        if upper_bound is None:
            upper_bound = self._manifold
        return upper_bound.subset_poset(open_covers=open_covers, points=points,
                                        lower_bound=self)

    def get_subset(self, name):
        r"""
        Get a subset by its name.

        The subset must have been previously created by the method
        :meth:`subset` (or
        :meth:`~sage.manifolds.manifold.TopologicalManifold.open_subset`)

        INPUT:

        - ``name`` -- string; name of the subset

        OUTPUT:

        - instance of :class:`~sage.manifolds.subset.ManifoldSubset` (or
          of the derived class
          :class:`~sage.manifolds.manifold.TopologicalManifold` for an open
          subset) representing the subset whose name is ``name``

        EXAMPLES::

            sage: M = Manifold(4, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: B = A.subset('B')
            sage: U = M.open_subset('U')
            sage: M.subset_family()
            Set {A, B, M, U} of subsets of the 4-dimensional topological manifold M
            sage: M.get_subset('A')
            Subset A of the 4-dimensional topological manifold M
            sage: M.get_subset('A') is A
            True
            sage: M.get_subset('B') is B
            True
            sage: A.get_subset('B') is B
            True
            sage: M.get_subset('U')
            Open subset U of the 4-dimensional topological manifold M
            sage: M.get_subset('U') is U
            True
        """
        for ss in self._subsets:
            if ss._name == name:
                return ss
        raise ValueError("no subset of name '{}' found".format(name))

    #### End of accessors

    def is_subset(self, other):
        r"""
        Return ``True`` if and only if ``self`` is included in ``other``.

        EXAMPLES:

        Subsets on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: a = M.subset('A')
            sage: b = a.subset('B')
            sage: c = M.subset('C')
            sage: a.is_subset(M)
            True
            sage: b.is_subset(a)
            True
            sage: b.is_subset(M)
            True
            sage: a.is_subset(b)
            False
            sage: c.is_subset(a)
            False
        """
        return self in other._subsets

    def declare_union(self, *subsets_or_families, disjoint=False):
        r"""
        Declare that the current subset is the union of two subsets.

        Suppose `U` is the current subset, then this method declares
        that `U = \bigcup_{S\in F} S`.

        INPUT:

        - ``subsets_or_families`` -- finitely many subsets or iterables of subsets
        - ``disjoint`` -- boolean (default: ``False``); whether to declare the subsets
          pairwise disjoint

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: AB = M.subset('AB')
            sage: A = AB.subset('A')
            sage: B = AB.subset('B')
            sage: def label(element):
            ....:     try:
            ....:         return element._name
            ....:     except AttributeError:
            ....:         return '[' + ', '.join(sorted(x._name for x in element)) + ']'
            sage: P = M.subset_poset(open_covers=True); P                               # needs sage.graphs
            Finite poset containing 4 elements
            sage: P.plot(element_labels={element: label(element) for element in P})     # needs sage.graphs sage.plot
            Graphics object consisting of 8 graphics primitives

            sage: AB.declare_union(A, B)
            sage: A.union(B)
            Subset AB of the 2-dimensional topological manifold M
            sage: P = M.subset_poset(open_covers=True); P                               # needs sage.graphs
            Finite poset containing 4 elements
            sage: P.plot(element_labels={element: label(element) for element in P})     # needs sage.graphs sage.plot
            Graphics object consisting of 8 graphics primitives

            sage: B1 = B.subset('B1', is_open=True)
            sage: B2 = B.subset('B2', is_open=True)
            sage: B.declare_union(B1, B2, disjoint=True)
            sage: P = M.subset_poset(open_covers=True); P                               # needs sage.graphs
            Finite poset containing 9 elements
            sage: P.plot(element_labels={element: label(element) for element in P})     # needs sage.graphs sage.plot
            Graphics object consisting of 19 graphics primitives

        .. PLOT::

            def label(element):
                try:
                    return element._name
                except AttributeError:
                    return '[' + ', '.join(sorted(x._name for x in element)) + ']'
            M = Manifold(2, 'M', structure='topological')
            AB = M.subset('AB')
            A = AB.subset('A')
            B = AB.subset('B')
            P = M.subset_poset(open_covers=True); P
            g1 = P.plot(element_labels={element: label(element) for element in P})
            AB.declare_union(A, B)
            A.union(B)
            P = M.subset_poset(open_covers=True); P
            g2 = P.plot(element_labels={element: label(element) for element in P})
            B1 = B.subset('B1', is_open=True)
            B2 = B.subset('B2', is_open=True)
            B.declare_union(B1, B2, disjoint=True)
            P = M.subset_poset(open_covers=True); P
            g3 = P.plot(element_labels={element: label(element) for element in P})
            sphinx_plot(graphics_array([g1, g2, g3]), figsize=(8, 3))
        """
        subsets = ManifoldSubsetFiniteFamily.from_subsets_or_families(*subsets_or_families)
        if disjoint:
            for U, V in itertools.combinations(subsets, 2):
                U.intersection(V).declare_empty()
        subsets = self._reduce_union_members(subsets)
        if not subsets:
            self.declare_empty()
        elif len(subsets) == 1:
            self.declare_equal(*subsets)
        else:
            subset_iter = iter(subsets)
            first = next(subset_iter)
            second = next(subset_iter)
            self._declare_union_2_subsets(first, second.union(subset_iter))

    def _declare_union_2_subsets(self, dom1, dom2):
        r"""
        Declare that the current subset is the union of two of its subsets.

        Suppose `U` is the current subset, then this method declares that

        .. MATH::

            U = U_1 \cup U_2,

        where `U_1 \subset U` and `U_2 \subset U`.

        INPUT:

        - ``dom1`` -- the subset `U_1`
        - ``dom2`` -- the subset `U_2`

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: B = M.subset('B')
            sage: M.declare_union(A, B)
            sage: A.union(B)
            2-dimensional topological manifold M
        """
        if dom1 == dom2:
            if dom1 != self:
                raise ValueError("the union of two identical sets must be " +
                                 "this set")
            return
        if not dom1.is_subset(self):
            raise TypeError("the {} is not a subset of ".format(dom1) +
                            "the {}".format(self))
        if not dom2.is_subset(self):
            raise TypeError("the {} is not a subset of ".format(dom2) +
                            "the {}".format(self))
        dom1._unions[dom2._name] = self
        dom2._unions[dom1._name] = self
        for oc1 in dom1._open_covers:
            for oc2 in dom2._open_covers:
                oc = oc1[:]
                for s in oc2:
                    if s not in oc:
                        oc.append(s)
                self._open_covers.append(oc)

    def declare_equal(self, *others):
        r"""
        Declare that ``self`` and ``others`` are the same sets.

        INPUT:

        - ``others`` -- finitely many subsets or iterables of subsets of the same
          manifold as ``self``

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U')
            sage: V = M.open_subset('V')
            sage: Vs = [M.open_subset(f'V{i}') for i in range(2)]
            sage: UV = U.intersection(V)
            sage: W = UV.open_subset('W')
            sage: P = M.subset_poset()                                                  # needs sage.graphs
            sage: def label(element):
            ....:     return element._name
            sage: P.plot(element_labels={element: label(element) for element in P})     # needs sage.graphs sage.plot
            Graphics object consisting of 15 graphics primitives
            sage: V.declare_equal(Vs)
            sage: P = M.subset_poset()                                                  # needs sage.graphs
            sage: P.plot(element_labels={element: label(element) for element in P})     # needs sage.graphs sage.plot
            Graphics object consisting of 11 graphics primitives
            sage: W.declare_equal(U)
            sage: P = M.subset_poset()                                                  # needs sage.graphs
            sage: P.plot(element_labels={element: label(element) for element in P})     # needs sage.graphs sage.plot
            Graphics object consisting of 6 graphics primitives

        .. PLOT::

            def label(element):
                return element._name
            M = Manifold(2, 'M')
            U = M.open_subset('U')
            V = M.open_subset('V')
            Vs = [M.open_subset(f'V{i}') for i in range(2)]
            UV = U.intersection(V)
            W = UV.open_subset('W')
            P = M.subset_poset()
            g1 = P.plot(element_labels={element: label(element) for element in P})
            V.declare_equal(Vs)
            P = M.subset_poset()
            g2 = P.plot(element_labels={element: label(element) for element in P})
            W.declare_equal(U)
            P = M.subset_poset()
            g3 = P.plot(element_labels={element: label(element) for element in P})
            sphinx_plot(graphics_array([g1, g2, g3]), figsize=(8, 3))
        """
        F = ManifoldSubsetFiniteFamily.from_subsets_or_families
        equal_sets = F(self, *others)
        all_supersets = F(*[S.supersets() for S in equal_sets])
        all_subsets = F(*[S.subsets() for S in equal_sets])
        for superset in all_supersets:
            superset._subsets.update(all_subsets)
        for subset in all_subsets:
            subset._supersets.update(all_supersets)

    def declare_subset(self, *supersets):
        r"""
        Declare ``self`` to be a subset of each of the given supersets.

        INPUT:

        - ``supersets`` -- other subsets of the same manifold

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: U1 = M.open_subset('U1')
            sage: U2 = M.open_subset('U2')
            sage: V = M.open_subset('V')
            sage: V.superset_family()
            Set {M, V} of open subsets of the 2-dimensional differentiable manifold M
            sage: U1.subset_family()
            Set {U1} of open subsets of the 2-dimensional differentiable manifold M
            sage: P = M.subset_poset()                                                  # needs sage.graphs
            sage: def label(element):
            ....:     return element._name
            sage: P.plot(element_labels={element: label(element) for element in P})     # needs sage.graphs sage.plot
            Graphics object consisting of 8 graphics primitives
            sage: V.declare_subset(U1, U2)
            sage: V.superset_family()
            Set {M, U1, U2, V} of open subsets of the 2-dimensional differentiable manifold M
            sage: P = M.subset_poset()                                                  # needs sage.graphs
            sage: P.plot(element_labels={element: label(element) for element in P})     # needs sage.graphs sage.plot
            Graphics object consisting of 9 graphics primitives

        Subsets in a directed cycle of inclusions are equal::

            sage: M.declare_subset(V)
            sage: M.superset_family()
            Set {M, U1, U2, V} of open subsets of the 2-dimensional differentiable manifold M
            sage: M.equal_subset_family()
            Set {M, U1, U2, V} of open subsets of the 2-dimensional differentiable manifold M
            sage: P = M.subset_poset()                                                  # needs sage.graphs
            sage: P.plot(element_labels={element: label(element) for element in P})     # needs sage.graphs sage.plot
            Graphics object consisting of 2 graphics primitives

        .. PLOT::

            def label(element):
                return element._name
            M = Manifold(2, 'M')
            U1 = M.open_subset('U1')
            U2 = M.open_subset('U2')
            V = M.open_subset('V')
            P = M.subset_poset()
            g1 = P.plot(element_labels={element: label(element) for element in P})
            V.declare_subset(U1, U2)
            P = M.subset_poset()
            g2 = P.plot(element_labels={element: label(element) for element in P})
            M.declare_subset(V)
            P = M.subset_poset()
            g3 = P.plot(element_labels={element: label(element) for element in P})
            sphinx_plot(graphics_array([g1, g2, g3]), figsize=(8, 3))
        """
        F = ManifoldSubsetFiniteFamily.from_subsets_or_families
        supersets = F(*supersets)
        all_supersets = F(*[S.supersets() for S in supersets])
        for superset in all_supersets:
            superset._subsets.update(self._subsets)
        for subset in self._subsets:
            subset._supersets.update(all_supersets)

    def declare_superset(self, *subsets):
        r"""
        Declare ``self`` to be a superset of each of the given subsets.

        INPUT:

        - ``subsets`` -- other subsets of the same manifold

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U')
            sage: V1 = M.open_subset('V1')
            sage: V2 = M.open_subset('V2')
            sage: W = V1.intersection(V2)
            sage: U.subset_family()
            Set {U} of open subsets of the 2-dimensional differentiable manifold M
            sage: P = M.subset_poset()                                                  # needs sage.graphs
            sage: def label(element):
            ....:     return element._name
            sage: P.plot(element_labels={element: label(element) for element in P})     # needs sage.graphs sage.plot
            Graphics object consisting of 11 graphics primitives
            sage: U.declare_superset(V1, V2)
            sage: U.subset_family()
            Set {U, V1, V1_inter_V2, V2} of open subsets of the 2-dimensional differentiable manifold M
            sage: P = M.subset_poset()                                                  # needs sage.graphs
            sage: P.plot(element_labels={element: label(element) for element in P})     # needs sage.graphs sage.plot
            Graphics object consisting of 11 graphics primitives

        Subsets in a directed cycle of inclusions are equal::

            sage: W.declare_superset(U)
            sage: W.subset_family()
            Set {U, V1, V1_inter_V2, V2} of open subsets of the 2-dimensional differentiable manifold M
            sage: W.equal_subset_family()
            Set {U, V1, V1_inter_V2, V2} of open subsets of the 2-dimensional differentiable manifold M
            sage: P = M.subset_poset()                                                  # needs sage.graphs
            sage: P.plot(element_labels={element: label(element) for element in P})     # needs sage.graphs sage.plot
            Graphics object consisting of 4 graphics primitives

        .. PLOT::

            def label(element):
                return element._name
            M = Manifold(2, 'M')
            U = M.open_subset('U')
            V1 = M.open_subset('V1')
            V2 = M.open_subset('V2')
            W = V1.intersection(V2)
            P = M.subset_poset()
            def label(element):
                return element._name
            g1 = P.plot(element_labels={element: label(element) for element in P})
            U.declare_superset(V1, V2)
            P = M.subset_poset()
            g2 = P.plot(element_labels={element: label(element) for element in P})
            W.declare_superset(U)
            P = M.subset_poset()
            g3 = P.plot(element_labels={element: label(element) for element in P})
            sphinx_plot(graphics_array([g1, g2, g3]), figsize=(8, 3))
        """
        F = ManifoldSubsetFiniteFamily.from_subsets_or_families
        subsets = F(*subsets)
        all_subsets = F(*[S.subsets() for S in subsets])
        for subset in all_subsets:
            subset._supersets.update(self._supersets)
        for superset in self._supersets:
            superset._subsets.update(all_subsets)

    def declare_empty(self):
        r"""
        Declare that ``self`` is the empty set.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A', is_open=True)
            sage: AA = A.subset('AA')
            sage: A
            Open subset A of the 2-dimensional topological manifold M
            sage: A.declare_empty()
            sage: A.is_empty()
            True

        Empty sets do not allow to define points on them::

            sage: A.point()
            Traceback (most recent call last):
            ...
            TypeError: cannot define a point on the
              Open subset A of the 2-dimensional topological manifold M
              because it has been declared empty

        Emptiness transfers to subsets::

            sage: AA.is_empty()
            True
            sage: AA.point()
            Traceback (most recent call last):
            ...
            TypeError: cannot define a point on the
              Subset AA of the 2-dimensional topological manifold M
              because it has been declared empty
            sage: AD = A.subset('AD')
            sage: AD.is_empty()
            True

        If points have already been defined on ``self`` (or its subsets),
        it is an error to declare it to be empty::

            sage: B = M.subset('B')
            sage: b = B.point(name='b'); b
            Point b on the 2-dimensional topological manifold M
            sage: B.declare_empty()
            Traceback (most recent call last):
            ...
            TypeError: cannot be empty because it has defined points

        Emptiness is recorded as empty open covers::

            sage: P = M.subset_poset(open_covers=True, points=[b])                      # needs sage.graphs
            sage: def label(element):
            ....:     if isinstance(element, str):
            ....:         return element
            ....:     try:
            ....:         return element._name
            ....:     except AttributeError:
            ....:         return '[' + ', '.join(sorted(x._name for x in element)) + ']'
            sage: P.plot(element_labels={element: label(element) for element in P})     # needs sage.graphs sage.plot
            Graphics object consisting of 10 graphics primitives

        .. PLOT::

            def label(element):
                if isinstance(element, str):
                    return element
                try:
                    return element._name
                except AttributeError:
                    return '[' + ', '.join(sorted(x._name for x in element)) + ']'
            M = Manifold(2, 'M', structure='topological')
            A = M.subset('A', is_open=True)
            AA = A.subset('AA')
            A.declare_empty()
            AD = A.subset('AD')
            B = M.subset('B')
            b = B.point(name='b')

            D = M.subset_digraph(open_covers=True, points=[b])
            g1 = D.relabel(label, inplace=False).plot(layout='spring')
            P = M.subset_poset(open_covers=True, points=[b])
            g2 = P.plot(element_labels={element: label(element) for element in P})
            sphinx_plot(graphics_array([g1, g2]), figsize=(8, 5))
        """
        if self.has_defined_points():
            raise TypeError('cannot be empty because it has defined points')
        if not self.is_empty():
            self._open_covers.append([])
            self.declare_equal(self.subsets())

    def is_empty(self):
        r"""
        Return whether the current subset is empty.

        By default, manifold subsets are considered nonempty: The method :meth:`point` can be
        used to define points on it, either with or without coordinates some chart.

        However, using :meth:`declare_empty`, a subset can be declared empty, and emptiness
        transfers to all of its subsets.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A', is_open=True)
            sage: AA = A.subset('AA')
            sage: A.is_empty()
            False
            sage: A.declare_empty()
            sage: A.is_empty()
            True
            sage: AA.is_empty()
            True
        """
        if self.has_defined_points(subsets=False):
            # Fast path, do not check subsets
            return False
        return any(not cover
                   for cover in self.open_covers(trivial=False, supersets=True))

    def declare_nonempty(self):
        r"""
        Declare that ``self`` is nonempty.

        Once declared nonempty, ``self`` (or any of its supersets) cannot be declared empty.

        This is equivalent to defining a point on ``self`` using :meth:`point`
        but is cheaper than actually creating a :class:`~sage.manifolds.point.ManifoldPoint`
        instance.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A', is_open=True)
            sage: AA = A.subset('AA')
            sage: AA.declare_nonempty()
            sage: A.has_defined_points()
            True
            sage: A.declare_empty()
            Traceback (most recent call last):
            ...
            TypeError: cannot be empty because it has defined points
        """
        if self.has_defined_points(subsets=False):
            # Fast path, do not check subsets
            return
        if self.is_empty():
            raise TypeError('cannot be nonempty because it has already been declared empty')
        self._has_defined_points = True

    def has_defined_points(self, subsets=True):
        r"""
        Return whether any points have been defined on ``self`` or any of its subsets.

        INPUT:

        - ``subsets`` -- boolean (default: ``True``); if ``False``, only consider points that have
          been defined directly on ``self``. If ``True``, also consider points on all subsets.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A', is_open=True)
            sage: AA = A.subset('AA')
            sage: AA.point()
            Point on the 2-dimensional topological manifold M
            sage: AA.has_defined_points()
            True
            sage: A.has_defined_points(subsets=False)
            False
            sage: A.has_defined_points()
            True
        """
        if subsets:
            return any(subset._has_defined_points for subset in self.subsets())
        else:
            return self._has_defined_points

    def point(self, coords=None, chart=None, name=None, latex_name=None):
        r"""
        Define a point in ``self``.

        See :class:`~sage.manifolds.point.ManifoldPoint` for a
        complete documentation.

        INPUT:

        - ``coords`` -- the point coordinates (as a tuple or a list) in the
          chart specified by ``chart``
        - ``chart`` -- (default: ``None``) chart in which the point
          coordinates are given; if ``None``, the coordinates are assumed
          to refer to the default chart of the current subset
        - ``name`` -- (default: ``None``) name given to the point
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          point; if ``None``, the LaTeX symbol is set to ``name``

        OUTPUT:

        - the declared point, as an instance of
          :class:`~sage.manifolds.point.ManifoldPoint`

        EXAMPLES:

        Points on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: c_xy.<x,y> = M.chart()
            sage: p = M.point((1,2), name='p'); p
            Point p on the 2-dimensional topological manifold M
            sage: p in M
            True
            sage: a = M.open_subset('A')
            sage: c_uv.<u,v> = a.chart()
            sage: q = a.point((-1,0), name='q'); q
            Point q on the 2-dimensional topological manifold M
            sage: q in a
            True
            sage: p._coordinates
            {Chart (M, (x, y)): (1, 2)}
            sage: q._coordinates
            {Chart (A, (u, v)): (-1, 0)}
        """
        return self.element_class(self, coords=coords, chart=chart,
                                  name=name, latex_name=latex_name)

    def declare_closed(self):
        r"""
        Declare ``self`` to be a closed subset of the manifold.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: B1 = A.subset('B1')
            sage: B1.is_closed()
            False
            sage: B1.declare_closed()
            sage: B1.is_closed()
            True

            sage: B2 = A.subset('B2')
            sage: cl_B2 = B2.closure()
            sage: A.declare_closed()
            sage: cl_B2.is_subset(A)
            True
        """
        if self.is_closed():
            return
        self.complement(is_open=True)
        from .subsets.closure import ManifoldSubsetClosure
        for closure in self.manifold().subsets():
            if isinstance(closure, ManifoldSubsetClosure):
                if closure._subset.is_subset(self):
                    closure.declare_subset(self)

    #### Construction of new sets from self:

    def subset(self, name, latex_name=None, is_open=False):
        r"""
        Create a subset of the current subset.

        INPUT:

        - ``name`` -- name given to the subset
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote
          the subset; if none are provided, it is set to ``name``
        - ``is_open`` -- boolean (default: ``False``); if ``True``, the created subset
          is assumed to be open with respect to the manifold's topology

        OUTPUT:

        - the subset, as an instance of :class:`ManifoldSubset`, or
          of the derived class
          :class:`~sage.manifolds.manifold.TopologicalManifold`
          if ``is_open`` is ``True``

        EXAMPLES:

        Creating a subset of a manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: a = M.subset('A'); a
            Subset A of the 2-dimensional topological manifold M

        Creating a subset of ``A``::

            sage: b = a.subset('B', latex_name=r'\mathcal{B}'); b
            Subset B of the 2-dimensional topological manifold M
            sage: latex(b)
            \mathcal{B}

        We have then::

            sage: b.is_subset(a)
            True
            sage: b in a.subsets()
            True
        """
        if is_open:
            return self.open_subset(name, latex_name=latex_name)
        res = ManifoldSubset(self._manifold, name, latex_name=latex_name)
        if self.is_empty():
            self.declare_equal(res)
        else:
            self.declare_superset(res)
        self._top_subsets.add(res)
        return res

    def open_subset(self, name, latex_name=None, coord_def={}, supersets=None):
        r"""
        Create an open subset of the manifold that is a subset of ``self``.

        An open subset is a set that is (i) included in the manifold and (ii)
        open with respect to the manifold's topology. It is a topological
        manifold by itself. Hence the returned object is an instance of
        :class:`~sage.manifolds.manifold.TopologicalManifold`.

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

        OUTPUT:

        - the open subset, as an instance of
          :class:`~sage.manifolds.manifold.TopologicalManifold`
          or one of its subclasses

        EXAMPLES::

            sage: M = Manifold(2, 'R^2', structure='topological')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: cl_D = M.subset('cl_D'); cl_D
            Subset cl_D of the 2-dimensional topological manifold R^2
            sage: D = cl_D.open_subset('D', coord_def={c_cart: x^2+y^2<1}); D
            Open subset D of the 2-dimensional topological manifold R^2
            sage: D.is_subset(cl_D)
            True
            sage: D.is_subset(M)
            True

            sage: M = Manifold(2, 'R^2', structure='differentiable')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: cl_D = M.subset('cl_D'); cl_D
            Subset cl_D of the 2-dimensional differentiable manifold R^2
            sage: D = cl_D.open_subset('D', coord_def={c_cart: x^2+y^2<1}); D
            Open subset D of the 2-dimensional differentiable manifold R^2
            sage: D.is_subset(cl_D)
            True
            sage: D.is_subset(M)
            True

            sage: M = Manifold(2, 'R^2', structure='Riemannian')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: cl_D = M.subset('cl_D'); cl_D
            Subset cl_D of the 2-dimensional Riemannian manifold R^2
            sage: D = cl_D.open_subset('D', coord_def={c_cart: x^2+y^2<1}); D
            Open subset D of the 2-dimensional Riemannian manifold R^2
            sage: D.is_subset(cl_D)
            True
            sage: D.is_subset(M)
            True
        """
        if supersets is None:
            supersets = set()
        else:
            supersets = set(supersets)
        supersets.update([self])
        # Delegate to the manifold's method.
        return self._manifold.open_subset(name, latex_name=latex_name,
                                          coord_def=coord_def,
                                          supersets=supersets)

    def _init_open_subset(self, resu, coord_def):
        r"""
        Initialize ``resu`` as an open subset of ``self``.

        INPUT:

        - ``resu`` -- an instance of :class:`TopologicalManifold` or
          a subclass

        - ``coord_def`` -- (default: {}) definition of the subset in
          terms of coordinates; ``coord_def`` must a be dictionary with keys
          charts on the manifold and values the symbolic expressions formed
          by the coordinates to define the subset

        EXAMPLES::

            sage: M = Manifold(2, 'R^2', structure='topological')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: cl_D = M.subset('cl_D')
            sage: coord_def = {c_cart: x^2+y^2<1}
            sage: D = M.open_subset('D', coord_def=coord_def)
            sage: D.is_subset(cl_D)
            False
            sage: cl_D._init_open_subset(D, coord_def)
            sage: D.is_subset(cl_D)
            True
        """
        resu._supersets.update(self._supersets)
        self._subsets.add(resu)
        # Recursively delegate to the supersets.
        for superset in self._supersets:
            if superset is not self:
                superset._init_open_subset(resu, coord_def=coord_def)

    def superset(self, name, latex_name=None, is_open=False):
        r"""
        Create a superset of the current subset.

        A *superset* is a manifold subset in which the current subset is
        included.

        INPUT:

        - ``name`` -- name given to the superset
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote
          the superset; if none are provided, it is set to ``name``
        - ``is_open`` -- boolean (default: ``False``); if ``True``, the created
          subset is assumed to be open with respect to the manifold's topology

        OUTPUT:

        - the superset, as an instance of :class:`ManifoldSubset` or
          of the derived class
          :class:`~sage.manifolds.manifold.TopologicalManifold`
          if ``is_open`` is ``True``

        EXAMPLES:

        Creating some superset of a given subset::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: a = M.subset('A')
            sage: b = a.superset('B'); b
            Subset B of the 2-dimensional topological manifold M
            sage: b.subset_family()
            Set {A, B} of subsets of the 2-dimensional topological manifold M
            sage: a.superset_family()
            Set {A, B, M} of subsets of the 2-dimensional topological manifold M

        The superset of the whole manifold is itself::

            sage: M.superset('SM') is M
            True

        Two supersets of a given subset are a priori different::

            sage: c = a.superset('C')
            sage: c == b
            False
        """
        if self is self._manifold:
            return self
        if is_open:
            res = self._manifold.open_subset(name, latex_name=latex_name)
        else:
            res = ManifoldSubset(self._manifold, name, latex_name=latex_name)
        res.declare_superset(self)
        if is_open and self._is_open:
            res._atlas = list(self._atlas)
            res._top_charts = list(self._top_charts)
            res._coord_changes = dict(self._coord_changes)
            res._def_chart = self._def_chart
        return res

    def intersection(self, *others: ManifoldSubset, name: Optional[str] = None, latex_name: Optional[str] = None) -> ManifoldSubset:
        r"""
        Return the intersection of the current subset with other subsets.

        This method may return a previously constructed intersection instead
        of creating a new subset.  In this case, ``name`` and ``latex_name``
        are not used.

        INPUT:

        - ``others`` -- other subsets of the same manifold
        - ``name`` -- (default: ``None``) name given to the intersection
          in the case the latter has to be created; the default is
          ``self._name`` inter ``other._name``
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          intersection in the case the latter has to be created; the default
          is built upon the symbol `\cap`

        OUTPUT:

        - instance of :class:`ManifoldSubset` representing the
          subset that is the intersection of the current subset with ``others``

        EXAMPLES:

        Intersection of two subsets::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: a = M.subset('A')
            sage: b = M.subset('B')
            sage: c = a.intersection(b); c
            Subset A_inter_B of the 2-dimensional topological manifold M
            sage: a.subset_family()
            Set {A, A_inter_B} of subsets of the 2-dimensional topological manifold M
            sage: b.subset_family()
            Set {A_inter_B, B} of subsets of the 2-dimensional topological manifold M
            sage: c.superset_family()
            Set {A, A_inter_B, B, M} of subsets of the 2-dimensional topological manifold M

        Intersection of six subsets::

            sage: T = Manifold(2, 'T', structure='topological')
            sage: S = [T.subset(f'S{i}') for i in range(6)]
            sage: [S[i].intersection(S[i+3]) for i in range(3)]
            [Subset S0_inter_S3 of the 2-dimensional topological manifold T,
             Subset S1_inter_S4 of the 2-dimensional topological manifold T,
             Subset S2_inter_S5 of the 2-dimensional topological manifold T]
            sage: inter_S_i = T.intersection(*S, name='inter_S_i'); inter_S_i
            Subset inter_S_i of the 2-dimensional topological manifold T
            sage: inter_S_i.superset_family()
            Set {S0, S0_inter_S3, S0_inter_S3_inter_S1_inter_S4, S1, S1_inter_S4,
                 S2, S2_inter_S5, S3, S4, S5, T, inter_S_i} of
             subsets of the 2-dimensional topological manifold T

        .. PLOT::

            def label(element):
                if isinstance(element, str):
                    return element
                try:
                    return element._name.replace('_inter_', '∩')
                except AttributeError:
                    return '[' + ', '.join(sorted(label(x) for x in element)) + ']'

            M = Manifold(2, 'M', structure='topological')
            a = M.subset('A')
            b = M.subset('B')
            c = a.intersection(b); c
            P = M.subset_poset(open_covers=True)
            g1 = P.plot(element_labels={element: label(element) for element in P})

            T = Manifold(2, 'T', structure='topological')
            from sage.typeset.unicode_art import unicode_subscript
            S = [T.subset(f'S{unicode_subscript(i)}') for i in range(6)]
            [S[i].intersection(S[i+3]) for i in range(3)]
            T.intersection(*S, name='⋂ᵢSᵢ')
            P = T.subset_poset(open_covers=True)
            g2 = P.plot(element_labels={element: label(element) for element in P})

            sphinx_plot(graphics_array([g1, g2]), figsize=(8, 3))

        TESTS::

            sage: (a.intersection(b)).is_subset(a)
            True
            sage: (a.intersection(b)).is_subset(a)
            True
            sage: a.intersection(b) is b.intersection(a)
            True
            sage: a.intersection(a.intersection(b)) is a.intersection(b)
            True
            sage: (a.intersection(b)).intersection(a) is a.intersection(b)
            True
            sage: M.intersection(a) is a
            True
            sage: a.intersection(M) is a
            True
        """
        subsets = ManifoldSubsetFiniteFamily.from_subsets_or_families(self, *others)
        subset_iter = iter(self._reduce_intersection_members(subsets))
        # _intersection_subset is able to build the intersection of several
        # subsets directly; but because we cache only pairwise intersections,
        # we build the intersection by a sequence of pairwise intersections.
        res = next(subset_iter)
        others = list(subset_iter)
        if not others:
            return res
        for other in others[:-1]:
            res = res._intersection_subset(other)
        # The last one gets the name
        return res._intersection_subset(others[-1], name=name, latex_name=latex_name)

    @staticmethod
    def _reduce_intersection_members(subsets):
        r"""
        Return a reduced set of subsets with the same intersection as the given subsets.

        It is reduced with respect to two operations:

        - replacing an inclusion chain by its minimal element

        - replacing a pair of subsets with a declared intersection by the intersection

        INPUT:

        - ``subsets`` -- a non-empty iterable of :class:`ManifoldSubset` instances
          of the same manifold

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: B1 = A.subset('B1')
            sage: B2 = A.subset('B2')
            sage: C = B1.intersection(B2)
            sage: M._reduce_intersection_members([A, M, A])
            Set {A} of subsets of the 2-dimensional topological manifold M
            sage: M._reduce_intersection_members([A, B1])
            Set {B1} of subsets of the 2-dimensional topological manifold M
            sage: M._reduce_intersection_members([B1, B2])
            Set {B1_inter_B2} of subsets of the 2-dimensional topological manifold M
            sage: M._reduce_intersection_members([])
            Traceback (most recent call last):
            ...
            TypeError: input set must be nonempty
        """
        subsets = set(subsets)
        if not subsets:
            raise TypeError('input set must be nonempty')

        def reduce():
            # Greedily replace inclusion chains by their minimal element
            # and pairs with declared intersections by their intersection
            for U, V in itertools.combinations(subsets, 2):
                if U.is_subset(V):
                    subsets.remove(V)
                    return True
                if V.is_subset(U):
                    subsets.remove(U)
                    return True
                try:
                    UV = U._intersections[V._name]
                except KeyError:
                    pass
                else:
                    subsets.difference_update([U, V])
                    subsets.add(UV)
                    return True
            return False
        while reduce():
            pass
        assert subsets   # there must be a survivor
        return ManifoldSubsetFiniteFamily(subsets)

    def _intersection_subset(self, *others, name=None, latex_name=None):
        r"""
        Return a subset that is the intersection of ``self`` and ``others``.

        The result is always a new subset of the manifold.  If the intersection
        involves two subsets only, the result is stored in the dictionaries
        of known intersections for later reuse by other methods.

        INPUT:

        - ``others`` -- an iterable of :class:`ManifoldSubset` instances
          of the same manifold
        - ``name`` -- (default: ``None``) name given to the intersection; the
          default is ``self._name`` inter [...] inter ``last_other._name``
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          intersection; the default is built upon the symbol `\cap`

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: B1 = M.subset('B1')
            sage: B2 = M.subset('B2')
            sage: B3 = M.subset('B3')
            sage: B1._intersection_subset(B2)
            Subset B1_inter_B2 of the 2-dimensional topological manifold M
            sage: B1._intersection_subset(B2, B3)
            Subset B1_inter_B2_inter_B3 of the 2-dimensional topological manifold M
        """
        subsets = ManifoldSubsetFiniteFamily.from_subsets_or_families(self, *others)
        if latex_name is None:
            if name is None:
                latex_name = r'\cap '.join(S._latex_name for S in subsets)
            else:
                latex_name = name
        if name is None:
            name = "_inter_".join(S._name for S in subsets)
        if all(S.is_open() for S in subsets):
            res = self.open_subset(name, latex_name=latex_name, supersets=subsets)
        else:
            res = self.subset(name, latex_name=latex_name)
            res.declare_subset(subsets)
            for S in subsets:
                S._top_subsets.add(res)
        if len(subsets) == 2:
            S1, S2 = subsets
            S1._intersections[S2._name] = S2._intersections[S1._name] = res
        return res

    def union(self, *others, name=None, latex_name=None):
        r"""
        Return the union of the current subset with other subsets.

        This method may return a previously constructed union instead
        of creating a new subset.  In this case, ``name`` and ``latex_name``
        are not used.

        INPUT:

        - ``others`` -- other subsets of the same manifold
        - ``name`` -- (default: ``None``) name given to the union in the
          case the latter has to be created; the default is
          ``self._name`` union ``other._name``
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          union in the case the latter has to be created; the default
          is built upon the symbol `\cup`

        OUTPUT:

        - instance of :class:`ManifoldSubset` representing the
          subset that is the union of the current subset with ``others``

        EXAMPLES:

        Union of two subsets::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: a = M.subset('A')
            sage: b = M.subset('B')
            sage: c = a.union(b); c
            Subset A_union_B of the 2-dimensional topological manifold M
            sage: a.superset_family()
            Set {A, A_union_B, M} of subsets of the 2-dimensional topological manifold M
            sage: b.superset_family()
            Set {A_union_B, B, M} of subsets of the 2-dimensional topological manifold M
            sage: c.superset_family()
            Set {A_union_B, M} of subsets of the 2-dimensional topological manifold M

        Union of six subsets::

            sage: T = Manifold(2, 'T', structure='topological')
            sage: S = [T.subset(f'S{i}') for i in range(6)]
            sage: [S[i].union(S[i+3]) for i in range(3)]
            [Subset S0_union_S3 of the 2-dimensional topological manifold T,
             Subset S1_union_S4 of the 2-dimensional topological manifold T,
             Subset S2_union_S5 of the 2-dimensional topological manifold T]
            sage: union_S_i = S[0].union(S[1:], name='union_S_i'); union_S_i
            Subset union_S_i of the 2-dimensional topological manifold T
            sage: T.subset_family()
            Set {S0, S0_union_S3, S0_union_S3_union_S1_union_S4, S1,
                 S1_union_S4, S2, S2_union_S5, S3, S4, S5, T, union_S_i}
             of subsets of the 2-dimensional topological manifold T

        .. PLOT::

            def label(element):
                if isinstance(element, str):
                    return element
                try:
                    return element._name.replace('_union_', '∪')
                except AttributeError:
                    return '[' + ', '.join(sorted(label(x) for x in element)) + ']'

            M = Manifold(2, 'M', structure='topological')
            a = M.subset('A')
            b = M.subset('B')
            c = a.union(b); c
            P = M.subset_poset(open_covers=True)
            g1 = P.plot(element_labels={element: label(element) for element in P})

            T = Manifold(2, 'T', structure='topological')
            from sage.typeset.unicode_art import unicode_subscript
            S = [T.subset(f'S{unicode_subscript(i)}') for i in range(6)]
            [S[i].union(S[i+3]) for i in range(3)]
            union_S_i = S[0].union(S[1:], name='⋃ᵢSᵢ'); union_S_i
            P = T.subset_poset(open_covers=True)
            g2 = P.plot(element_labels={element: label(element) for element in P})

            sphinx_plot(graphics_array([g1, g2]), figsize=(8, 3))

        TESTS::

            sage: a.is_subset(a.union(b))
            True
            sage: b.is_subset(a.union(b))
            True
            sage: a.union(b) is b.union(a)
            True
            sage: a.union(a.union(b)) is a.union(b)
            True
            sage: (a.union(b)).union(a) is a.union(b)
            True
            sage: a.union(M) is M
            True
            sage: M.union(a) is M
            True

        Check that :issue:`30401` is fixed::

            sage: d = a.subset('D')
            sage: e = a.subset('E')
            sage: d.union(e).is_subset(a)
            True
        """
        subsets = ManifoldSubsetFiniteFamily.from_subsets_or_families(self, *others)
        subsets = self._reduce_union_members(subsets)
        assert subsets
        subset_iter = iter(subsets)
        res = next(subset_iter)
        others = list(subset_iter)
        if not others:
            return res
        for other in others[:-1]:
            res = res._union_subset(other)
        # The last one gets the name
        return res._union_subset(others[-1], name=name, latex_name=latex_name)

    @staticmethod
    def _reduce_union_members(subsets):
        r"""
        Return a reduced set of subsets with the same union as the given subsets.

        It is reduced with respect to two operations:

        - replacing an inclusion chain by its maximal element

        - replacing a pair of subsets with a declared union by the union

        INPUT:

        - ``subsets`` -- an iterable of :class:`ManifoldSubset` instances
          of the same manifold

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: B1 = A.subset('B1')
            sage: B2 = A.subset('B2')
            sage: B = B1.union(B2)
            sage: M._reduce_union_members([])
            {}
            sage: M._reduce_union_members([B1, B])
            Set {B1_union_B2} of subsets of the 2-dimensional topological manifold M
            sage: M._reduce_union_members([A, B1, B2])
            Set {A} of subsets of the 2-dimensional topological manifold M
        """
        subsets = set(subsets)

        def reduce():
            # Greedily replace inclusion chains by their maximal element
            # and pairs with declared unions by their union
            for U, V in itertools.combinations(subsets, 2):
                if U.is_subset(V):
                    subsets.remove(U)
                    return True
                if V.is_subset(U):
                    subsets.remove(V)
                    return True
                try:
                    UV = U._unions[V._name]
                except KeyError:
                    pass
                else:
                    subsets.difference_update([U, V])
                    subsets.add(UV)
                    return True
            return False
        while reduce():
            pass
        return ManifoldSubsetFiniteFamily(subsets)

    def _union_subset(self, other, name=None, latex_name=None):
        r"""
        Return a subset of the manifold that is the union of ``self`` and ``other``.

        The result is always a new subset of the manifold and is also
        stored in ``self`` and ``other``'s dictionaries of known unions.

        INPUT:

        - ``other`` -- an instance of :class:`ManifoldSubset`
        - ``name`` -- (default: ``None``) name given to the union; the default is
          ``self._name`` union ``other._name``
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          union; the default is built upon the symbol `\cup`

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: B1 = M.subset('B1')
            sage: B2 = M.subset('B2')
            sage: B1._union_subset(B2)
            Subset B1_union_B2 of the 2-dimensional topological manifold M
        """
        if latex_name is None:
            if name is None:
                latex_name = r'\cup '.join(S._latex_name for S in (self, other))
            else:
                latex_name = name
        if name is None:
            name = "_union_".join(S._name for S in (self, other))
        res_open = all(S.is_open() for S in (self, other))
        res = self.superset(name, latex_name, is_open=res_open)
        res.declare_superset(other)
        res._top_subsets.add(self)
        res._top_subsets.add(other)
        self._unions[other._name] = other._unions[self._name] = res
        for sp in self._supersets:
            if sp in other._supersets:
                sp._subsets.add(res)
                res._supersets.add(sp)
        if res._is_open:
            for chart in other._atlas:
                if chart not in res._atlas:
                    res._atlas.append(chart)
            for chart in other._top_charts:
                if chart not in res._top_charts:
                    res._top_charts.append(chart)
            res._coord_changes.update(other._coord_changes)
        # Open covers of the union:
        for oc1 in self._open_covers:
            for oc2 in other._open_covers:
                oc = oc1[:]
                for s in oc2:
                    if s not in oc:
                        oc.append(s)
            res._open_covers.append(oc)
        return res

    def complement(self, superset=None, name=None, latex_name=None, is_open=False):
        r"""
        Return the complement of ``self`` in the manifold or in ``superset``.

        INPUT:

        - ``superset`` -- (default: ``self.manifold()``) a superset of ``self``
        - ``name`` -- (default: ``None``) name given to the complement in the
          case the latter has to be created; the default is
          ``superset._name`` minus ``self._name``
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          complement in the case the latter has to be created; the default
          is built upon the symbol `\setminus`
        - ``is_open`` -- boolean (default: ``False``); if ``True``, the created
          subset is assumed to be open with respect to the manifold's topology

        OUTPUT:

        - instance of :class:`ManifoldSubset` representing the subset that
          is ``superset`` minus ``self``

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: B1 = A.subset('B1')
            sage: B2 = A.subset('B2')
            sage: B1.complement()
            Subset M_minus_B1 of the 2-dimensional topological manifold M
            sage: B1.complement(A)
            Subset A_minus_B1 of the 2-dimensional topological manifold M
            sage: B1.complement(B2)
            Traceback (most recent call last):
            ...
            TypeError: superset must be a superset of self

        Demanding that the complement is open makes ``self`` a closed subset::

            sage: A.is_closed()  # False a priori
            False
            sage: A.complement(is_open=True)
            Open subset M_minus_A of the 2-dimensional topological manifold M
            sage: A.is_closed()
            True
        """
        if superset is None:
            superset = self.manifold()
        elif not self.is_subset(superset):
            raise TypeError("superset must be a superset of self")
        return superset.difference(self,
                                   name=name, latex_name=latex_name,
                                   is_open=is_open)

    def difference(self, other, name=None, latex_name=None, is_open=False):
        r"""
        Return the set difference of ``self`` minus ``other``.

        INPUT:

        - ``other`` -- another subset of the same manifold
        - ``name`` -- (default: ``None``) name given to the difference in the
          case the latter has to be created; the default is
          ``self._name`` minus ``other._name``
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          difference in the case the latter has to be created; the default
          is built upon the symbol `\setminus`
        - ``is_open`` -- boolean (default: ``False``); if ``True``, the created
          subset is assumed to be open with respect to the manifold's topology

        OUTPUT:

        - instance of :class:`ManifoldSubset` representing the subset that is
          ``self`` minus ``other``

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.subset('A')
            sage: CA = M.difference(A); CA
            Subset M_minus_A of the 2-dimensional topological manifold M
            sage: latex(CA)
            M\setminus A
            sage: A.intersection(CA).is_empty()
            True
            sage: A.union(CA)
            2-dimensional topological manifold M

            sage: O = M.open_subset('O')
            sage: CO = M.difference(O); CO
            Subset M_minus_O of the 2-dimensional topological manifold M
            sage: M.difference(O) is CO
            True

            sage: CO2 = M.difference(O, is_open=True, name='CO2'); CO2
            Open subset CO2 of the 2-dimensional topological manifold M
            sage: CO is CO2
            False
            sage: CO.is_subset(CO2) and CO2.is_subset(CO)
            True
            sage: M.difference(O, is_open=True)
            Open subset CO2 of the 2-dimensional topological manifold M

        Since `O` is open and we have asked `M\setminus O` to be open, `O`
        is a clopen set (if `O\neq M` and `O\neq\emptyset`, this implies that
        `M` is not connected)::

            sage: O.is_closed() and O.is_open()
            True
        """
        # See if it has been created already
        diffs = []
        for diff_name, intersection in other._intersections.items():
            if intersection.is_empty():
                try:
                    union = other._unions[diff_name]
                except KeyError:
                    pass
                else:
                    if union == self:
                        diff = self.subset_family()[diff_name]
                        if not is_open:
                            return diff
                        if diff.is_open():
                            return diff
                        # is_open=True but we found a subset that
                        # is not known to be open - and we cannot
                        # declare it open.
                        diffs.append(diff)

        if latex_name is None:
            if name is None:
                latex_name = r'\setminus '.join(S._latex_name for S in (self, other))
            else:
                latex_name = name
        if name is None:
            name = "_minus_".join(S._name for S in (self, other))

        is_open = is_open or (self.is_open() and other.is_closed())

        diff = self.subset(name=name, latex_name=latex_name, is_open=is_open)
        diff.declare_equal(diffs)
        self.declare_union(other, diff, disjoint=True)
        return diff

    def closure(self, name=None, latex_name=None):
        r"""
        Return the topological closure of ``self`` as a subset of the manifold.

        INPUT:

        - ``name`` -- (default: ``None``) name given to the difference in the
          case the latter has to be created; the default prepends ``cl_``
          to ``self._name``
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          difference in the case the latter has to be created; the default
          is built upon the operator `\mathrm{cl}`

        OUTPUT:

        - if ``self`` is already known to be closed (see :meth:`is_closed`),
          ``self``; otherwise, an instance of
          :class:`~sage.manifolds.subsets.closure.ManifoldSubsetClosure`

        EXAMPLES::

            sage: M = Manifold(2, 'R^2', structure='topological')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: M.closure() is M
            True
            sage: D2 = M.open_subset('D2', coord_def={c_cart: x^2+y^2<2}); D2
            Open subset D2 of the 2-dimensional topological manifold R^2
            sage: cl_D2 = D2.closure(); cl_D2
            Topological closure cl_D2 of the
             Open subset D2 of the 2-dimensional topological manifold R^2
            sage: cl_D2.is_closed()
            True
            sage: cl_D2 is cl_D2.closure()
            True

            sage: D1 = D2.open_subset('D1'); D1
            Open subset D1 of the 2-dimensional topological manifold R^2
            sage: D1.closure().is_subset(D2.closure())
            True
        """
        if self.is_closed():
            return self
        from .subsets.closure import ManifoldSubsetClosure
        return ManifoldSubsetClosure(self, name=name, latex_name=latex_name)

    #### End of construction of new sets from self
