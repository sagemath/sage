# sage_setup: distribution = sagemath-symbolics
r"""
Topological Manifolds

Given a topological field `K` (in most applications, `K = \RR` or
`K = \CC`) and a nonnegative integer `n`, a *topological manifold of
dimension* `n` *over K* is a topological space `M` such that

- `M` is a Hausdorff space,
- `M` is second countable,
- every point in `M` has a neighborhood homeomorphic to `K^n`.

Topological manifolds are implemented via the class
:class:`TopologicalManifold`. Open subsets of topological manifolds
are also implemented via :class:`TopologicalManifold`, since they are
topological manifolds by themselves.

In the current setting, topological manifolds are mostly described by
means of charts (see :class:`~sage.manifolds.chart.Chart`).

:class:`TopologicalManifold` serves as a base class for more specific
manifold classes.

The user interface is provided by the generic function
:func:`~sage.manifolds.manifold.Manifold`, with
with the argument ``structure`` set to ``'topological'``.

.. RUBRIC:: Example 1: the 2-sphere as a topological manifold of dimension
  2 over `\RR`

One starts by declaring `S^2` as a 2-dimensional topological manifold::

    sage: M = Manifold(2, 'S^2', structure='topological')
    sage: M
    2-dimensional topological manifold S^2

Since the base topological field has not been specified in the argument list
of ``Manifold``, `\RR` is assumed::

    sage: M.base_field()
    Real Field with 53 bits of precision
    sage: dim(M)
    2

Let us consider the complement of a point, the "North pole" say; this is an
open subset of `S^2`, which we call `U`::

    sage: U = M.open_subset('U'); U
    Open subset U of the 2-dimensional topological manifold S^2

A standard chart on `U` is provided by the stereographic projection from the
North pole to the equatorial plane::

    sage: stereoN.<x,y> = U.chart(); stereoN
    Chart (U, (x, y))

Thanks to the operator ``<x,y>`` on the left-hand side, the coordinates
declared in a chart (here `x` and `y`), are accessible by their names;
they are Sage's symbolic variables::

    sage: y
    y
    sage: type(y)
    <class 'sage.symbolic.expression.Expression'>

The South pole is the point of coordinates `(x, y) = (0, 0)` in the above
chart::

    sage: S = U.point((0,0), chart=stereoN, name='S'); S
    Point S on the 2-dimensional topological manifold S^2

Let us call `V` the open subset that is the complement of the South pole and
let us introduce on it the chart induced by the stereographic projection from
the South pole to the equatorial plane::

    sage: V = M.open_subset('V'); V
    Open subset V of the 2-dimensional topological manifold S^2
    sage: stereoS.<u,v> = V.chart(); stereoS
    Chart (V, (u, v))

The North pole is the point of coordinates `(u, v) = (0, 0)` in this chart::

    sage: N = V.point((0,0), chart=stereoS, name='N'); N
    Point N on the 2-dimensional topological manifold S^2

To fully construct the manifold, we declare that it is the union of `U`
and `V`::

    sage: M.declare_union(U,V)

and we provide the transition map between the charts ``stereoN`` =
`(U, (x, y))` and ``stereoS`` = `(V, (u, v))`, denoting by `W` the
intersection of `U` and `V` (`W` is the subset of `U` defined by
`x^2 + y^2 \neq 0`, as well as the subset of `V` defined by
`u^2 + v^2 \neq 0`)::

    sage: stereoN_to_S = stereoN.transition_map(stereoS, [x/(x^2+y^2), y/(x^2+y^2)],
    ....:                          intersection_name='W', restrictions1= x^2+y^2!=0,
    ....:                                                 restrictions2= u^2+v^2!=0)
    sage: stereoN_to_S
    Change of coordinates from Chart (W, (x, y)) to Chart (W, (u, v))
    sage: stereoN_to_S.display()
    u = x/(x^2 + y^2)
    v = y/(x^2 + y^2)

We give the name ``W`` to the Python variable representing `W = U \cap V`::

    sage: W = U.intersection(V)

The inverse of the transition map is computed by the method
:meth:`sage.manifolds.chart.CoordChange.inverse`::

    sage: stereoN_to_S.inverse()
    Change of coordinates from Chart (W, (u, v)) to Chart (W, (x, y))
    sage: stereoN_to_S.inverse().display()
    x = u/(u^2 + v^2)
    y = v/(u^2 + v^2)

At this stage, we have four open subsets on `S^2`::

    sage: M.subset_family()
    Set {S^2, U, V, W} of open subsets of the 2-dimensional topological manifold S^2

`W` is the open subset that is the complement of the two poles::

    sage: N in W or S in W
    False

The North pole lies in `V` and the South pole in `U`::

    sage: N in V, N in U
    (True, False)
    sage: S in U, S in V
    (True, False)

The manifold's (user) atlas contains four charts, two of them
being restrictions of charts to a smaller domain::

    sage: M.atlas()
    [Chart (U, (x, y)), Chart (V, (u, v)),
     Chart (W, (x, y)), Chart (W, (u, v))]

Let us consider the point of coordinates `(1, 2)` in the chart ``stereoN``::

    sage: p = M.point((1,2), chart=stereoN, name='p'); p
    Point p on the 2-dimensional topological manifold S^2
    sage: p.parent()
    2-dimensional topological manifold S^2
    sage: p in W
    True

The coordinates of `p` in the chart ``stereoS`` are computed by letting
the chart act on the point::

    sage: stereoS(p)
    (1/5, 2/5)

Given the definition of `p`, we have of course::

    sage: stereoN(p)
    (1, 2)

Similarly::

    sage: stereoS(N)
    (0, 0)
    sage: stereoN(S)
    (0, 0)

A continuous map `S^2 \to \RR` (scalar field)::

    sage: f = M.scalar_field({stereoN: atan(x^2+y^2), stereoS: pi/2-atan(u^2+v^2)},
    ....:                    name='f')
    sage: f
    Scalar field f on the 2-dimensional topological manifold S^2
    sage: f.display()
    f: S^2 → ℝ
    on U: (x, y) ↦ arctan(x^2 + y^2)
    on V: (u, v) ↦ 1/2*pi - arctan(u^2 + v^2)
    sage: f(p)
    arctan(5)
    sage: f(N)
    1/2*pi
    sage: f(S)
    0
    sage: f.parent()
    Algebra of scalar fields on the 2-dimensional topological manifold S^2
    sage: f.parent().category()
    Join of Category of commutative algebras over Symbolic Ring and Category of homsets of topological spaces


.. RUBRIC:: Example 2: the Riemann sphere as a topological manifold of
  dimension 1 over `\CC`

We declare the Riemann sphere `\CC^*` as a 1-dimensional topological manifold
over `\CC`::

    sage: M = Manifold(1, 'ℂ*', structure='topological', field='complex'); M
    Complex 1-dimensional topological manifold ℂ*

We introduce a first open subset, which is actually
`\CC = \CC^*\setminus\{\infty\}` if we interpret `\CC^*` as the
Alexandroff one-point compactification of `\CC`::

    sage: U = M.open_subset('U')

A natural chart on `U` is then nothing but the identity map of `\CC`, hence
we denote the associated coordinate by `z`::

    sage: Z.<z> = U.chart()

The origin of the complex plane is the point of coordinate `z = 0`::

    sage: O = U.point((0,), chart=Z, name='O'); O
    Point O on the Complex 1-dimensional topological manifold ℂ*

Another open subset of `\CC^*` is `V = \CC^*\setminus\{O\}`::

    sage: V = M.open_subset('V')

We define a chart on `V` such that the point at infinity is the point of
coordinate `0` in this chart::

    sage: W.<w> = V.chart(); W
    Chart (V, (w,))
    sage: inf = M.point((0,), chart=W, name='inf', latex_name=r'\infty')
    sage: inf
    Point inf on the Complex 1-dimensional topological manifold ℂ*

To fully construct the Riemann sphere, we declare that it is the union
of `U` and `V`::

    sage: M.declare_union(U,V)

and we provide the transition map between the two charts as `w = 1 / z`
on `A = U \cap V`::

    sage: Z_to_W = Z.transition_map(W, 1/z, intersection_name='A',
    ....:                           restrictions1= z!=0, restrictions2= w!=0)
    sage: Z_to_W
    Change of coordinates from Chart (A, (z,)) to Chart (A, (w,))
    sage: Z_to_W.display()
    w = 1/z
    sage: Z_to_W.inverse()
    Change of coordinates from Chart (A, (w,)) to Chart (A, (z,))
    sage: Z_to_W.inverse().display()
    z = 1/w

Let consider the complex number `i` as a point of the Riemann sphere::

    sage: i = M((I,), chart=Z, name='i'); i
    Point i on the Complex 1-dimensional topological manifold ℂ*

Its coordinates w.r.t. the charts ``Z`` and ``W`` are::

    sage: Z(i)
    (I,)
    sage: W(i)
    (-I,)

and we have::

    sage: i in U
    True
    sage: i in V
    True

The following subsets and charts have been defined::

    sage: M.subset_family()
    Set {A, U, V, ℂ*} of open subsets of the Complex 1-dimensional topological manifold ℂ*
    sage: M.atlas()
    [Chart (U, (z,)), Chart (V, (w,)), Chart (A, (z,)), Chart (A, (w,))]

A constant map `\CC^* \rightarrow \CC`::

    sage: f = M.constant_scalar_field(3+2*I, name='f'); f
    Scalar field f on the Complex 1-dimensional topological manifold ℂ*
    sage: f.display()
    f: ℂ* → ℂ
    on U: z ↦ 2*I + 3
    on V: w ↦ 2*I + 3
    sage: f(O)
    2*I + 3
    sage: f(i)
    2*I + 3
    sage: f(inf)
    2*I + 3
    sage: f.parent()
    Algebra of scalar fields on the Complex 1-dimensional topological
     manifold ℂ*
    sage: f.parent().category()
    Join of Category of commutative algebras over Symbolic Ring and Category of homsets of topological spaces

AUTHORS:

- Eric Gourgoulhon (2015): initial version
- Travis Scrimshaw (2015): structure described via
  :class:`~sage.manifolds.structure.TopologicalStructure` or
  :class:`~sage.manifolds.structure.RealTopologicalStructure`
- Michael Jung (2020): topological vector bundles and orientability


REFERENCES:

- [Lee2011]_
- [Lee2013]_
- [KN1963]_
- [Huy2005]_
"""

# ****************************************************************************
#       Copyright (C) 2015-2020 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015      Travis Scrimshaw <tscrimsh@umn.edu>
#       Copyright (C) 2016      Andrew Mathas
#       Copyright (C) 2018      Florentin Jaffredo
#       Copyright (C) 2019      Hans Fotsing Tetsing
#       Copyright (C) 2019-2020 Michael Jung
#       Copyright (C) 2021      Matthias Koeppe <mkoeppe@math.ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from __future__ import annotations

from typing import TYPE_CHECKING, Optional, Union, overload

import sage.rings.abc
from sage.categories.fields import Fields
from sage.categories.homset import Hom
from sage.categories.manifolds import Manifolds
from sage.manifolds.structure import (
    DifferentialStructure,
    RealDifferentialStructure,
    RealTopologicalStructure,
    TopologicalStructure,
)
from sage.manifolds.subset import ManifoldSubset
from sage.misc.cachefunc import cached_method
from sage.misc.prandom import getrandbits
from sage.rings.cc import CC
from sage.rings.integer import Integer
from sage.rings.real_mpfr import RR
from sage.structure.global_options import GlobalOptions

if TYPE_CHECKING:
    from sage.manifolds.chart import Chart
    from sage.manifolds.continuous_map import ContinuousMap
    from sage.manifolds.differentiable.diff_map import DiffMap
    from sage.manifolds.differentiable.manifold import DifferentiableManifold
    from sage.manifolds.scalarfield import ScalarField


#############################################################################
# Global options

#############################################################################
# Class

class TopologicalManifold(ManifoldSubset):
    r"""
    Topological manifold over a topological field `K`.

    Given a topological field `K` (in most applications, `K = \RR` or
    `K = \CC`) and a nonnegative integer `n`, a *topological manifold of
    dimension* `n` *over K* is a topological space `M` such that

    - `M` is a Hausdorff space,
    - `M` is second countable, and
    - every point in `M` has a neighborhood homeomorphic to `K^n`.

    This is a Sage *parent* class, the corresponding *element*
    class being :class:`~sage.manifolds.point.ManifoldPoint`.

    INPUT:

    - ``n`` -- positive integer; dimension of the manifold
    - ``name`` -- string; name (symbol) given to the manifold
    - ``field`` -- field `K` on which the manifold is
      defined; allowed values are

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

    A 4-dimensional topological manifold (over `\RR`)::

        sage: M = Manifold(4, 'M', latex_name=r'\mathcal{M}', structure='topological')
        sage: M
        4-dimensional topological manifold M
        sage: latex(M)
        \mathcal{M}
        sage: type(M)
        <class 'sage.manifolds.manifold.TopologicalManifold_with_category'>
        sage: M.base_field()
        Real Field with 53 bits of precision
        sage: dim(M)
        4

    The input parameter ``start_index`` defines the range of indices
    on the manifold::

        sage: M = Manifold(4, 'M', structure='topological')
        sage: list(M.irange())
        [0, 1, 2, 3]
        sage: M = Manifold(4, 'M', structure='topological', start_index=1)
        sage: list(M.irange())
        [1, 2, 3, 4]
        sage: list(Manifold(4, 'M', structure='topological', start_index=-2).irange())
        [-2, -1, 0, 1]

    A complex manifold::

        sage: N = Manifold(3, 'N', structure='topological', field='complex'); N
        Complex 3-dimensional topological manifold N

    A manifold over `\QQ`::

        sage: N = Manifold(6, 'N', structure='topological', field=QQ); N
        6-dimensional topological manifold N over the Rational Field

    A manifold over `\QQ_5`, the field of 5-adic numbers::

        sage: N = Manifold(2, 'N', structure='topological', field=Qp(5)); N             # needs sage.rings.padics
        2-dimensional topological manifold N over the 5-adic Field with capped
         relative precision 20

    A manifold is a Sage *parent* object, in the category of topological
    manifolds over a given topological field (see
    :class:`~sage.categories.manifolds.Manifolds`)::

        sage: isinstance(M, Parent)
        True
        sage: M.category()
        Category of manifolds over Real Field with 53 bits of precision
        sage: from sage.categories.manifolds import Manifolds
        sage: M.category() is Manifolds(RR)
        True
        sage: M.category() is Manifolds(M.base_field())
        True
        sage: M in Manifolds(RR)
        True
        sage: N in Manifolds(Qp(5))                                                     # needs sage.rings.padics
        True

    The corresponding Sage *elements* are points::

        sage: X.<t, x, y, z> = M.chart()
        sage: p = M.an_element(); p
        Point on the 4-dimensional topological manifold M
        sage: p.parent()
        4-dimensional topological manifold M
        sage: M.is_parent_of(p)
        True
        sage: p in M
        True

    The manifold's points are instances of class
    :class:`~sage.manifolds.point.ManifoldPoint`::

        sage: isinstance(p, sage.manifolds.point.ManifoldPoint)
        True

    Since an open subset of a topological manifold `M` is itself a
    topological manifold, open subsets of `M` are instances of the class
    :class:`TopologicalManifold`::

        sage: U = M.open_subset('U'); U
        Open subset U of the 4-dimensional topological manifold M
        sage: isinstance(U, sage.manifolds.manifold.TopologicalManifold)
        True
        sage: U.base_field() == M.base_field()
        True
        sage: dim(U) == dim(M)
        True
        sage: U.category()
        Join of Category of subobjects of sets and Category of manifolds over
         Real Field with 53 bits of precision

    The manifold passes all the tests of the test suite relative to its
    category::

        sage: TestSuite(M).run()

    .. SEEALSO::

        :mod:`sage.manifolds.manifold`
    """
    _dim: int

    def __init__(self, n, name, field, structure, base_manifold=None,
                 latex_name=None, start_index=0, category=None,
                 unique_tag=None):
        r"""
        Construct a topological manifold.

        TESTS::

            sage: M = Manifold(3, 'M', latex_name=r'\mathbb{M}',
            ....:              structure='topological', start_index=1)
            sage: M
            3-dimensional topological manifold M
            sage: latex(M)
            \mathbb{M}
            sage: dim(M)
            3
            sage: X.<x,y,z> = M.chart()
            sage: TestSuite(M).run()

        Tests for open subsets::

            sage: U = M.open_subset('U', coord_def={X: x>0})
            sage: TestSuite(U).run()
            sage: U.category() is M.category().Subobjects()
            True
        """
        # Initialization of the attributes _dim, _field, _field_type:
        if not isinstance(n, (int, Integer)):
            raise TypeError("the manifold dimension must be an integer")
        if n < 1:
            raise ValueError("the manifold dimension must be strictly positive")
        self._dim = n
        if field == 'real':
            self._field = RR
            self._field_type = 'real'
        elif field == 'complex':
            self._field = CC
            self._field_type = 'complex'
        else:
            if field not in Fields():
                raise TypeError("the argument 'field' must be a field")
            self._field = field
            if isinstance(field, sage.rings.abc.RealField):
                self._field_type = 'real'
            elif isinstance(field, sage.rings.abc.ComplexField):
                self._field_type = 'complex'
            else:
                self._field_type = 'neither_real_nor_complex'
        # Structure and category:
        self._structure = structure
        if base_manifold is None:
            base_manifold = self
            category = Manifolds(self._field).or_subcategory(category)
            category = self._structure.subcategory(category)
        else:
            category = base_manifold.category().Subobjects()
        # Initialization as a manifold set:
        ManifoldSubset.__init__(self, base_manifold, name, latex_name=latex_name,
                                category=category)
        self._is_open = True
        self._open_covers.append([self])  # list of open covers of self

        if not isinstance(start_index, (int, Integer)):
            raise TypeError("the starting index must be an integer")
        self._sindex = start_index

        self._atlas = []  # list of charts defined on subsets of self
        self._top_charts = []  # list of charts defined on subsets of self
                        # that are not subcharts of charts on larger subsets
        self._def_chart = None  # default chart
        self._orientation = [] # set no orientation a priori
        self._charts_by_coord = {} # dictionary of charts whose domain is self
                                   # (key: string formed by the coordinate
                                   #  symbols separated by a white space)
        self._coord_changes = {} # dictionary of transition maps (key: pair of
                                 # of charts)
        # List of charts that individually cover self, i.e. whose
        # domains are self (if non-empty, self is a coordinate domain):
        self._covering_charts = []
        # Algebra of scalar fields defined on self:
        self._scalar_field_algebra = self._structure.scalar_field_algebra(self)
        # The zero scalar field:
        self._zero_scalar_field = self._scalar_field_algebra.zero()
        # The unit scalar field:
        self._one_scalar_field = self._scalar_field_algebra.one()
        # The current calculus method on the manifold
        #   (to be changed by set_calculus_method)
        self._calculus_method = 'SR'

    def _repr_(self):
        r"""
        Return a string representation of the manifold.

        TESTS::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: M._repr_()
            '3-dimensional topological manifold M'
            sage: repr(M)  # indirect doctest
            '3-dimensional topological manifold M'
            sage: M  # indirect doctest
            3-dimensional topological manifold M
            sage: M = Manifold(3, 'M', structure='topological', field='complex')
            sage: M._repr_()
            'Complex 3-dimensional topological manifold M'
            sage: M = Manifold(3, 'M', structure='topological', field=QQ)
            sage: M._repr_()
            '3-dimensional topological manifold M over the Rational Field'

        If the manifold is actually an open subset of a larger manifold, the
        string representation is different::

            sage: U = M.open_subset('U')
            sage: U._repr_()
            'Open subset U of the 3-dimensional topological manifold M
             over the Rational Field'
        """
        if self is self._manifold:
            if self._field_type == 'real':
                return "{}-dimensional {} manifold {}".format(self._dim,
                                                          self._structure.name,
                                                          self._name)
            elif self._field_type == 'complex':
                if isinstance(self._structure, DifferentialStructure):
                    return "{}-dimensional complex manifold {}".format(
                                                                    self._dim,
                                                                    self._name)
                else:
                    return "Complex {}-dimensional {} manifold {}".format(
                                                          self._dim,
                                                          self._structure.name,
                                                          self._name)
            return "{}-dimensional {} manifold {} over the {}".format(
                                                          self._dim,
                                                          self._structure.name,
                                                          self._name,
                                                          self._field)
        else:
            return "Open subset {} of the {}".format(self._name, self._manifold)

    def _an_element_(self):
        r"""
        Construct some point on the manifold.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: p = M._an_element_(); p
            Point on the 2-dimensional topological manifold M
            sage: p.coord()
            (0, 0)
            sage: U = M.open_subset('U', coord_def={X: y>1}); U
            Open subset U of the 2-dimensional topological manifold M
            sage: p = U._an_element_(); p
            Point on the 2-dimensional topological manifold M
            sage: p in U
            True
            sage: p.coord()
            (0, 2)
            sage: V = U.open_subset('V', coord_def={X.restrict(U): x<-pi})
            sage: p = V._an_element_(); p
            Point on the 2-dimensional topological manifold M
            sage: p in V
            True
            sage: p.coord()
            (-pi - 1, 2)
        """
        from sage.rings.infinity import Infinity
        if self._def_chart is None:
            return self.element_class(self)
        # Attempt to construct a point in the domain of the default chart
        chart = self._def_chart
        if self._field_type == 'real':
            coords = []
            for coord_range in chart._bounds:
                xmin = coord_range[0][0]
                xmax = coord_range[1][0]
                if xmin == -Infinity:
                    if xmax == Infinity:
                        x = 0
                    else:
                        x = xmax - 1
                else:
                    if xmax == Infinity:
                        x = xmin + 1
                    else:
                        x = (xmin + xmax)/2
                coords.append(x)
        else:
            coords = self._dim*[0]
        if not chart.valid_coordinates(*coords):
            # Attempt to construct a point in the domain of other charts
            if self._field_type == 'real':
                for ch in self._atlas:
                    if ch is self._def_chart:
                        continue # since this case has already been attempted
                    coords = []
                    for coord_range in ch._bounds:
                        xmin = coord_range[0][0]
                        xmax = coord_range[1][0]
                        if xmin == -Infinity:
                            if xmax == Infinity:
                                x = 0
                            else:
                                x = xmax - 1
                        else:
                            if xmax == Infinity:
                                x = xmin + 1
                            else:
                                x = (xmin + xmax)/2
                        coords.append(x)
                    if ch.valid_coordinates(*coords):
                        chart = ch
                        break
                else:
                    # A generic element with specific coordinates could not be
                    # automatically generated, due to too complex coordinate
                    # conditions. An element without any coordinate set is
                    # returned instead:
                    return self.element_class(self)
            else:
                # Case of manifolds over a field different from R
                for ch in self._atlas:
                    if ch is self._def_chart:
                        continue # since this case has already been attempted
                    if ch.valid_coordinates(*coords):
                        chart = ch
                        break
                else:
                    return self.element_class(self)
        # The point is constructed with check_coords=False since the check
        # has just been performed above:
        return self.element_class(self, coords=coords, chart=chart,
                                  check_coords=False)

    def __contains__(self, point):
        r"""
        Check whether a point is contained in the manifold.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((1,2), chart=X)
            sage: M.__contains__(p)
            True
            sage: p in M  # indirect doctest
            True
            sage: U = M.open_subset('U', coord_def={X: x>0})
            sage: U.__contains__(p)
            True
            sage: p in U  # indirect doctest
            True
            sage: V = U.open_subset('V', coord_def={X.restrict(U): y<0})
            sage: V.__contains__(p)
            False
            sage: p in V  # indirect doctest
            False
        """
        # for efficiency, a quick test first:
        if point.parent() is self:
            return True
        if point.parent().is_subset(self):
            return True
        for chart in self._atlas:
            if chart in point._coordinates:
                if chart.valid_coordinates( *(point._coordinates[chart]) ):
                    return True
        for chart in point._coordinates:
            for schart in chart._subcharts:
                if schart in self._atlas and schart.valid_coordinates(
                                          *(point._coordinates[chart]) ):
                    return True
        return False

    def open_subset(self, name, latex_name=None, coord_def={}, supersets=None):
        r"""
        Create an open subset of the manifold.

        An open subset is a set that is (i) included in the manifold and (ii)
        open with respect to the manifold's topology. It is a topological
        manifold by itself. Hence the returned object is an instance of
        :class:`TopologicalManifold`.

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

        OUTPUT: the open subset, as an instance of :class:`TopologicalManifold`

        EXAMPLES:

        Creating an open subset of a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: A = M.open_subset('A'); A
            Open subset A of the 2-dimensional topological manifold M

        As an open subset of a topological manifold, ``A`` is itself a
        topological manifold, on the same topological field and of the same
        dimension as ``M``::

            sage: isinstance(A, sage.manifolds.manifold.TopologicalManifold)
            True
            sage: A.base_field() == M.base_field()
            True
            sage: dim(A) == dim(M)
            True
            sage: A.category() is M.category().Subobjects()
            True

        Creating an open subset of ``A``::

            sage: B = A.open_subset('B'); B
            Open subset B of the 2-dimensional topological manifold M

        We have then::

            sage: frozenset(A.subsets())  # random (set output)
            {Open subset B of the 2-dimensional topological manifold M,
             Open subset A of the 2-dimensional topological manifold M}
            sage: B.is_subset(A)
            True
            sage: B.is_subset(M)
            True

        Defining an open subset by some coordinate restrictions: the open
        unit disk in `\RR^2`::

            sage: M = Manifold(2, 'R^2', structure='topological')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: U = M.open_subset('U', coord_def={c_cart: x^2+y^2<1}); U
            Open subset U of the 2-dimensional topological manifold R^2

        Since the argument ``coord_def`` has been set, ``U`` is automatically
        provided with a chart, which is the restriction of the Cartesian one
        to ``U``::

            sage: U.atlas()
            [Chart (U, (x, y))]

        Therefore, one can immediately check whether a point belongs
        to ``U``::

            sage: M.point((0,0)) in U
            True
            sage: M.point((1/2,1/3)) in U
            True
            sage: M.point((1,2)) in U
            False
        """
        resu = TopologicalManifold(self._dim, name, self._field,
                                   self._structure,
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

            sage: M = Manifold(2, 'R^2', structure='topological')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: from sage.manifolds.manifold import TopologicalManifold
            sage: U = TopologicalManifold(2, 'U', field=M._field, structure=M._structure, base_manifold=M)
            sage: M._init_open_subset(U, coord_def={c_cart: x^2+y^2<1})
            sage: U
            Open subset U of the 2-dimensional topological manifold R^2
        """
        resu._calculus_method = self._calculus_method
        if self.is_empty():
            self.declare_equal(resu)
        else:
            self.declare_superset(resu)
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

    def get_chart(self, coordinates, domain=None):
        r"""
        Get a chart from its coordinates.

        The chart must have been previously created by the method
        :meth:`chart`.

        INPUT:

        - ``coordinates`` -- single string composed of the coordinate symbols
          separated by a space
        - ``domain`` -- (default: ``None``) string containing the name of the
          chart's domain, which must be a subset of the current manifold; if
          ``None``, the current manifold is assumed

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.chart.Chart` (or of the subclass
          :class:`~sage.manifolds.chart.RealChart`) representing the chart
          corresponding to the above specifications

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: M.get_chart('x y')
            Chart (M, (x, y))
            sage: M.get_chart('x y') is X
            True
            sage: U = M.open_subset('U', coord_def={X: (y!=0,x<0)})
            sage: Y.<r, ph> = U.chart(r'r:(0,+oo) ph:(0,2*pi):\phi')
            sage: M.atlas()
            [Chart (M, (x, y)), Chart (U, (x, y)), Chart (U, (r, ph))]
            sage: M.get_chart('x y', domain='U')
            Chart (U, (x, y))
            sage: M.get_chart('x y', domain='U') is X.restrict(U)
            True
            sage: U.get_chart('r ph')
            Chart (U, (r, ph))
            sage: M.get_chart('r ph', domain='U')
            Chart (U, (r, ph))
            sage: M.get_chart('r ph', domain='U') is Y
            True
        """
        if domain is None:
            dom = self
        else:
            dom = self.get_subset(domain)
        try:
            return dom._charts_by_coord[coordinates]
        except KeyError:
            raise KeyError("the coordinates '{}' ".format(coordinates) +
                           "do not correspond to any chart with " +
                           "the {} as domain".format(dom))

    def dimension(self):
        r"""
        Return the dimension of the manifold over its base field.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: M.dimension()
            2

        A shortcut is ``dim()``::

            sage: M.dim()
            2

        The Sage global function ``dim`` can also be used::

            sage: dim(M)
            2
        """
        return self._dim

    dim = dimension

    def base_field(self):
        r"""
        Return the field on which the manifold is defined.

        OUTPUT: a topological field

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: M.base_field()
            Real Field with 53 bits of precision
            sage: M = Manifold(3, 'M', structure='topological', field='complex')
            sage: M.base_field()
            Complex Field with 53 bits of precision
            sage: M = Manifold(3, 'M', structure='topological', field=QQ)
            sage: M.base_field()
            Rational Field
        """
        return self._field

    def base_field_type(self):
        r"""
        Return the type of topological field on which the manifold is defined.

        OUTPUT:

        A string describing the field, with three possible values:

        - ``'real'`` for the real field `\RR`
        - ``'complex'`` for the complex field `\CC`
        - ``'neither_real_nor_complex'`` for a field different from `\RR` and `\CC`

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: M.base_field_type()
            'real'
            sage: M = Manifold(3, 'M', structure='topological', field='complex')
            sage: M.base_field_type()
            'complex'
            sage: M = Manifold(3, 'M', structure='topological', field=QQ)
            sage: M.base_field_type()
            'neither_real_nor_complex'
        """
        return self._field_type

    def start_index(self):
        r"""
        Return the first value of the index range used on the manifold.

        This is the parameter ``start_index`` passed at the construction of
        the manifold.

        OUTPUT:

        - the integer `i_0` such that all indices of indexed objects on the
          manifold range from `i_0` to `i_0 + n - 1`, where `n` is the
          manifold's dimension

        EXAMPLES::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: M.start_index()
            0
            sage: M = Manifold(3, 'M', structure='topological', start_index=1)
            sage: M.start_index()
            1
        """
        return self._sindex

    def irange(self, start=None, end=None):
        r"""
        Single index generator.

        INPUT:

        - ``start`` -- (default: ``None``) initial value `i_0` of the index;
          if ``None``, the value returned by :meth:`start_index()` is assumed
        - ``end`` -- (default: ``None``) final value `i_n` of the index;
          if ``None``, the value returned by :meth:`start_index()` plus
          `n - 1`, where `n` is the manifold dimension, is assumed

        OUTPUT: an iterable index, starting from `i_0` and ending at `i_0 + i_n`

        EXAMPLES:

        Index range on a 4-dimensional manifold::

            sage: M = Manifold(4, 'M', structure='topological')
            sage: list(M.irange())
            [0, 1, 2, 3]
            sage: list(M.irange(start=2))
            [2, 3]
            sage: list(M.irange(end=2))
            [0, 1, 2]
            sage: list(M.irange(start=1, end=2))
            [1, 2]

        Index range on a 4-dimensional manifold with starting index=1::

            sage: M = Manifold(4, 'M', structure='topological', start_index=1)
            sage: list(M.irange())
            [1, 2, 3, 4]
            sage: list(M.irange(start=2))
            [2, 3, 4]
            sage: list(M.irange(end=2))
            [1, 2]
            sage: list(M.irange(start=2, end=3))
            [2, 3]

        In general, one has always::

            sage: next(M.irange()) == M.start_index()
            True
        """
        si = self._sindex
        if start is None:
            i = si
        else:
            i = start
        if end is None:
            imax = self._dim + si
        else:
            imax = end + 1
        while i < imax:
            yield i
            i += 1

    def index_generator(self, nb_indices):
        r"""
        Generator of index series.

        INPUT:

        - ``nb_indices`` -- number of indices in a series

        OUTPUT:

        - an iterable index series for a generic component with the specified
          number of indices

        EXAMPLES:

        Indices on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological', start_index=1)
            sage: list(M.index_generator(2))
            [(1, 1), (1, 2), (2, 1), (2, 2)]

        Loops can be nested::

            sage: for ind1 in M.index_generator(2):
            ....:     print("{} : {}".format(ind1, list(M.index_generator(2))))
            (1, 1) : [(1, 1), (1, 2), (2, 1), (2, 2)]
            (1, 2) : [(1, 1), (1, 2), (2, 1), (2, 2)]
            (2, 1) : [(1, 1), (1, 2), (2, 1), (2, 2)]
            (2, 2) : [(1, 1), (1, 2), (2, 1), (2, 2)]
        """
        si = self._sindex
        imax = self._dim - 1 + si
        ind = [si for k in range(nb_indices)]
        ind_end = [si for k in range(nb_indices)]
        ind_end[0] = imax+1
        while ind != ind_end:
            yield tuple(ind)
            ret = 1
            for pos in range(nb_indices-1, -1, -1):
                if ind[pos] != imax:
                    ind[pos] += ret
                    ret = 0
                elif ret == 1:
                    if pos == 0:
                        ind[pos] = imax + 1 # end point reached
                    else:
                        ind[pos] = si
                        ret = 1

    def atlas(self) -> list[Chart]:
        r"""
        Return the list of charts that have been defined on the manifold.

        EXAMPLES:

        Let us consider `\RR^2` as a 2-dimensional manifold::

            sage: M = Manifold(2, 'R^2', structure='topological')

        Immediately after the manifold creation, the atlas is empty, since no
        chart has been defined yet::

            sage: M.atlas()
            []

        Let us introduce the chart of Cartesian coordinates::

            sage: c_cart.<x,y> = M.chart()
            sage: M.atlas()
            [Chart (R^2, (x, y))]

        The complement of the half line `\{y = 0, x \geq 0\}`::

            sage: U = M.open_subset('U', coord_def={c_cart: (y!=0,x<0)})
            sage: U.atlas()
            [Chart (U, (x, y))]
            sage: M.atlas()
            [Chart (R^2, (x, y)), Chart (U, (x, y))]

        Spherical (polar) coordinates on ``U``::

            sage: c_spher.<r, ph> = U.chart(r'r:(0,+oo) ph:(0,2*pi):\phi')
            sage: U.atlas()
            [Chart (U, (x, y)), Chart (U, (r, ph))]
            sage: M.atlas()
            [Chart (R^2, (x, y)), Chart (U, (x, y)), Chart (U, (r, ph))]

        .. SEEALSO::

            :meth:`top_charts`
        """
        return list(self._atlas) # Make a (shallow) copy

    def top_charts(self):
        r"""
        Return the list of charts defined on subsets of the current manifold
        that are not subcharts of charts on larger subsets.

        OUTPUT:

        - list of charts defined on open subsets of the manifold but not on
          larger subsets

        EXAMPLES:

        Charts on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: U = M.open_subset('U', coord_def={X: x>0})
            sage: Y.<u,v> = U.chart()
            sage: M.top_charts()
            [Chart (M, (x, y)), Chart (U, (u, v))]

        Note that the (user) atlas contains one more chart: ``(U, (x,y))``,
        which is not a "top" chart::

            sage: M.atlas()
            [Chart (M, (x, y)), Chart (U, (x, y)), Chart (U, (u, v))]

        .. SEEALSO::

            :meth:`atlas` for the complete list of charts defined on the
            manifold.
        """
        return list(self._top_charts) # Make a (shallow) copy

    def default_chart(self):
        r"""
        Return the default chart defined on the manifold.

        Unless changed via :meth:`set_default_chart`, the *default chart*
        is the first one defined on a subset of the manifold (possibly itself).

        OUTPUT:

        - instance of :class:`~sage.manifolds.chart.Chart`
          representing the default chart

        EXAMPLES:

        Default chart on a 2-dimensional manifold and on some subsets::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: M.chart('x y')
            Chart (M, (x, y))
            sage: M.chart('u v')
            Chart (M, (u, v))
            sage: M.default_chart()
            Chart (M, (x, y))
            sage: A = M.open_subset('A')
            sage: A.chart('t z')
            Chart (A, (t, z))
            sage: A.default_chart()
            Chart (A, (t, z))
        """
        return self._def_chart

    def set_default_chart(self, chart):
        r"""
        Changing the default chart on ``self``.

        INPUT:

        - ``chart`` -- a chart (must be defined on some subset ``self``)

        EXAMPLES:

        Charts on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: c_xy.<x,y> = M.chart()
            sage: c_uv.<u,v> = M.chart()
            sage: M.default_chart()
            Chart (M, (x, y))
            sage: M.set_default_chart(c_uv)
            sage: M.default_chart()
            Chart (M, (u, v))
        """
        from .chart import Chart
        if not isinstance(chart, Chart):
            raise TypeError("{} is not a chart".format(chart))
        if chart not in self._atlas:
            raise ValueError("the chart must be defined on the {}".format(self))
        self._def_chart = chart

    def coord_change(self, chart1, chart2):
        r"""
        Return the change of coordinates (transition map) between two charts
        defined on the manifold.

        The change of coordinates must have been defined previously, for
        instance by the method
        :meth:`~sage.manifolds.chart.Chart.transition_map`.

        INPUT:

        - ``chart1`` -- chart 1
        - ``chart2`` -- chart 2

        OUTPUT:

        - instance of :class:`~sage.manifolds.chart.CoordChange`
          representing the transition map from chart 1 to chart 2

        EXAMPLES:

        Change of coordinates on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: c_xy.<x,y> = M.chart()
            sage: c_uv.<u,v> = M.chart()
            sage: c_xy.transition_map(c_uv, (x+y, x-y)) # defines the coord. change
            Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))
            sage: M.coord_change(c_xy, c_uv) # returns the coord. change defined above
            Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))
        """
        if (chart1, chart2) not in self._coord_changes:
            raise TypeError("the change of coordinates from " +
                            "{} to {}".format(chart1, chart2) + " has not " +
                            "been defined on the {}".format(self))
        return self._coord_changes[(chart1, chart2)]

    def coord_changes(self):
        r"""
        Return the changes of coordinates (transition maps) defined on
        subsets of the manifold.

        OUTPUT: dictionary of changes of coordinates, with pairs of charts as keys

        EXAMPLES:

        Various changes of coordinates on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: c_xy.<x,y> = M.chart()
            sage: c_uv.<u,v> = M.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, [x+y, x-y])
            sage: M.coord_changes()
            {(Chart (M, (x, y)),
              Chart (M, (u, v))): Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))}
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: M.coord_changes()  # random (dictionary output)
            {(Chart (M, (u, v)),
              Chart (M, (x, y))): Change of coordinates from Chart (M, (u, v)) to Chart (M, (x, y)),
             (Chart (M, (x, y)),
              Chart (M, (u, v))): Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))}
            sage: c_rs.<r,s> = M.chart()
            sage: uv_to_rs = c_uv.transition_map(c_rs, [-u+2*v, 3*u-v])
            sage: M.coord_changes()  # random (dictionary output)
            {(Chart (M, (u, v)),
              Chart (M, (r, s))): Change of coordinates from Chart (M, (u, v)) to Chart (M, (r, s)),
             (Chart (M, (u, v)),
              Chart (M, (x, y))): Change of coordinates from Chart (M, (u, v)) to Chart (M, (x, y)),
             (Chart (M, (x, y)),
              Chart (M, (u, v))): Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))}
            sage: xy_to_rs = uv_to_rs * xy_to_uv
            sage: M.coord_changes()  # random (dictionary output)
            {(Chart (M, (u, v)),
              Chart (M, (r, s))): Change of coordinates from Chart (M, (u, v)) to Chart (M, (r, s)),
             (Chart (M, (u, v)),
              Chart (M, (x, y))): Change of coordinates from Chart (M, (u, v)) to Chart (M, (x, y)),
             (Chart (M, (x, y)),
              Chart (M, (u, v))): Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v)),
             (Chart (M, (x, y)),
              Chart (M, (r, s))): Change of coordinates from Chart (M, (x, y)) to Chart (M, (r, s))}
        """
        return self._coord_changes.copy()

    def is_manifestly_coordinate_domain(self):
        r"""
        Return ``True`` if the manifold is known to be the domain of some
        coordinate chart and ``False`` otherwise.

        If ``False`` is returned, either the manifold cannot be the domain of
        some coordinate chart or no such chart has been declared yet.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: U = M.open_subset('U')
            sage: X.<x,y> = U.chart()
            sage: U.is_manifestly_coordinate_domain()
            True
            sage: M.is_manifestly_coordinate_domain()
            False
            sage: Y.<u,v> = M.chart()
            sage: M.is_manifestly_coordinate_domain()
            True
        """
        return bool(self._covering_charts)

    def chart(
        self,
        coordinates: str = "",
        names=None,
        calc_method=None,
        coord_restrictions=None,
    ) -> Chart:
        r"""
        Define a chart, the domain of which is the manifold.

        A *chart* is a pair `(U, \varphi)`, where `U` is the current
        manifold and `\varphi: U \rightarrow V \subset K^n`
        is a homeomorphism from `U` to an open subset `V` of `K^n`, `K`
        being the field on which the manifold is defined.

        The components `(x^1, \ldots, x^n)` of `\varphi`, defined by
        `\varphi(p) = (x^1(p), \ldots, x^n(p)) \in K^n` for any point
        `p \in U`, are called the *coordinates* of the chart `(U, \varphi)`.

        See :class:`~sage.manifolds.chart.Chart` for a complete
        documentation.

        INPUT:

        - ``coordinates`` -- (default: ``''`` (empty string)) string
          defining the coordinate symbols, ranges and possible periodicities,
          see below
        - ``names`` -- (default: ``None``) unused argument, except if
          ``coordinates`` is not provided; it must then be a tuple containing
          the coordinate symbols (this is guaranteed if the shortcut operator
          ``<,>`` is used)
        - ``calc_method`` -- (default: ``None``) string defining the calculus
          method to be used on this chart; must be one of

          - ``'SR'``: Sage's default symbolic engine (Symbolic Ring)
          - ``'sympy'``: SymPy
          - ``None``: the current calculus method defined on the manifold is
            used (cf. :meth:`set_calculus_method`)
        - ``coord_restrictions`` -- additional restrictions on the coordinates.
          See below.

        The coordinates declared in the string ``coordinates`` are
        separated by ``' '`` (whitespace) and each coordinate has at most four
        fields, separated by a colon (``':'``):

        1. The coordinate symbol (a letter or a few letters).
        2. (optional, only for manifolds over `\RR`) The interval `I`
           defining the coordinate range: if not provided, the coordinate
           is assumed to span all `\RR`; otherwise `I` must be provided
           in the form ``(a,b)`` (or equivalently ``]a,b[``)
           The bounds ``a`` and ``b`` can be ``+/-Infinity``, ``Inf``,
           ``infinity``, ``inf`` or ``oo``. For *singular* coordinates,
           non-open intervals such as ``[a,b]`` and
           ``(a,b]`` (or equivalently ``]a,b]``) are allowed. Note that
           the interval declaration must not contain any space character.
        3. (optional) Indicator of the periodic character of the coordinate,
           either as ``period=T``, where ``T`` is the period, or, for manifolds
           over `\RR` only, as the keyword ``periodic`` (the value of the
           period is then deduced from the interval `I` declared in field 2;
           see the example below)
        4. (optional) The LaTeX spelling of the coordinate; if not provided
           the coordinate symbol given in the first field will be used.

        The order of fields 2 to 4 does not matter and each of them can
        be omitted. If it contains any LaTeX expression, the string
        ``coordinates`` must be declared with the prefix 'r' (for "raw") to
        allow for a proper treatment of the backslash character (see
        examples below). If no interval range, no period and no LaTeX spelling
        is to be set for any coordinate, the argument ``coordinates`` can be
        omitted when the shortcut operator ``<,>`` is used to declare the
        chart (see examples below).

        Additional restrictions on the coordinates can be set using the
        argument ``coord_restrictions``.

        A restriction can be any symbolic equality or inequality involving
        the coordinates, such as ``x > y`` or ``x^2 + y^2 != 0``. The items
        of the list (or set or frozenset) ``coord_restrictions`` are combined
        with the ``and`` operator; if some restrictions are to be combined with
        the ``or`` operator instead, they have to be passed as a tuple in some
        single item of the list (or set or frozenset) ``coord_restrictions``.
        For example::

          coord_restrictions=[x > y, (x != 0, y != 0), z^2 < x]

        means ``(x > y) and ((x != 0) or (y != 0)) and (z^2 < x)``.
        If the list ``coord_restrictions`` contains only one item, this
        item can be passed as such, i.e. writing ``x > y`` instead
        of the single element list ``[x > y]``.  If the chart variables have
        not been declared as variables yet, ``coord_restrictions`` must
        be ``lambda``-quoted.

        OUTPUT:

        - the created chart, as an instance of
          :class:`~sage.manifolds.chart.Chart` or one of its subclasses, like
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart` for
          differentiable manifolds over `\RR`.

        EXAMPLES:

        Chart on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X = M.chart('x y'); X
            Chart (M, (x, y))
            sage: X[0]
            x
            sage: X[1]
            y
            sage: X[:]
            (x, y)

        The declared coordinates are not known at the global level::

            sage: y
            Traceback (most recent call last):
            ...
            NameError: name 'y' is not defined

        They can be recovered by the operator ``[:]`` applied to the chart::

            sage: (x, y) = X[:]
            sage: y
            y
            sage: type(y)
            <class 'sage.symbolic.expression.Expression'>

        But a shorter way to proceed is to use the operator ``<,>`` in the
        left-hand side of the chart declaration (there is then no need to
        pass the string 'x y' to chart())::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart(); X
            Chart (M, (x, y))

        Indeed, the declared coordinates are then known at the global level::

            sage: y
            y
            sage: (x,y) == X[:]
            True

        Actually the instruction ``X.<x,y> = M.chart()`` is
        equivalent to the combination of the two instructions
        ``X = M.chart('x y')`` and ``(x,y) = X[:]``.

        As an example of coordinate ranges and LaTeX symbols passed via the
        string ``coordinates`` to ``chart()``, let us introduce polar
        coordinates::

            sage: U = M.open_subset('U', coord_def={X: x^2+y^2 != 0})
            sage: P.<r,ph> = U.chart(r'r:(0,+oo) ph:(0,2*pi):periodic:\phi'); P
            Chart (U, (r, ph))
            sage: P.coord_range()
            r: (0, +oo); ph: [0, 2*pi] (periodic)
            sage: latex(P)
            \left(U,(r, {\phi})\right)

        Using ``coord_restrictions``::

            sage: D = Manifold(2, 'D', structure='topological')
            sage: X.<x,y> = D.chart(coord_restrictions=lambda x,y: [x^2+y^2<1, x>0]); X
            Chart (D, (x, y))
            sage: X.valid_coordinates(0, 0)
            False
            sage: X.valid_coordinates(1/2, 0)
            True

        See the documentation of classes
        :class:`~sage.manifolds.chart.Chart` and
        :class:`~sage.manifolds.chart.RealChart` for more examples,
        especially regarding the coordinates ranges and restrictions.
        """
        if calc_method is None:
            calc_method = self._calculus_method
        return self._structure.chart(self, coordinates=coordinates,
                                     names=names, calc_method=calc_method,
                                     coord_restrictions=coord_restrictions)

    def is_open(self):
        """
        Return if ``self`` is an open set.

        In the present case (manifold or open subset of it), always
        return ``True``.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: M.is_open()
            True
        """
        return True

    def set_orientation(self, orientation):
        r"""
        Set the preferred orientation of ``self``.

        INPUT:

        - ``orientation`` -- a chart or a list of charts

        .. WARNING::

            It is the user's responsibility that the orientation set here
            is indeed an orientation. There is no check going on in the
            background. See :meth:`orientation` for the definition of an
            orientation.

        EXAMPLES:

        Set an orientation on a manifold::

            sage: M = Manifold(2, 'M', structure='top')
            sage: c_xy.<x,y> = M.chart(); c_uv.<u,v> = M.chart()
            sage: M.set_orientation(c_uv)
            sage: M.orientation()
            [Chart (M, (u, v))]

        Set an orientation in the non-trivial case::

            sage: M = Manifold(2, 'M', structure='top')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: M.declare_union(U, V)
            sage: c_xy.<x,y> = U.chart(); c_uv.<u,v> = V.chart()
            sage: M.set_orientation([c_xy, c_uv])
            sage: M.orientation()
            [Chart (U, (x, y)), Chart (V, (u, v))]
        """
        chart_type = self._structure.chart
        if isinstance(orientation, chart_type):
            orientation = [orientation]
        elif isinstance(orientation, (tuple, list)):
            orientation = list(orientation)
        else:
            raise TypeError("orientation must be a chart or a list/tuple of "
                            "charts")
        dom_union = None
        for c in orientation:
            if not isinstance(c, chart_type):
                raise ValueError("orientation must consist of charts")
            dom = c._domain
            if not dom.is_subset(self):
                raise ValueError("{} must be defined ".format(c) +
                                 "on a subset of {}".format(self))
            if dom_union is not None:
                dom_union = dom.union(dom_union)
            else:
                dom_union = dom
        if dom_union != self:
            raise ValueError("chart domains must cover {}".format(self))
        self._orientation = orientation

    def orientation(self):
        r"""
        Get the preferred orientation of ``self`` if available.

        An *orientation* of an `n`-dimensional topologial manifold is an
        atlas of charts whose transition maps are orientation preserving. A
        homeomorphism `f \colon U \to V` for open subsets `U, V \subset \RR^n`
        is called *orientation preserving* if for each `x \in U` the
        following map between singular homologies is the identity:

        .. MATH::

            H_n(\RR^n, \RR^n - 0; \ZZ) \cong H_n(U, U - x; \ZZ)
            \xrightarrow{f_*} H_n(V, V - f(x)) \cong H_n(\RR^n, \RR^n - 0; \ZZ)

        See `this link
        <http://www.map.mpim-bonn.mpg.de/Orientation_of_manifolds>`_
        for details.

        .. NOTE::

            Notice that for differentiable manifolds, the notion of
            orientability does not need homology theory at all. See
            :meth:`~sage.manifolds.differentiable.manifold.DifferentiableManifold.orientation`
            for details

        The trivial case corresponds to the manifold being covered by
        one chart. In that case, if no preferred orientation has been manually
        set before, one of those charts (usually the default chart) is
        set to the preferred orientation and returned here.

        EXAMPLES:

        If the manifold is covered by only one chart, it certainly admits an
        orientation::

            sage: M = Manifold(3, 'M', structure='top')
            sage: c.<x,y,z> = M.chart()
            sage: M.orientation()
            [Chart (M, (x, y, z))]

        Usually, an orientation cannot be obtained so easily::

            sage: M = Manifold(2, 'M', structure='top')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: M.declare_union(U, V)
            sage: c_xy.<x,y> = U.chart(); c_uv.<u,v> = V.chart()
            sage: M.orientation()
            []

        In that case, the orientation can be set by the user manually::

            sage: M.set_orientation([c_xy, c_uv])
            sage: M.orientation()
            [Chart (U, (x, y)), Chart (V, (u, v))]

        The orientation on submanifolds are inherited from the ambient
        manifold::

            sage: W = U.intersection(V, name='W')
            sage: W.orientation()
            [Chart (W, (x, y))]
        """
        if not self._orientation:
            # try to get an orientation from super domains:
            for sdom in self._supersets:
                sorient = sdom._orientation
                if sorient:
                    rst_orient = [c.restrict(self) for c in sorient]
                    # clear duplicated domains:
                    rst_orient = list(self._get_min_covering(rst_orient))
                    self._orientation = rst_orient
                    break
            else:
                # Trivial case:
                if self.is_manifestly_coordinate_domain():
                    # Try the default chart:
                    def_chart = self._def_chart
                    if def_chart is not None:
                        if def_chart.domain() is self:
                            self._orientation = [self._def_chart]
                    # Still no orientation? Choose arbitrary chart:
                    if not self._orientation:
                        for chart in self._covering_charts:
                            self._orientation = [chart]
                            break
        return list(self._orientation)

    def has_orientation(self):
        r"""
        Check whether ``self`` admits an obvious or by user set orientation.

        .. SEEALSO::

            Consult :meth:`orientation` for details about orientations.

        .. NOTE::

            Notice that if :meth:`has_orientation` returns ``False`` this
            does not necessarily mean that the manifold admits no orientation.
            It just means that the user has to set an orientation manually
            in that case, see :meth:`set_orientation`.

        EXAMPLES:

        The trivial case::

            sage: M = Manifold(3, 'M', structure='top')
            sage: c.<x,y,z> = M.chart()
            sage: M.has_orientation()
            True

        The non-trivial case::

            sage: M = Manifold(2, 'M', structure='top')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: M.declare_union(U, V)
            sage: c_xy.<x,y> = U.chart(); c_uv.<u,v> = V.chart()
            sage: M.has_orientation()
            False
            sage: M.set_orientation([c_xy, c_uv])
            sage: M.has_orientation()
            True
        """
        return bool(self.orientation())

    def _get_min_covering(self, object_list):
        r"""
        Helper method to return the minimal amount of objects necessary to
        cover the union of all their domains.

        INPUT:

        - list of objects having an `domain` method

        OUTPUT: set of objects

        TESTS::

            sage: M = Manifold(1, 'M', structure='top')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: c1.<x> = U.chart(); c2.<y> = V.chart()
            sage: c3.<z> = M.chart()
            sage: M._get_min_covering([c1, c2, c3])
            {Chart (M, (z,))}
        """
        min_obj_set = set()
        for obj in object_list:
            redund_obj_set = set()
            for oobj in min_obj_set:
                if obj.domain().is_subset(oobj.domain()):
                    break
                elif oobj.domain().is_subset(obj.domain()):
                    redund_obj_set.add(oobj)
            else:
                min_obj_set.add(obj)
            min_obj_set.difference_update(redund_obj_set)
        return min_obj_set

    def vector_bundle(self, rank, name, field='real', latex_name=None):
        r"""
        Return a topological vector bundle over the given field with given rank
        over this topological manifold.

        INPUT:

        - ``rank`` -- rank of the vector bundle
        - ``name`` -- name given to the total space
        - ``field`` -- (default: ``'real'``) topological field giving the
          vector space structure to the fibers
        - ``latex_name`` -- (optional) LaTeX name for the total space

        OUTPUT:

        - a topological vector bundle as an instance of
          :class:`~sage.manifolds.vector_bundle.TopologicalVectorBundle`

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='top')
            sage: M.vector_bundle(2, 'E')
            Topological real vector bundle E -> M of rank 2 over the base space
             2-dimensional topological manifold M
        """
        from sage.manifolds.vector_bundle import TopologicalVectorBundle
        return TopologicalVectorBundle(rank, name, self, field=field,
                                       latex_name=latex_name)

    def scalar_field_algebra(self):
        r"""
        Return the algebra of scalar fields defined the manifold.

        See :class:`~sage.manifolds.scalarfield_algebra.ScalarFieldAlgebra`
        for a complete documentation.

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.scalarfield_algebra.ScalarFieldAlgebra`
          representing the algebra `C^0(U)` of all scalar fields defined
          on `U` = ``self``

        EXAMPLES:

        Scalar algebra of a 3-dimensional open subset::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: U = M.open_subset('U')
            sage: CU = U.scalar_field_algebra() ; CU
            Algebra of scalar fields on the Open subset U of the 3-dimensional topological manifold M
            sage: CU.category()
            Join of Category of commutative algebras over Symbolic Ring and Category of homsets of topological spaces
            sage: CU.zero()
            Scalar field zero on the Open subset U of the 3-dimensional topological manifold M

        The output is cached::

            sage: U.scalar_field_algebra() is CU
            True
        """
        return self._scalar_field_algebra

    def scalar_field(
        self, coord_expression=None, chart=None, name=None, latex_name=None
    ) -> ScalarField:
        r"""
        Define a scalar field on the manifold.

        See :class:`~sage.manifolds.scalarfield.ScalarField` (or
        :class:`~sage.manifolds.differentiable.scalarfield.DiffScalarField`
        if the manifold is differentiable) for a complete documentation.

        INPUT:

        - ``coord_expression`` -- (default: ``None``) coordinate expression(s)
          of the scalar field; this can be either

          * a single coordinate expression; if the argument ``chart`` is
            ``'all'``, this expression is set to all the charts defined
            on the open set; otherwise, the expression is set in the
            specific chart provided by the argument ``chart``
          * a dictionary of coordinate expressions, with the charts as keys

        - ``chart`` -- (default: ``None``) chart defining the coordinates
          used in ``coord_expression`` when the latter is a single
          coordinate expression; if ``None``, the default chart of the
          open set is assumed; if ``chart=='all'``, ``coord_expression`` is
          assumed to be independent of the chart (constant scalar field)

        - ``name`` -- (default: ``None``) name given to the scalar field

        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          scalar field; if ``None``, the LaTeX symbol is set to ``name``

        If ``coord_expression`` is ``None`` or does not fully specified the
        scalar field, other coordinate expressions can be added subsequently
        by means of the methods
        :meth:`~sage.manifolds.scalarfield.ScalarField.add_expr`,
        :meth:`~sage.manifolds.scalarfield.ScalarField.add_expr_by_continuation`,
        or :meth:`~sage.manifolds.scalarfield.ScalarField.set_expr`

        OUTPUT:

        - instance of :class:`~sage.manifolds.scalarfield.ScalarField`
          (or of the subclass
          :class:`~sage.manifolds.differentiable.scalarfield.DiffScalarField`
          if the manifold is differentiable) representing the defined scalar
          field

        EXAMPLES:

        A scalar field defined by its coordinate expression in the open
        set's default chart::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: U = M.open_subset('U')
            sage: c_xyz.<x,y,z> = U.chart()
            sage: f = U.scalar_field(sin(x)*cos(y) + z, name='F'); f
            Scalar field F on the Open subset U of the 3-dimensional topological manifold M
            sage: f.display()
            F: U → ℝ
               (x, y, z) ↦ cos(y)*sin(x) + z
            sage: f.parent()
            Algebra of scalar fields on the Open subset U of the 3-dimensional topological manifold M
            sage: f in U.scalar_field_algebra()
            True

        Equivalent definition with the chart specified::

            sage: f = U.scalar_field(sin(x)*cos(y) + z, chart=c_xyz, name='F')
            sage: f.display()
            F: U → ℝ
               (x, y, z) ↦ cos(y)*sin(x) + z

        Equivalent definition with a dictionary of coordinate expression(s)::

            sage: f = U.scalar_field({c_xyz: sin(x)*cos(y) + z}, name='F')
            sage: f.display()
            F: U → ℝ
               (x, y, z) ↦ cos(y)*sin(x) + z

        See the documentation of class
        :class:`~sage.manifolds.scalarfield.ScalarField` for more
        examples.

        .. SEEALSO::

            :meth:`constant_scalar_field`, :meth:`zero_scalar_field`,
            :meth:`one_scalar_field`
        """
        if isinstance(coord_expression, dict):
            # check validity of entry
            for chart in coord_expression:
                if not chart.domain().is_subset(self):
                    raise ValueError("the {} is not defined ".format(chart) +
                                     "on some subset of the " + str(self))
        alg = self.scalar_field_algebra()
        return alg.element_class(alg, coord_expression=coord_expression,
                                 name=name, latex_name=latex_name, chart=chart)

    def constant_scalar_field(self, value, name=None, latex_name=None):
        r"""
        Define a constant scalar field on the manifold.

        INPUT:

        - ``value`` -- constant value of the scalar field, either a numerical
          value or a symbolic expression not involving any chart coordinates
        - ``name`` -- (default: ``None``) name given to the scalar field
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          scalar field; if ``None``, the LaTeX symbol is set to ``name``

        OUTPUT:

        - instance of :class:`~sage.manifolds.scalarfield.ScalarField`
          representing the scalar field whose constant value is ``value``

        EXAMPLES:

        A constant scalar field on the 2-sphere::

            sage: M = Manifold(2, 'M', structure='topological') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            ....:                                intersection_name='W',
            ....:                                restrictions1= x^2+y^2!=0,
            ....:                                restrictions2= u^2+v^2!=0)
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: f = M.constant_scalar_field(-1) ; f
            Scalar field on the 2-dimensional topological manifold M
            sage: f.display()
            M → ℝ
            on U: (x, y) ↦ -1
            on V: (u, v) ↦ -1

        We have::

            sage: f.restrict(U) == U.constant_scalar_field(-1)
            True
            sage: M.constant_scalar_field(0) is M.zero_scalar_field()
            True

        .. SEEALSO::

            :meth:`zero_scalar_field`, :meth:`one_scalar_field`
        """
        if value == 0:
            return self.zero_scalar_field()
        alg = self.scalar_field_algebra()
        return alg.element_class(alg, coord_expression=value, name=name,
                                 latex_name=latex_name, chart='all')

    def zero_scalar_field(self):
        r"""
        Return the zero scalar field defined on ``self``.

        OUTPUT:

        - a :class:`~sage.manifolds.scalarfield.ScalarField`
          representing the constant scalar field with value `0`

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.zero_scalar_field() ; f
            Scalar field zero on the 2-dimensional topological manifold M
            sage: f.display()
            zero: M → ℝ
               (x, y) ↦ 0
            sage: f.parent()
            Algebra of scalar fields on the 2-dimensional topological manifold M
            sage: f is M.scalar_field_algebra().zero()
            True
        """
        return self._zero_scalar_field

    def one_scalar_field(self):
        r"""
        Return the constant scalar field with value the unit element
        of the base field of ``self``.

        OUTPUT:

        - a :class:`~sage.manifolds.scalarfield.ScalarField` representing
          the constant scalar field with value the unit element
          of the base field of ``self``

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.one_scalar_field(); f
            Scalar field 1 on the 2-dimensional topological manifold M
            sage: f.display()
            1: M → ℝ
               (x, y) ↦ 1
            sage: f.parent()
            Algebra of scalar fields on the 2-dimensional topological manifold M
            sage: f is M.scalar_field_algebra().one()
            True
        """
        return self._one_scalar_field

    class options(GlobalOptions):
        r"""
        Set and displays the options for manifolds. If no parameters
        are set, then the function returns a copy of the options dictionary.

        The ``options`` to manifolds can be accessed as the method
        :obj:`Manifold.options`.

        @OPTIONS@

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: g = function('g')(x, y)

        For coordinate functions, the display is more "textbook" like::

            sage: f = X.function(diff(g, x) + diff(g, y))
            sage: f
            d(g)/dx + d(g)/dy
            sage: latex(f)
            \frac{\partial\,g}{\partial x} + \frac{\partial\,g}{\partial y}

        One can switch to Pynac notation by changing ``textbook_output``
        to ``False``::

            sage: Manifold.options.textbook_output=False
            sage: f
            diff(g(x, y), x) + diff(g(x, y), y)
            sage: latex(f)
            \frac{\partial}{\partial x}g\left(x, y\right)
             + \frac{\partial}{\partial y}g\left(x, y\right)
            sage: Manifold.options._reset()

        If there is a clear understanding that `u` and `v` are functions of
        `(x,y)`, the explicit mention of the latter can be cumbersome in lengthy
        tensor expressions::

            sage: f = X.function(function('u')(x, y) * function('v')(x, y))
            sage: f
            u(x, y)*v(x, y)

        We can switch it off by::

            sage: M.options.omit_function_arguments=True
            sage: f
            u*v
            sage: M.options._reset()
        """
        NAME = 'manifolds'
        module = 'sage.manifolds.manifold'
        option_class = 'TopologicalManifold'
        textbook_output = dict(default=True,
                             description='textbook-like output instead of the Pynac output for derivatives',
                             checker=lambda x: isinstance(x, bool))
        omit_function_arguments = dict(default=False,
                                     description='Determine whether the arguments of symbolic functions are printed',
                                     checker=lambda x: isinstance(x, bool))

    def _Hom_(self, other, category=None):
        r"""
        Construct the set of morphisms (i.e. continuous maps)
        ``self`` to ``other``.

        INPUT:

        - ``other`` -- an open subset of some topological manifold over the
          same field as ``self``
        - ``category`` -- (default: ``None``) not used here (to ensure
          compatibility with generic hook ``_Hom_``)

        OUTPUT:

        - the homset `\mathrm{Hom}(U,V)`, where `U` is ``self``
          and `V` is ``other``

        .. SEEALSO::

            For more documentation, see
            :class:`~sage.manifolds.manifold_homset.TopologicalManifoldHomset`.

        EXAMPLES::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: N = Manifold(3, 'N', structure='topological')
            sage: H = M._Hom_(N); H
            Set of Morphisms from 2-dimensional topological manifold M to
             3-dimensional topological manifold N in Category of manifolds over
             Real Field with 53 bits of precision

        The result is cached::

            sage: H is Hom(M, N)
            True
        """
        return self._structure.homset(self, other)

    def continuous_map(self, codomain, coord_functions=None, chart1=None,
                       chart2=None, name=None, latex_name=None):
        r"""
        Define a continuous map from ``self`` to ``codomain``.

        INPUT:

        - ``codomain`` -- :class:`TopologicalManifold`; the map's codomain
        - ``coord_functions`` -- (default: ``None``) if not ``None``,
          must be either

          - (i) a dictionary of the coordinate expressions (as lists
            (or tuples) of the coordinates of the image expressed in
            terms of the coordinates of the considered point) with the
            pairs of charts ``(chart1, chart2)`` as keys (``chart1`` being
            a chart on ``self`` and ``chart2`` a chart on ``codomain``);
          - (ii) a single coordinate expression in a given pair of charts, the
            latter being provided by the arguments ``chart1`` and ``chart2``;

          in both cases, if the dimension of the codomain is `1`, a single
          coordinate expression can be passed instead of a tuple with
          a single element
        - ``chart1`` -- (default: ``None``; used only in case (ii) above)
          chart on ``self`` defining the start coordinates involved in
          ``coord_functions`` for case (ii); if ``None``, the coordinates
          are assumed to refer to the default chart of ``self``
        - ``chart2`` -- (default: ``None``; used only in case (ii) above)
          chart on ``codomain`` defining the target coordinates involved in
          ``coord_functions`` for case (ii); if ``None``, the coordinates
          are assumed to refer to the default chart of ``codomain``
        - ``name`` -- (default: ``None``) name given to the continuous map
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          continuous map; if ``None``, the LaTeX symbol is set to ``name``

        OUTPUT:

        - the continuous map as an instance of
          :class:`~sage.manifolds.continuous_map.ContinuousMap`

        EXAMPLES:

        A continuous map between an open subset of `S^2` covered by regular
        spherical coordinates and `\RR^3`::

            sage: M = Manifold(2, 'S^2', structure='topological')
            sage: U = M.open_subset('U')
            sage: c_spher.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
            sage: N = Manifold(3, 'R^3', latex_name=r'\RR^3', structure='topological')
            sage: c_cart.<x,y,z> = N.chart()  # Cartesian coord. on R^3
            sage: Phi = U.continuous_map(N, (sin(th)*cos(ph), sin(th)*sin(ph), cos(th)),
            ....:                        name='Phi', latex_name=r'\Phi')
            sage: Phi
            Continuous map Phi from the Open subset U of the 2-dimensional topological manifold S^2 to the 3-dimensional topological manifold R^3

        The same definition, but with a dictionary with pairs of charts as
        keys (case (i) above)::

            sage: Phi1 = U.continuous_map(N,
            ....:        {(c_spher, c_cart): (sin(th)*cos(ph), sin(th)*sin(ph), cos(th))},
            ....:        name='Phi', latex_name=r'\Phi')
            sage: Phi1 == Phi
            True

        The continuous map acting on a point::

            sage: p = U.point((pi/2, pi)) ; p
            Point on the 2-dimensional topological manifold S^2
            sage: Phi(p)
            Point on the 3-dimensional topological manifold R^3
            sage: Phi(p).coord(c_cart)
            (-1, 0, 0)
            sage: Phi1(p) == Phi(p)
            True

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: M.continuous_map(ZZ)
            Traceback (most recent call last):
            ...
            ValueError: Integer Ring is not a manifold
             over Real Field with 53 bits of precision

        .. SEEALSO::

            See :class:`~sage.manifolds.continuous_map.ContinuousMap`
            for the complete documentation and more examples.

        .. TODO::

            Allow the construction of continuous maps from ``self`` to the
            base field (considered as a trivial 1-dimensional manifold).
        """
        if (not isinstance(codomain, TopologicalManifold)
            or codomain.base_field() != self.base_field()):
            raise ValueError("{} is not a manifold over {}".format(codomain, self.base_field()))
        homset = Hom(self, codomain)
        if coord_functions is None:
            coord_functions = {}
        if not isinstance(coord_functions, dict):
            # Turn coord_functions into a dictionary:
            if chart1 is None:
                chart1 = self._def_chart
            elif chart1 not in self._atlas:
                raise ValueError("{} is not a chart ".format(chart1) +
                                 "defined on the {}".format(self))
            if chart2 is None:
                chart2 = codomain._def_chart
            elif chart2 not in codomain._atlas:
                raise ValueError("{} is not a chart ".format(chart2) +
                                 " defined on the {}".format(codomain))
            coord_functions = {(chart1, chart2): coord_functions}
        return homset(coord_functions, name=name, latex_name=latex_name)

    def homeomorphism(self, codomain, coord_functions=None, chart1=None,
                       chart2=None, name=None, latex_name=None):
        r"""
        Define a homeomorphism between the current manifold and another one.

        See :class:`~sage.manifolds.continuous_map.ContinuousMap` for a
        complete documentation.

        INPUT:

        - ``codomain`` -- :class:`TopologicalManifold`; codomain of
          the homeomorphism
        - ``coord_functions`` -- (default: ``None``) if not ``None``,
          must be either

          - (i) a dictionary of the coordinate expressions (as lists
            (or tuples) of the coordinates of the image expressed in
            terms of the coordinates of the considered point) with the
            pairs of charts ``(chart1, chart2)`` as keys (``chart1`` being
            a chart on ``self`` and ``chart2`` a chart on ``codomain``);
          - (ii) a single coordinate expression in a given pair of charts, the
            latter being provided by the arguments ``chart1`` and ``chart2``;

          in both cases, if the dimension of the codomain is `1`, a single
          coordinate expression can be passed instead of a tuple with
          a single element
        - ``chart1`` -- (default: ``None``; used only in case (ii) above)
          chart on ``self`` defining the start coordinates involved in
          ``coord_functions`` for case (ii); if ``None``, the coordinates
          are assumed to refer to the default chart of ``self``
        - ``chart2`` -- (default: ``None``; used only in case (ii) above)
          chart on ``codomain`` defining the target coordinates involved in
          ``coord_functions`` for case (ii); if ``None``, the coordinates
          are assumed to refer to the default chart of ``codomain``
        - ``name`` -- (default: ``None``) name given to the homeomorphism
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          homeomorphism; if ``None``, the LaTeX symbol is set to ``name``

        OUTPUT:

        - the homeomorphism, as an instance of
          :class:`~sage.manifolds.continuous_map.ContinuousMap`

        EXAMPLES:

        Homeomorphism between the open unit disk in `\RR^2` and `\RR^2`::

            sage: forget()  # for doctests only
            sage: M = Manifold(2, 'M', structure='topological')  # the open unit disk
            sage: c_xy.<x,y> = M.chart('x:(-1,1) y:(-1,1)', coord_restrictions=lambda x,y: x^2+y^2<1)
            ....:    # Cartesian coord on M
            sage: N = Manifold(2, 'N', structure='topological')  # R^2
            sage: c_XY.<X,Y> = N.chart()  # canonical coordinates on R^2
            sage: Phi = M.homeomorphism(N, [x/sqrt(1-x^2-y^2), y/sqrt(1-x^2-y^2)],
            ....:                       name='Phi', latex_name=r'\Phi')
            sage: Phi
            Homeomorphism Phi from the 2-dimensional topological manifold M to
             the 2-dimensional topological manifold N
            sage: Phi.display()
            Phi: M → N
               (x, y) ↦ (X, Y) = (x/sqrt(-x^2 - y^2 + 1), y/sqrt(-x^2 - y^2 + 1))

        The inverse homeomorphism::

            sage: Phi^(-1)
            Homeomorphism Phi^(-1) from the 2-dimensional topological
             manifold N to the 2-dimensional topological manifold M
            sage: (Phi^(-1)).display()
            Phi^(-1): N → M
               (X, Y) ↦ (x, y) = (X/sqrt(X^2 + Y^2 + 1), Y/sqrt(X^2 + Y^2 + 1))

        See the documentation of
        :class:`~sage.manifolds.continuous_map.ContinuousMap` for more
        examples.
        """
        homset = Hom(self, codomain)
        if coord_functions is None:
            coord_functions = {}
        if not isinstance(coord_functions, dict):
            # Turn coord_functions into a dictionary:
            if chart1 is None:
                chart1 = self._def_chart
            elif chart1 not in self._atlas:
                raise ValueError("{} is not a chart ".format(chart1) +
                                 "defined on the {}".format(self))
            if chart2 is None:
                chart2 = codomain._def_chart
            elif chart2 not in codomain._atlas:
                raise ValueError("{} is not a chart ".format(chart2) +
                                 " defined on the {}".format(codomain))
            coord_functions = {(chart1, chart2): coord_functions}
        return homset(coord_functions, name=name, latex_name=latex_name,
                      is_isomorphism=True)

    @overload
    def identity_map(self: TopologicalManifold) -> ContinuousMap: ...
    @overload
    def identity_map(self: DifferentiableManifold) -> DiffMap: ...
    @cached_method
    def identity_map(self):
        r"""
        Identity map of ``self``.

        The identity map of a topological manifold `M` is the trivial
        homeomorphism:

        .. MATH::

            \begin{array}{cccc}
            \mathrm{Id}_M: & M & \longrightarrow & M \\
                & p & \longmapsto & p
            \end{array}

        OUTPUT:

        - the identity map as an instance of
          :class:`~sage.manifolds.continuous_map.ContinuousMap`

        EXAMPLES:

        Identity map of a complex manifold::

            sage: M = Manifold(2, 'M', structure='topological', field='complex')
            sage: X.<x,y> = M.chart()
            sage: id = M.identity_map(); id
            Identity map Id_M of the Complex 2-dimensional topological manifold M
            sage: id.parent()
            Set of Morphisms from Complex 2-dimensional topological manifold M
             to Complex 2-dimensional topological manifold M in Category of
             manifolds over Complex Field with 53 bits of precision
            sage: id.display()
            Id_M: M → M
               (x, y) ↦ (x, y)

        The identity map acting on a point::

            sage: p = M((1+I, 3-I), name='p'); p
            Point p on the Complex 2-dimensional topological manifold M
            sage: id(p)
            Point p on the Complex 2-dimensional topological manifold M
            sage: id(p) == p
            True

        .. SEEALSO::

            See :class:`~sage.manifolds.continuous_map.ContinuousMap`
            for the complete documentation.
        """
        return Hom(self, self).one()

    def set_calculus_method(self, method):
        r"""
        Set the calculus method to be used for coordinate computations on this
        manifold.

        The provided method is transmitted to all coordinate charts defined on
        the manifold.

        INPUT:

        - ``method`` -- string specifying the method to be used for
          coordinate computations on this manifold; one of

          - ``'SR'``: Sage's default symbolic engine (Symbolic Ring)
          - ``'sympy'``: SymPy

        EXAMPLES:

        Let us consider a scalar field ``f`` on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field(x^2 + cos(y)*sin(x), name='F')

        By default, the coordinate expression of ``f`` returned by
        :meth:`~sage.manifolds.scalarfield.ScalarField.expr` is a Sage's
        symbolic expression::

            sage: f.expr()
            x^2 + cos(y)*sin(x)
            sage: type(f.expr())
            <class 'sage.symbolic.expression.Expression'>
            sage: parent(f.expr())
            Symbolic Ring
            sage: f.display()
            F: M → ℝ
               (x, y) ↦ x^2 + cos(y)*sin(x)

        If we change the calculus method to SymPy, it becomes a SymPy object
        instead::

            sage: M.set_calculus_method('sympy')
            sage: f.expr()
            x**2 + sin(x)*cos(y)
            sage: type(f.expr())
            <class 'sympy.core.add.Add'>
            sage: parent(f.expr())
            <class 'sympy.core.add.Add'>
            sage: f.display()
            F: M → ℝ
               (x, y) ↦ x**2 + sin(x)*cos(y)

        Back to the Symbolic Ring::

            sage: M.set_calculus_method('SR')
            sage: f.display()
            F: M → ℝ
               (x, y) ↦ x^2 + cos(y)*sin(x)

        The calculus method chosen via ``set_calculus_method()`` applies to any
        chart defined subsequently on the manifold::

            sage: M.set_calculus_method('sympy')
            sage: Y.<u,v> = M.chart()  # a new chart
            sage: Y.calculus_method()
            Available calculus methods (* = current):
             - SR (default)
             - sympy (*)

        .. SEEALSO::

            :meth:`~sage.manifolds.chart.Chart.calculus_method` for a
            control of the calculus method chart by chart
        """
        self._calculus_method = method
        for chart in self._atlas:
            chart.calculus_method().set(method)

    def set_simplify_function(self, simplifying_func, method=None):
        r"""
        Set the simplifying function associated to a given coordinate
        calculus method in all the charts defined on ``self``.

        INPUT:

        - ``simplifying_func`` -- either the string ``'default'`` for restoring
          the default simplifying function or a function ``f`` of a single
          argument ``expr`` such that ``f(expr)`` returns an object of the same
          type as ``expr`` (hopefully the simplified version of ``expr``), this
          type being

          - :class:`~sage.symbolic.expression.Expression` if ``method`` = ``'SR'``
          - a SymPy type if ``method`` = ``'sympy'``

        - ``method`` -- (default: ``None``) string defining the calculus method
          for which ``simplifying_func`` is provided; must be one of

          - ``'SR'``: Sage's default symbolic engine (Symbolic Ring)
          - ``'sympy'``: SymPy
          - ``None``: the currently active calculus method on each chart is
            assumed

        .. SEEALSO::

            :meth:`~sage.manifolds.chart.Chart.calculus_method`
            and :meth:`sage.manifolds.calculus_method.CalculusMethod.simplify`
            for a control of the calculus method chart by chart

        EXAMPLES:

        Les us add two scalar fields on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: f = M.scalar_field((x+y)^2 + cos(x)^2)
            sage: g = M.scalar_field(-x^2-2*x*y-y^2 + sin(x)^2)
            sage: f.expr()
            (x + y)^2 + cos(x)^2
            sage: g.expr()
            -x^2 - 2*x*y - y^2 + sin(x)^2
            sage: s = f + g

        The outcome is automatically simplified::

            sage: s.expr()
            1

        The simplification is performed thanks to the default simplifying
        function on chart ``X``, which is
        :func:`~sage.manifolds.utilities.simplify_chain_real` in the present
        case (real manifold and ``SR`` calculus)::

            sage: X.calculus_method().simplify_function() is \
            ....: sage.manifolds.utilities.simplify_chain_real
            True

        Let us change it to the generic Sage function
        :func:`~sage.calculus.functional.simplify`::

            sage: M.set_simplify_function(simplify)
            sage: X.calculus_method().simplify_function() is simplify
            True

        :func:`~sage.calculus.functional.simplify` is faster, but it does not
        do much::

            sage: s = f + g
            sage: s.expr()
            (x + y)^2 - x^2 - 2*x*y - y^2 + cos(x)^2 + sin(x)^2

        We can replaced it by any user defined function, for instance::

            sage: def simpl_trig(a):
            ....:     return a.simplify_trig()
            sage: M.set_simplify_function(simpl_trig)
            sage: s = f + g
            sage: s.expr()
            1

        The default simplifying function is restored via::

            sage: M.set_simplify_function('default')

        Then we are back to::

            sage: X.calculus_method().simplify_function() is \
            ....: sage.manifolds.utilities.simplify_chain_real
            True

        Thanks to the argument ``method``, one can specify a simplifying
        function for a calculus method distinct from the current one. For
        instance, let us define a simplifying function for SymPy (note that
        ``trigsimp()`` is a SymPy method only)::

            sage: def simpl_trig_sympy(a):
            ....:     return a.trigsimp()
            sage: M.set_simplify_function(simpl_trig_sympy, method='sympy')

        Then, it becomes active as soon as we change the calculus engine to
        SymPy::

            sage: M.set_calculus_method('sympy')
            sage: X.calculus_method().simplify_function() is simpl_trig_sympy
            True

        We have then::

            sage: s = f + g
            sage: s.expr()
            1
            sage: type(s.expr())
            <class 'sympy.core.numbers.One'>
        """
        for chart in self._atlas:
            chart.calculus_method().set_simplify_function(simplifying_func,
                                                          method=method)


###########################################################
# Constructor function

_manifold_id = Integer(0)


def Manifold(
    dim: int,
    name: Optional[str],
    latex_name: Optional[str] = None,
    field: str = "real",
    structure: Optional[str] = None,
    start_index: int = 0,
    **extra_kwds,
) -> Union[TopologicalManifold, DifferentiableManifold]:
    r"""
    Construct a manifold of a given type over a topological field.

    Given a topological field `K` (in most applications, `K = \RR` or
    `K = \CC`) and a nonnegative integer `n`, a *topological manifold of
    dimension* `n` *over K* is a topological space `M` such that

    - `M` is a Hausdorff space,
    - `M` is second countable, and
    - every point in `M` has a neighborhood homeomorphic to `K^n`.

    A *real manifold* is a manifold over `\RR`. A *differentiable* (resp.
    *smooth*, resp. *analytic*) *manifold* is a manifold such that all
    transition maps are *differentiable* (resp. *smooth*, resp. *analytic*). A
    *pseudo-Riemannian manifold* is a real differentiable manifold equipped
    with a metric tensor `g` (i.e. a field of non-degenerate symmetric bilinear
    forms), with the two subcases of *Riemannian manifold* (`g`
    positive-definite) and *Lorentzian manifold* (`g` has signature `n-2` or
    `2-n`).

    INPUT:

    - ``dim`` -- positive integer; dimension of the manifold
    - ``name`` -- string; name (symbol) given to the manifold
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the manifold; if none are provided, it is set to ``name``
    - ``field`` -- (default: ``'real'``) field `K` on which the
      manifold is defined; allowed values are

      - ``'real'`` or an object of type ``RealField`` (e.g. ``RR``) for a
        manifold over `\RR`
      - ``'complex'`` or an object of type ``ComplexField`` (e.g. ``CC``)
        for a manifold over `\CC`
      - an object in the category of topological fields (see
        :class:`~sage.categories.fields.Fields` and
        :class:`~sage.categories.topological_spaces.TopologicalSpaces`)
        for other types of manifolds

    - ``structure`` -- (default: ``'smooth'``) to specify the structure or
      type of manifold; allowed values are

      - ``'topological'`` or ``'top'`` for a topological manifold
      - ``'differentiable'`` or ``'diff'`` for a differentiable manifold
      - ``'smooth'`` for a smooth manifold
      - ``'analytic'`` for an analytic manifold
      - ``'pseudo-Riemannian'`` for a real differentiable manifold equipped
        with a pseudo-Riemannian metric; the signature is specified via the
        keyword argument ``signature`` (see below)
      - ``'Riemannian'`` for a real differentiable manifold equipped with a
        Riemannian (i.e. positive definite) metric
      - ``'Lorentzian'`` for a real differentiable manifold equipped with a
        Lorentzian metric; the signature convention is specified by the
        keyword argument ``signature='positive'`` (default) or ``'negative'``

    - ``start_index`` -- (default: 0) integer; lower value of the range of
      indices used for "indexed objects" on the manifold, e.g. coordinates
      in a chart
    - ``extra_kwds`` -- keywords meaningful only for some specific types
      of manifolds:

      - ``diff_degree`` -- (only for differentiable manifolds; default:
        ``infinity``): the degree of differentiability
      - ``ambient`` -- (only to construct a submanifold): the ambient manifold
      - ``metric_name`` -- (only for pseudo-Riemannian manifolds; default:
        ``'g'``) string; name (symbol) given to the metric
      - ``metric_latex_name`` -- (only for pseudo-Riemannian manifolds;
        default: ``None``) string; LaTeX symbol to denote the metric; if none
        is provided, the symbol is set to ``metric_name``
      - ``signature`` -- (only for pseudo-Riemannian manifolds; default:
        ``None``) signature `S` of the metric as a single integer:
        `S = n_+ - n_-`, where `n_+` (resp. `n_-`) is the number of positive
        terms (resp. negative terms) in any diagonal writing of the
        metric components; if ``signature`` is not provided, `S` is set to the
        manifold's dimension (Riemannian signature); for Lorentzian manifolds
        the values ``signature='positive'`` (default) or
        ``signature='negative'`` are allowed to indicate the chosen signature
        convention.

    OUTPUT:

    - a manifold of the specified type, as an instance of
      :class:`~sage.manifolds.manifold.TopologicalManifold` or one of its
      subclasses
      :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`
      or
      :class:`~sage.manifolds.differentiable.pseudo_riemannian.PseudoRiemannianManifold`,
      or, if the keyword ``ambient`` is used, one of the subclasses
      :class:`~sage.manifolds.topological_submanifold.TopologicalSubmanifold`,
      :class:`~sage.manifolds.differentiable.differentiable_submanifold.DifferentiableSubmanifold`,
      or
      :class:`~sage.manifolds.differentiable.pseudo_riemannian_submanifold.PseudoRiemannianSubmanifold`.

    EXAMPLES:

    A 3-dimensional real topological manifold::

        sage: M = Manifold(3, 'M', structure='topological'); M
        3-dimensional topological manifold M

    Given the default value of the parameter ``field``, the above is
    equivalent to::

        sage: M = Manifold(3, 'M', structure='topological', field='real'); M
        3-dimensional topological manifold M

    A complex topological manifold::

        sage: M = Manifold(3, 'M', structure='topological', field='complex'); M
        Complex 3-dimensional topological manifold M

    A topological manifold over `\QQ`::

        sage: M = Manifold(3, 'M', structure='topological', field=QQ); M
        3-dimensional topological manifold M over the Rational Field

    A 3-dimensional real differentiable manifold of class `C^4`::

        sage: M = Manifold(3, 'M', field='real', structure='differentiable',
        ....:              diff_degree=4); M
        3-dimensional differentiable manifold M

    Since the default value of the parameter ``field`` is ``'real'``, the above
    is equivalent to::

        sage: M = Manifold(3, 'M', structure='differentiable', diff_degree=4)
        sage: M
        3-dimensional differentiable manifold M
        sage: M.base_field_type()
        'real'

    A 3-dimensional real smooth manifold::

        sage: M = Manifold(3, 'M', structure='differentiable', diff_degree=+oo)
        sage: M
        3-dimensional differentiable manifold M

    Instead of ``structure='differentiable', diff_degree=+oo``, it suffices to
    use ``structure='smooth'`` to get the same result::

        sage: M = Manifold(3, 'M', structure='smooth'); M
        3-dimensional differentiable manifold M
        sage: M.diff_degree()
        +Infinity

    Actually, since ``'smooth'`` is the default value of the parameter
    ``structure``, the creation of a real smooth manifold can be shortened to::

        sage: M = Manifold(3, 'M'); M
        3-dimensional differentiable manifold M
        sage: M.diff_degree()
        +Infinity

    Other parameters can change the default of the parameter ``structure``::

        sage: M = Manifold(3, 'M', diff_degree=0); M
        3-dimensional topological manifold M
        sage: M = Manifold(3, 'M', diff_degree=2); M
        3-dimensional differentiable manifold M
        sage: M = Manifold(3, 'M', metric_name='g'); M
        3-dimensional Riemannian manifold M

    For a complex smooth manifold, we have to set the parameter ``field``::

        sage: M = Manifold(3, 'M', field='complex'); M
        3-dimensional complex manifold M
        sage: M.diff_degree()
        +Infinity

    Submanifolds are constructed by means of the keyword ``ambient``::

        sage: N = Manifold(2, 'N', field='complex', ambient=M); N
        2-dimensional differentiable submanifold N immersed in the
         3-dimensional complex manifold M

    The immersion `N\to M` has to be specified in a second stage, via the
    method
    :meth:`~sage.manifolds.topological_submanifold.TopologicalSubmanifold.set_immersion`
    or
    :meth:`~sage.manifolds.topological_submanifold.TopologicalSubmanifold.set_embedding`.

    For more detailed examples, see the documentation of
    :class:`~sage.manifolds.manifold.TopologicalManifold`,
    :class:`~sage.manifolds.differentiable.manifold.DifferentiableManifold`
    and
    :class:`~sage.manifolds.differentiable.pseudo_riemannian.PseudoRiemannianManifold`,
    or the documentation of
    :class:`~sage.manifolds.topological_submanifold.TopologicalSubmanifold`,
    :class:`~sage.manifolds.differentiable.differentiable_submanifold.DifferentiableSubmanifold`
    and
    :class:`~sage.manifolds.differentiable.pseudo_riemannian_submanifold.PseudoRiemannianSubmanifold`
    for submanifolds.

    .. RUBRIC:: Uniqueness of manifold objects

    Suppose we construct a manifold named `M`::

        sage: M = Manifold(2, 'M', structure='topological')
        sage: X.<x,y> = M.chart()

    At some point, we change our mind and would like to restart with a new
    manifold, using the same name `M` and keeping the previous manifold for
    reference::

        sage: M_old = M  # for reference
        sage: M = Manifold(2, 'M', structure='topological')

    This results in a brand new object::

        sage: M.atlas()
        []

    The object ``M_old`` is intact::

        sage: M_old.atlas()
        [Chart (M, (x, y))]

    Both objects have the same display::

        sage: M
        2-dimensional topological manifold M
        sage: M_old
        2-dimensional topological manifold M

    but they are different::

        sage: M != M_old
        True

    Let us introduce a chart on ``M``, using the same coordinate symbols
    as for ``M_old``::

        sage: X.<x,y> = M.chart()

    The charts are displayed in the same way::

        sage: M.atlas()
        [Chart (M, (x, y))]
        sage: M_old.atlas()
        [Chart (M, (x, y))]

    but they are actually different::

        sage: M.atlas()[0] != M_old.atlas()[0]
        True

    Moreover, the two manifolds ``M`` and ``M_old`` are still considered
    distinct::

        sage: M != M_old
        True

    This reflects the fact that the equality of manifold objects holds only
    for identical objects, i.e. one has ``M1 == M2`` if, and only if,
    ``M1 is M2``. Actually, the manifold classes inherit from
    :class:`~sage.misc.fast_methods.WithEqualityById`::

        sage: isinstance(M, sage.misc.fast_methods.WithEqualityById)
        True
    """
    from sage.rings.infinity import infinity
    from sage.manifolds.differentiable.manifold import DifferentiableManifold
    from sage.manifolds.differentiable.pseudo_riemannian import PseudoRiemannianManifold
    from sage.manifolds.differentiable.degenerate import DegenerateManifold
    from sage.manifolds.topological_submanifold import TopologicalSubmanifold
    from sage.manifolds.differentiable.differentiable_submanifold import DifferentiableSubmanifold
    from sage.manifolds.differentiable.pseudo_riemannian_submanifold import PseudoRiemannianSubmanifold
    from sage.manifolds.differentiable.degenerate_submanifold import DegenerateSubmanifold

    global _manifold_id

    _manifold_id += 1
    unique_tag = lambda: getrandbits(128)*_manifold_id

    if structure is None:
        if any(extra_kwds.get(x, None) is not None
               for x in ('metric_name', 'metric_latex_name', 'signature')):
            structure = 'pseudo-Riemannian'

    if structure is None:
        diff_degree = extra_kwds.get('diff_degree', infinity)
        if diff_degree == infinity:
            structure = 'smooth'
        elif diff_degree > 0:
            structure = 'differentiable'
        else:
            structure = 'topological'

    if structure in ['topological', 'top']:
        if field == 'real' or isinstance(field, sage.rings.abc.RealField):
            structure = RealTopologicalStructure()
        else:
            structure = TopologicalStructure()
        if 'ambient' in extra_kwds:
            ambient = extra_kwds['ambient']
            return TopologicalSubmanifold(dim, name, field, structure,
                                          ambient=ambient,
                                          latex_name=latex_name,
                                          start_index=start_index,
                                          unique_tag=unique_tag())
        return TopologicalManifold(dim, name, field, structure,
                                   latex_name=latex_name,
                                   start_index=start_index,
                                   unique_tag=unique_tag())
    elif structure in ['differentiable', 'diff', 'smooth']:
        if 'diff_degree' in extra_kwds:
            diff_degree = extra_kwds['diff_degree']
            if structure == 'smooth' and diff_degree != infinity:
                raise ValueError("diff_degree = {} is ".format(diff_degree) +
                                 "not compatible with a smooth structure")
        else:
            diff_degree = infinity
        if field == 'real' or isinstance(field, sage.rings.abc.RealField):
            structure = RealDifferentialStructure()
        else:
            structure = DifferentialStructure()
        if 'ambient' in extra_kwds:
            ambient = extra_kwds['ambient']
            return DifferentiableSubmanifold(dim, name, field, structure,
                                             ambient=ambient,
                                             diff_degree=diff_degree,
                                             latex_name=latex_name,
                                             start_index=start_index,
                                             unique_tag=unique_tag())
        return DifferentiableManifold(dim, name, field, structure,
                                      diff_degree=diff_degree,
                                      latex_name=latex_name,
                                      start_index=start_index,
                                      unique_tag=unique_tag())
    elif structure in ['pseudo-Riemannian', 'Riemannian', 'Lorentzian','degenerate_metric']:
        diff_degree = extra_kwds.get('diff_degree', infinity)
        metric_name = extra_kwds.get('metric_name', None)
        metric_latex_name = extra_kwds.get('metric_latex_name', None)
        if structure == 'pseudo-Riemannian':
            signature = extra_kwds.get('signature', None)
        elif structure == 'Riemannian':
            signature = dim
        elif structure == 'degenerate_metric':
            signature = (0, dim-1, 1)
        elif structure == 'Lorentzian':
            if 'signature' in extra_kwds:
                signat = extra_kwds['signature']
                if signat == 'positive' or signat == dim - 2:
                    signature = dim - 2
                elif signat == 'negative' or signat == 2 - dim:
                    signature = 2 - dim
                else:
                    raise ValueError("signature {} not ".format(signat) +
                                     "compatible with a Lorentzian " +
                                     "manifold of dimension {}".format(dim))
            else:
                signature = dim - 2  # default value for a Lorentzian manifold
        if 'ambient' in extra_kwds:
            ambient = extra_kwds['ambient']
            if structure == 'degenerate_metric':
                return DegenerateSubmanifold(dim, name, ambient=ambient,
                                               metric_name=metric_name,
                                               signature=signature,
                                               diff_degree=diff_degree,
                                               latex_name=latex_name,
                                               metric_latex_name=metric_latex_name,
                                               start_index=start_index,
                                               unique_tag=unique_tag())
            return PseudoRiemannianSubmanifold(dim, name, ambient=ambient,
                                               metric_name=metric_name,
                                               signature=signature,
                                               diff_degree=diff_degree,
                                               latex_name=latex_name,
                                               metric_latex_name=metric_latex_name,
                                               start_index=start_index,
                                               unique_tag=unique_tag())
        if structure == 'degenerate_metric':
                return DegenerateManifold(dim, name, metric_name=metric_name,
                                               signature=signature,
                                               diff_degree=diff_degree,
                                               latex_name=latex_name,
                                               metric_latex_name=metric_latex_name,
                                               start_index=start_index,
                                               unique_tag=unique_tag())
        return PseudoRiemannianManifold(dim, name, metric_name=metric_name,
                                        signature=signature,
                                        diff_degree=diff_degree,
                                        latex_name=latex_name,
                                        metric_latex_name=metric_latex_name,
                                        start_index=start_index,
                                        unique_tag=unique_tag())
    raise NotImplementedError("manifolds of type {} are ".format(structure) +
                              "not implemented")


Manifold.options = TopologicalManifold.options
