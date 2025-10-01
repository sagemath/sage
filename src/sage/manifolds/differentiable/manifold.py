r"""
Differentiable Manifolds

Given a non-discrete topological field `K` (in most applications, `K = \RR` or
`K = \CC`; see however [Ser1992]_ for `K = \QQ_p` and [Ber2008]_ for other fields),
a *differentiable manifold over* `K` is a topological manifold `M` over `K`
equipped with an atlas whose transitions maps are of class `C^k` (i.e.
`k`-times  continuously differentiable) for a fixed positive integer `k`
(possibly `k=\infty`). `M` is then called a `C^k`-*manifold over* `K`.

Note that

- if the mention of `K` is omitted, then `K=\RR` is assumed;
- if `K=\CC`, any `C^k`-manifold with `k\geq 1` is actually a
  `C^\infty`-manifold (even an analytic manifold);
- if `K=\RR`, any `C^k`-manifold with `k\geq 1` admits a compatible
  `C^\infty`-structure (Whitney's smoothing theorem).

Differentiable manifolds are implemented via the class
:class:`DifferentiableManifold`.
Open subsets of differentiable manifolds are also implemented via
:class:`DifferentiableManifold`, since they are differentiable manifolds by
themselves.

The user interface is provided by the generic function
:func:`~sage.manifolds.manifold.Manifold`, with
the argument ``structure`` set to ``'differentiable'`` and the argument
``diff_degree`` set to `k`,  or the argument ``structure`` set to ``'smooth'``
(the default value).

.. RUBRIC:: Example 1: the 2-sphere as a differentiable manifold of dimension
  2 over `\RR`

One starts by declaring `S^2` as a 2-dimensional differentiable manifold::

    sage: M = Manifold(2, 'S^2')
    sage: M
    2-dimensional differentiable manifold S^2

Since the base topological field has not been specified in the argument list
of ``Manifold``, `\RR` is assumed::

    sage: M.base_field()
    Real Field with 53 bits of precision
    sage: dim(M)
    2

By default, the created object is a smooth manifold::

   sage: M.diff_degree()
    +Infinity

Let us consider the complement of a point, the "North pole" say; this is an
open subset of `S^2`, which we call `U`::

    sage: U = M.open_subset('U'); U
    Open subset U of the 2-dimensional differentiable manifold S^2

A standard chart on `U` is provided by the stereographic projection from the
North pole to the equatorial plane::

    sage: stereoN.<x,y> = U.chart(); stereoN
    Chart (U, (x, y))

Thanks to the operator ``<x,y>`` on the left-hand side, the coordinates
declared in a chart (here `x` and `y`), are accessible by their names; they are
Sage's symbolic variables::

    sage: y
    y
    sage: type(y)
    <class 'sage.symbolic.expression.Expression'>

The South pole is the point of coordinates `(x,y)=(0,0)` in the above
chart::

    sage: S = U.point((0,0), chart=stereoN, name='S'); S
    Point S on the 2-dimensional differentiable manifold S^2

Let us call `V` the open subset that is the complement of the South pole and
let us introduce on it the chart induced by the stereographic projection from
the South pole to the equatorial plane::

    sage: V = M.open_subset('V'); V
    Open subset V of the 2-dimensional differentiable manifold S^2
    sage: stereoS.<u,v> = V.chart(); stereoS
    Chart (V, (u, v))

The North pole is the point of coordinates `(u,v)=(0,0)` in this chart::

    sage: N = V.point((0,0), chart=stereoS, name='N'); N
    Point N on the 2-dimensional differentiable manifold S^2

To fully construct the manifold, we declare that it is the union of `U`
and `V`::

    sage: M.declare_union(U,V)

and we provide the transition map between the charts ``stereoN`` = `(U, (x, y))`
and ``stereoS`` = `(V, (u, v))`, denoting by `W` the intersection of `U` and
`V` (`W` is the subset of `U` defined by `x^2+y^2\not=0`, as well as the subset
of `V` defined by `u^2+v^2\not=0`)::

    sage: stereoN_to_S = stereoN.transition_map(stereoS,
    ....:                [x/(x^2+y^2), y/(x^2+y^2)], intersection_name='W',
    ....:                restrictions1= x^2+y^2!=0, restrictions2= u^2+v^2!=0)
    sage: stereoN_to_S
    Change of coordinates from Chart (W, (x, y)) to Chart (W, (u, v))
    sage: stereoN_to_S.display()
    u = x/(x^2 + y^2)
    v = y/(x^2 + y^2)

We give the name ``W`` to the Python variable representing `W=U\cap V`::

    sage: W = U.intersection(V)

The inverse of the transition map is computed by the method ``inverse()``::

    sage: stereoN_to_S.inverse()
    Change of coordinates from Chart (W, (u, v)) to Chart (W, (x, y))
    sage: stereoN_to_S.inverse().display()
    x = u/(u^2 + v^2)
    y = v/(u^2 + v^2)

At this stage, we have four open subsets on `S^2`::

    sage: M.subset_family()
    Set {S^2, U, V, W} of open subsets of the 2-dimensional differentiable manifold S^2

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
    [Chart (U, (x, y)), Chart (V, (u, v)), Chart (W, (x, y)), Chart (W, (u, v))]

Let us consider the point of coordinates (1,2) in the chart ``stereoN``::

    sage: p = M.point((1,2), chart=stereoN, name='p'); p
    Point p on the 2-dimensional differentiable manifold S^2
    sage: p.parent()
    2-dimensional differentiable manifold S^2
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

A differentiable scalar field on the sphere::

    sage: f = M.scalar_field({stereoN: atan(x^2+y^2), stereoS: pi/2-atan(u^2+v^2)},
    ....:                    name='f')
    sage: f
    Scalar field f on the 2-dimensional differentiable manifold S^2
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
    Algebra of differentiable scalar fields on the 2-dimensional differentiable
     manifold S^2
    sage: f.parent().category()
    Join of Category of commutative algebras over Symbolic Ring and Category of homsets of topological spaces

A differentiable manifold has a default vector frame, which, unless otherwise
specified, is the coordinate frame associated with the first defined chart::

    sage: M.default_frame()
    Coordinate frame (U, (∂/∂x,∂/∂y))
    sage: latex(M.default_frame())
    \left(U, \left(\frac{\partial}{\partial x },\frac{\partial}{\partial y }\right)\right)
    sage: M.default_frame() is stereoN.frame()
    True

A vector field on the sphere::

    sage: w = M.vector_field(name='w')
    sage: w[stereoN.frame(), :] = [x, y]
    sage: w.add_comp_by_continuation(stereoS.frame(), W, stereoS)
    sage: w.display() # display in the default frame (stereoN.frame())
    w = x ∂/∂x + y ∂/∂y
    sage: w.display(stereoS.frame())
    w = -u ∂/∂u - v ∂/∂v
    sage: w.parent()
    Module X(S^2) of vector fields on the 2-dimensional differentiable
     manifold S^2
    sage: w.parent().category()
    Category of modules over Algebra of differentiable scalar fields on the
     2-dimensional differentiable manifold S^2

Vector fields act on scalar fields::

    sage: w(f)
    Scalar field w(f) on the 2-dimensional differentiable manifold S^2
    sage: w(f).display()
    w(f): S^2 → ℝ
    on U: (x, y) ↦ 2*(x^2 + y^2)/(x^4 + 2*x^2*y^2 + y^4 + 1)
    on V: (u, v) ↦ 2*(u^2 + v^2)/(u^4 + 2*u^2*v^2 + v^4 + 1)
    sage: w(f) == f.differential()(w)
    True

The value of the vector field at point `p` is a vector tangent to the sphere::

    sage: w.at(p)
    Tangent vector w at Point p on the 2-dimensional differentiable manifold S^2
    sage: w.at(p).display()
    w = ∂/∂x + 2 ∂/∂y
    sage: w.at(p).parent()
    Tangent space at Point p on the 2-dimensional differentiable manifold S^2

A 1-form on the sphere::

    sage: df = f.differential() ; df
    1-form df on the 2-dimensional differentiable manifold S^2
    sage: df.display()
    df = 2*x/(x^4 + 2*x^2*y^2 + y^4 + 1) dx + 2*y/(x^4 + 2*x^2*y^2 + y^4 + 1) dy
    sage: df.display(stereoS.frame())
    df = -2*u/(u^4 + 2*u^2*v^2 + v^4 + 1) du - 2*v/(u^4 + 2*u^2*v^2 + v^4 + 1) dv
    sage: df.parent()
    Module Omega^1(S^2) of 1-forms on the 2-dimensional differentiable
     manifold S^2
    sage: df.parent().category()
    Category of modules over Algebra of differentiable scalar fields on the
     2-dimensional differentiable manifold S^2

The value of the 1-form at point `p` is a linear form on the tangent space
at `p`::

    sage: df.at(p)
    Linear form df on the Tangent space at Point p on the 2-dimensional
     differentiable manifold S^2
    sage: df.at(p).display()
    df = 1/13 dx + 2/13 dy
    sage: df.at(p).parent()
    Dual of the Tangent space at Point p on the 2-dimensional differentiable
     manifold S^2


.. RUBRIC:: Example 2: the Riemann sphere as a differentiable manifold of
  dimension 1 over `\CC`

We declare the Riemann sphere `\CC^*` as a 1-dimensional differentiable
manifold over `\CC`::

    sage: M = Manifold(1, 'ℂ*', field='complex'); M
    1-dimensional complex manifold ℂ*

We introduce a first open subset, which is actually
`\CC = \CC^*\setminus\{\infty\}` if we interpret `\CC^*` as the Alexandroff
one-point compactification of `\CC`::

    sage: U = M.open_subset('U')

A natural chart on `U` is then nothing but the identity map of `\CC`, hence
we denote the associated coordinate by `z`::

    sage: Z.<z> = U.chart()

The origin of the complex plane is the point of coordinate `z=0`::

    sage: O = U.point((0,), chart=Z, name='O'); O
    Point O on the 1-dimensional complex manifold ℂ*

Another open subset of `\CC^*` is `V = \CC^*\setminus\{O\}`::

    sage: V = M.open_subset('V')

We define a chart on `V` such that the point at infinity is the point of
coordinate 0 in this chart::

    sage: W.<w> = V.chart(); W
    Chart (V, (w,))
    sage: inf = M.point((0,), chart=W, name='inf', latex_name=r'\infty')
    sage: inf
    Point inf on the 1-dimensional complex manifold ℂ*

To fully construct the Riemann sphere, we declare that it is the union of `U`
and `V`::

    sage: M.declare_union(U,V)

and we provide the transition map between the two charts as `w=1/z` on
on `A = U\cap V`::

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
    Point i on the 1-dimensional complex manifold ℂ*

Its coordinates with respect to the charts ``Z`` and ``W`` are::

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
    Set {A, U, V, ℂ*} of open subsets of the 1-dimensional complex manifold ℂ*
    sage: M.atlas()
    [Chart (U, (z,)), Chart (V, (w,)), Chart (A, (z,)), Chart (A, (w,))]

A constant map `\CC^* \rightarrow \CC`::

    sage: f = M.constant_scalar_field(3+2*I, name='f'); f
    Scalar field f on the 1-dimensional complex manifold ℂ*
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
    Algebra of differentiable scalar fields on the 1-dimensional complex
     manifold ℂ*
    sage: f.parent().category()
    Join of Category of commutative algebras over Symbolic Ring and Category of homsets of topological spaces

A vector field on the Riemann sphere::

    sage: v = M.vector_field(name='v')
    sage: v[Z.frame(), 0] = z^2
    sage: v.add_comp_by_continuation(W.frame(), U.intersection(V), W)
    sage: v.display(Z.frame())
    v = z^2 ∂/∂z
    sage: v.display(W.frame())
    v = -∂/∂w
    sage: v.parent()
    Module X(ℂ*) of vector fields on the 1-dimensional complex manifold ℂ*

The vector field `v` acting on the scalar field `f`::

    sage: v(f)
    Scalar field zero on the 1-dimensional complex manifold ℂ*

Since `f` is constant, `v(f)` is vanishing::

    sage: v(f).display()
    zero: ℂ* → ℂ
    on U: z ↦ 0
    on V: w ↦ 0

The value of the vector field `v` at the point `\infty` is a vector tangent to
the Riemann sphere::

    sage: v.at(inf)
    Tangent vector v at Point inf on the 1-dimensional complex manifold ℂ*
    sage: v.at(inf).display()
    v = -∂/∂w
    sage: v.at(inf).parent()
    Tangent space at Point inf on the 1-dimensional complex manifold ℂ*

AUTHORS:

- Eric Gourgoulhon (2015): initial version
- Travis Scrimshaw (2016): review tweaks
- Michael Jung (2020): tensor bundles and orientability
- Matthias Koeppe (2021): refactoring of subsets code

REFERENCES:

- [Lee2013]_
- [KN1963]_
- [Huy2005]_
- [Ser1992]_
- [Ber2008]_
- [BG1988]_
"""

# ****************************************************************************
#       Copyright (C) 2015-2019 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015      Michal Bejger <bejger@camk.edu.pl>
#       Copyright (C) 2015-2016 Travis Scrimshaw <tscrimsh@umn.edu>
#       Copyright (C) 2017      Karim Van Aelst
#       Copyright (C) 2019      Hans Fotsing Tetsing
#       Copyright (C) 2019-2020 Michael Jung
#       Copyright (C) 2021      Matthias Koeppe
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from __future__ import annotations

from typing import TYPE_CHECKING, Optional, Union

from sage.categories.homset import Hom
from sage.categories.manifolds import Manifolds
from sage.manifolds.differentiable.mixed_form_algebra import MixedFormAlgebra
from sage.manifolds.manifold import TopologicalManifold
from sage.rings.cc import CC
from sage.rings.infinity import infinity, minus_infinity
from sage.rings.integer import Integer
from sage.rings.real_mpfr import RR

if TYPE_CHECKING:
    from sage.manifolds.differentiable.diff_form import DiffForm
    from sage.manifolds.differentiable.diff_map import DiffMap
    from sage.manifolds.differentiable.metric import PseudoRiemannianMetric
    from sage.manifolds.differentiable.vectorfield_module import (
        VectorFieldFreeModule,
        VectorFieldModule,
    )
    from sage.manifolds.differentiable.vectorframe import VectorFrame

###############################################################################


class DifferentiableManifold(TopologicalManifold):
    r"""
    Differentiable manifold over a topological field `K`.

    Given a non-discrete topological field `K` (in most applications,
    `K = \RR` or `K = \CC`; see however [Ser1992]_ for `K = \QQ_p` and
    [Ber2008]_ for other fields), a *differentiable manifold over* `K` is a
    topological manifold `M` over `K` equipped with an atlas whose transitions
    maps are of class `C^k` (i.e. `k`-times  continuously differentiable) for
    a fixed positive integer `k` (possibly `k=\infty`). `M` is then called a
    `C^k`-*manifold over* `K`.

    Note that

    - if the mention of `K` is omitted, then `K=\RR` is assumed;
    - if `K=\CC`, any `C^k`-manifold with `k\geq 1` is actually a
      `C^\infty`-manifold (even an analytic manifold);
    - if `K=\RR`, any `C^k`-manifold with `k\geq 1` admits a compatible
      `C^\infty`-structure (Whitney's smoothing theorem).

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
      :class:`~sage.manifolds.structure.DifferentialStructure` or
      :class:`~sage.manifolds.structure.RealDifferentialStructure`)
    - ``base_manifold`` -- (default: ``None``) if not ``None``, must be a
      differentiable manifold; the created object is then an open subset of
      ``base_manifold``
    - ``diff_degree`` -- (default: ``infinity``) degree `k` of
      differentiability
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the manifold; if none is provided, it is set to ``name``
    - ``start_index`` -- (default: 0) integer; lower value of the range of
      indices used for "indexed objects" on the manifold, e.g. coordinates
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
      :class:`~sage.manifolds.subset.ManifoldSubset`,
      via :class:`~sage.manifolds.manifold.TopologicalManifold`,
      would return the previously constructed object corresponding to these
      arguments).

    EXAMPLES:

    A 4-dimensional differentiable manifold (over `\RR`)::

        sage: M = Manifold(4, 'M', latex_name=r'\mathcal{M}'); M
        4-dimensional differentiable manifold M
        sage: type(M)
        <class 'sage.manifolds.differentiable.manifold.DifferentiableManifold_with_category'>
        sage: latex(M)
        \mathcal{M}
        sage: dim(M)
        4

    Since the base field has not been specified, `\RR` has been assumed::

        sage: M.base_field()
        Real Field with 53 bits of precision

    Since the degree of differentiability has not been specified, the default
    value, `C^\infty`, has been assumed::

        sage: M.diff_degree()
        +Infinity

    The input parameter ``start_index`` defines the range of indices on the
    manifold::

        sage: M = Manifold(4, 'M')
        sage: list(M.irange())
        [0, 1, 2, 3]
        sage: M = Manifold(4, 'M', start_index=1)
        sage: list(M.irange())
        [1, 2, 3, 4]
        sage: list(Manifold(4, 'M', start_index=-2).irange())
        [-2, -1, 0, 1]

    A complex manifold::

        sage: N = Manifold(3, 'N', field='complex'); N
        3-dimensional complex manifold N

    A differentiable manifold over `\QQ_5`, the field of 5-adic numbers::

        sage: N = Manifold(2, 'N', field=Qp(5)); N
        2-dimensional differentiable manifold N over the 5-adic Field with
         capped relative precision 20

    A differentiable manifold is of course a topological manifold::

        sage: isinstance(M, sage.manifolds.manifold.TopologicalManifold)
        True
        sage: isinstance(N, sage.manifolds.manifold.TopologicalManifold)
        True

    A differentiable manifold is a Sage *parent* object, in the category of
    differentiable (here smooth) manifolds over a given topological field (see
    :class:`~sage.categories.manifolds.Manifolds`)::

        sage: isinstance(M, Parent)
        True
        sage: M.category()
        Category of smooth manifolds over Real Field with 53 bits of precision
        sage: from sage.categories.manifolds import Manifolds
        sage: M.category() is Manifolds(RR).Smooth()
        True
        sage: M.category() is Manifolds(M.base_field()).Smooth()
        True
        sage: M in Manifolds(RR).Smooth()
        True
        sage: N in Manifolds(Qp(5)).Smooth()
        True

    The corresponding Sage *elements* are points::

        sage: X.<t, x, y, z> = M.chart()
        sage: p = M.an_element(); p
        Point on the 4-dimensional differentiable manifold M
        sage: p.parent()
        4-dimensional differentiable manifold M
        sage: M.is_parent_of(p)
        True
        sage: p in M
        True

    The manifold's points are instances of class
    :class:`~sage.manifolds.point.ManifoldPoint`::

        sage: isinstance(p, sage.manifolds.point.ManifoldPoint)
        True

    Since an open subset of a differentiable manifold `M` is itself a
    differentiable manifold, open subsets of `M` have all attributes of
    manifolds::

        sage: U = M.open_subset('U', coord_def={X: t>0}); U
        Open subset U of the 4-dimensional differentiable manifold M
        sage: U.category()
        Join of Category of subobjects of sets and Category of smooth manifolds
         over Real Field with 53 bits of precision
        sage: U.base_field() == M.base_field()
        True
        sage: dim(U) == dim(M)
        True

    The manifold passes all the tests of the test suite relative to its
    category::

        sage: TestSuite(M).run()
    """
    def __init__(self, n, name, field, structure, base_manifold=None,
                 diff_degree=infinity, latex_name=None, start_index=0,
                 category=None, unique_tag=None):
        r"""
        Construct a differentiable manifold.

        TESTS::

            sage: M = Manifold(3, 'M', latex_name=r'\mathbb{M}',
            ....:                  start_index=1)
            sage: M
            3-dimensional differentiable manifold M
            sage: latex(M)
            \mathbb{M}
            sage: dim(M)
            3
            sage: X.<x,y,z> = M.chart()
            sage: TestSuite(M).run()

        Tests for open subsets, which are constructed as differentiable
        manifolds::

            sage: U = M.open_subset('U', coord_def={X: x>0})
            sage: type(U)
            <class 'sage.manifolds.differentiable.manifold.DifferentiableManifold_with_category'>
            sage: U.category() is M.category().Subobjects()
            True
            sage: TestSuite(U).run()
        """
        if base_manifold is None:
            if category is None:
                if field == 'real':
                    field_c = RR
                elif field == 'complex':
                    field_c = CC
                else:
                    field_c = field
                if diff_degree == infinity:
                    category = Manifolds(field_c).Smooth()
                else:
                    category = Manifolds(field_c).Differentiable()
        elif not isinstance(base_manifold, DifferentiableManifold):
            raise TypeError("the argument 'base_manifold' must be a " +
                            "differentiable manifold")
        TopologicalManifold.__init__(self, n, name, field, structure,
                                     base_manifold=base_manifold,
                                     latex_name=latex_name,
                                     start_index=start_index,
                                     category=category)
        # The degree of differentiability:
        if diff_degree == infinity:
            self._diff_degree = infinity
        elif not isinstance(diff_degree, (int, Integer)):
            raise TypeError("the argument 'diff_degree' must be an integer")
        elif diff_degree < 1:
            raise ValueError("the argument 'diff_degree' must be a positive " +
                             "integer")
        else:
            self._diff_degree = diff_degree
        # Vector frames:
        self._frames = []  # list of vector frames defined on subsets of self
        # list of vector frames defined on subsets of self that are
        # not subframes of frames on larger subsets
        self._top_frames = []
        self._def_frame = None  # default frame
        self._coframes = []  # list of coframes defined on subsets of self
        # List of vector frames that individually cover self, i.e. whose
        # domains are self (if non-empty, self is parallelizable):
        self._covering_frames = []
        self._parallelizable_parts = set() # parallelizable subsets contained in self
        self._frame_changes = {} # dictionary of changes of frames
        # Dictionary of vector field modules along self
        # (keys = diff. map from self to an open set (possibly the identity map))
        self._vector_field_modules = {} # dict of all established vector field
                                        # modules
        self._tensor_bundles = {} # dict of dict of all established tensor
                                  # bundles

    def diff_degree(self):
        r"""
        Return the manifold's degree of differentiability.

        The degree of differentiability is the integer `k` (possibly
        `k=\infty`) such that the manifold is a `C^k`-manifold over its base
        field.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: M.diff_degree()
            +Infinity
            sage: M = Manifold(2, 'M', structure='differentiable', diff_degree=3)
            sage: M.diff_degree()
            3
        """
        return self._diff_degree

    def open_subset(self, name, latex_name=None, coord_def={}, supersets=None):
        r"""
        Create an open subset of the manifold.

        An open subset is a set that is (i) included in the manifold and (ii)
        open with respect to the manifold's topology. It is a differentiable
        manifold by itself. Hence the returned object is an instance of
        :class:`DifferentiableManifold`.

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

        OUTPUT: the open subset, as an instance of :class:`DifferentiableManifold`

        EXAMPLES:

        Creating an open subset of a differentiable manifold::

            sage: M = Manifold(2, 'M')
            sage: A = M.open_subset('A'); A
            Open subset A of the 2-dimensional differentiable manifold M

        As an open subset of a differentiable manifold, ``A`` is itself a
        differentiable manifold, on the same topological field and of the same
        dimension as ``M``::

            sage: A.category()
            Join of Category of subobjects of sets and Category of smooth
             manifolds over Real Field with 53 bits of precision
            sage: A.base_field() == M.base_field()
            True
            sage: dim(A) == dim(M)
            True

        Creating an open subset of ``A``::

            sage: B = A.open_subset('B'); B
            Open subset B of the 2-dimensional differentiable manifold M

        We have then::

            sage: A.subset_family()
            Set {A, B} of open subsets of the 2-dimensional differentiable manifold M
            sage: B.is_subset(A)
            True
            sage: B.is_subset(M)
            True

        Defining an open subset by some coordinate restrictions: the open
        unit disk in of the Euclidean plane::

            sage: X.<x,y> = M.chart() # Cartesian coordinates on M
            sage: U = M.open_subset('U', coord_def={X: x^2+y^2<1}); U
            Open subset U of the 2-dimensional differentiable manifold M

        Since the argument ``coord_def`` has been set, ``U`` is automatically
        endowed with a chart, which is the restriction of ``X``
        to ``U``::

            sage: U.atlas()
            [Chart (U, (x, y))]
            sage: U.default_chart()
            Chart (U, (x, y))
            sage: U.default_chart() is X.restrict(U)
            True

        A point in ``U``::

            sage: p = U.an_element(); p
            Point on the 2-dimensional differentiable manifold M
            sage: X(p)  # the coordinates (x,y) of p
            (0, 0)
            sage: p in U
            True

        Checking whether various points, defined by their coordinates
        with respect to chart ``X``,  are in ``U``::

            sage: M((0,1/2)) in U
            True
            sage: M((0,1)) in U
            False
            sage: M((1/2,1)) in U
            False
            sage: M((-1/2,1/3)) in U
            True
        """
        resu = DifferentiableManifold(self._dim, name, self._field,
                                      self._structure, base_manifold=self._manifold,
                                      diff_degree=self._diff_degree,
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

            sage: M = Manifold(2, 'R^2', structure='differentiable')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: from sage.manifolds.differentiable.manifold import DifferentiableManifold
            sage: U = DifferentiableManifold(2, 'U', field=M._field, structure=M._structure, base_manifold=M)
            sage: M._init_open_subset(U, coord_def={c_cart: x^2+y^2<1})
            sage: U
            Open subset U of the 2-dimensional differentiable manifold R^2
        """
        super()._init_open_subset(resu, coord_def=coord_def)
        #!# update vector frames and change of frames

    def diff_map(self, codomain, coord_functions=None, chart1=None,
                       chart2=None, name=None, latex_name=None):
        r"""
        Define a differentiable map between the current differentiable manifold
        and a differentiable manifold over the same topological field.

        See :class:`~sage.manifolds.differentiable.diff_map.DiffMap` for a
        complete documentation.

        INPUT:

        - ``codomain`` -- the map codomain (a differentiable manifold over the
          same topological field as the current differentiable manifold)
        - ``coord_functions`` -- (default: ``None``) if not ``None``, must be
          either

          - (i) a dictionary of
            the coordinate expressions (as lists (or tuples) of the
            coordinates of the image expressed in terms of the coordinates of
            the considered point) with the pairs of charts (chart1, chart2)
            as keys (chart1 being a chart on the current manifold and chart2 a
            chart on ``codomain``)
          - (ii) a single coordinate expression in a given pair of charts, the
            latter being provided by the arguments ``chart1`` and ``chart2``

          In both cases, if the dimension of the arrival manifold is 1,
          a single coordinate expression can be passed instead of a tuple with
          a single element
        - ``chart1`` -- (default: ``None``; used only in case (ii) above) chart
          on the current manifold defining the start coordinates involved in
          ``coord_functions`` for case (ii); if none is provided, the
          coordinates are assumed to refer to the manifold's default chart
        - ``chart2`` -- (default: ``None``; used only in case (ii) above) chart
          on ``codomain`` defining the arrival coordinates involved in
          ``coord_functions`` for case (ii); if none is provided, the
          coordinates are assumed to refer to the default chart of ``codomain``
        - ``name`` -- (default: ``None``) name given to the differentiable
          map
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          differentiable map; if none is provided, the LaTeX symbol is set to
          ``name``

        OUTPUT:

        - the differentiable map, as an instance of
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`

        EXAMPLES:

        A differentiable map between an open subset of `S^2` covered by regular
        spherical coordinates and `\RR^3`::

            sage: M = Manifold(2, 'S^2')
            sage: U = M.open_subset('U')
            sage: c_spher.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
            sage: N = Manifold(3, 'R^3', r'\RR^3')
            sage: c_cart.<x,y,z> = N.chart()  # Cartesian coord. on R^3
            sage: Phi = U.diff_map(N, (sin(th)*cos(ph), sin(th)*sin(ph), cos(th)),
            ....:                  name='Phi', latex_name=r'\Phi')
            sage: Phi
            Differentiable map Phi from the Open subset U of the 2-dimensional
             differentiable manifold S^2 to the 3-dimensional differentiable
             manifold R^3

        The same definition, but with a dictionary with pairs of charts as
        keys (case (i) above)::

            sage: Phi1 = U.diff_map(N,
            ....:        {(c_spher, c_cart): (sin(th)*cos(ph), sin(th)*sin(ph),
            ....:         cos(th))}, name='Phi', latex_name=r'\Phi')
            sage: Phi1 == Phi
            True

        The differentiable map acting on a point::

            sage: p = U.point((pi/2, pi)) ; p
            Point on the 2-dimensional differentiable manifold S^2
            sage: Phi(p)
            Point on the 3-dimensional differentiable manifold R^3
            sage: Phi(p).coord(c_cart)
            (-1, 0, 0)
            sage: Phi1(p) == Phi(p)
            True

        See the documentation of class
        :class:`~sage.manifolds.differentiable.diff_map.DiffMap` for more
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
        return homset(coord_functions, name=name, latex_name=latex_name)

    def diffeomorphism(self, codomain=None, coord_functions=None, chart1=None,
                       chart2=None, name=None, latex_name=None):
        r"""
        Define a diffeomorphism between the current manifold and another one.

        See :class:`~sage.manifolds.differentiable.diff_map.DiffMap` for a
        complete documentation.

        INPUT:

        - ``codomain`` -- (default: ``None``) codomain of the diffeomorphism (the arrival manifold
          or some subset of it). If ``None``, the current manifold is taken.
        - ``coord_functions`` -- (default: ``None``) if not ``None``, must be
          either

          - (i) a dictionary of
            the coordinate expressions (as lists (or tuples) of the
            coordinates of the image expressed in terms of the coordinates of
            the considered point) with the pairs of charts (chart1, chart2)
            as keys (chart1 being a chart on the current manifold and chart2
            a chart on ``codomain``)
          - (ii) a single coordinate expression in a given pair of charts, the
            latter being provided by the arguments ``chart1`` and ``chart2``

          In both cases, if the dimension of the arrival manifold is 1,
          a single coordinate expression can be passed instead of a tuple with
          a single element
        - ``chart1`` -- (default: ``None``; used only in case (ii) above) chart
          on the current manifold defining the start coordinates involved in
          ``coord_functions`` for case (ii); if none is provided, the
          coordinates are assumed to refer to the manifold's default chart
        - ``chart2`` -- (default: ``None``; used only in case (ii) above) chart
          on ``codomain`` defining the arrival coordinates involved in
          ``coord_functions`` for case (ii); if none is provided, the
          coordinates are assumed to refer to the default chart of ``codomain``
        - ``name`` -- (default: ``None``) name given to the diffeomorphism
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          diffeomorphism; if none is provided, the LaTeX symbol is set to
          ``name``

        OUTPUT:

        - the diffeomorphism, as an instance of
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`

        EXAMPLES:

        Diffeomorphism between the open unit disk in `\RR^2` and `\RR^2`::

            sage: M = Manifold(2, 'M')  # the open unit disk
            sage: forget()  # for doctests only
            sage: c_xy.<x,y> = M.chart('x:(-1,1) y:(-1,1)', coord_restrictions=lambda x,y: x^2+y^2<1)
            ....:    # Cartesian coord on M
            sage: N = Manifold(2, 'N')  # R^2
            sage: c_XY.<X,Y> = N.chart()  # canonical coordinates on R^2
            sage: Phi = M.diffeomorphism(N, [x/sqrt(1-x^2-y^2), y/sqrt(1-x^2-y^2)],
            ....:                        name='Phi', latex_name=r'\Phi')
            sage: Phi
            Diffeomorphism Phi from the 2-dimensional differentiable manifold M
             to the 2-dimensional differentiable manifold N
            sage: Phi.display()
            Phi: M → N
               (x, y) ↦ (X, Y) = (x/sqrt(-x^2 - y^2 + 1), y/sqrt(-x^2 - y^2 + 1))

        The inverse diffeomorphism::

            sage: Phi^(-1)
            Diffeomorphism Phi^(-1) from the 2-dimensional differentiable
             manifold N to the 2-dimensional differentiable manifold M
            sage: (Phi^(-1)).display()
            Phi^(-1): N → M
               (X, Y) ↦ (x, y) = (X/sqrt(X^2 + Y^2 + 1), Y/sqrt(X^2 + Y^2 + 1))

        See the documentation of class
        :class:`~sage.manifolds.differentiable.diff_map.DiffMap` for more
        examples.
        """
        if codomain is None:
            codomain = self

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

    def vector_bundle(self, rank, name, field='real', latex_name=None):
        r"""
        Return a differentiable vector bundle over the given field with given
        rank over this differentiable manifold of the same differentiability
        class as the manifold.

        INPUT:

        - ``rank`` -- rank of the vector bundle
        - ``name`` -- name given to the total space
        - ``field`` -- (default: ``'real'``) topological field giving the
          vector space structure to the fibers
        - ``latex_name`` -- (optional) LaTeX name for the total space

        OUTPUT:

        - a differentiable vector bundle as an instance of
          :class:`~sage.manifolds.differentiable.vector_bundle.DifferentiableVectorBundle`

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: M.vector_bundle(2, 'E')
            Differentiable real vector bundle E -> M of rank 2 over the base
             space 2-dimensional differentiable manifold M
        """
        from sage.manifolds.differentiable.vector_bundle import (
            DifferentiableVectorBundle,
        )
        return DifferentiableVectorBundle(rank, name, self, field=field,
                                          latex_name=latex_name)

    def tangent_bundle(self, dest_map=None):
        r"""
        Return the tangent bundle possibly along a destination map with base
        space ``self``.

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.vector_bundle.TensorBundle`
            for complete documentation.

        INPUT:

        - ``dest_map`` -- (default: ``None``) destination map
          `\Phi:\ M \rightarrow N`
          (type: :class:`~sage.manifolds.differentiable.diff_map.DiffMap`) from
          which the tangent bundle is pulled back; if
          ``None``, it is assumed that `N=M` and `\Phi` is the identity map of
          `M` (case of the standard tangent bundle over `M`)

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: TM = M.tangent_bundle(); TM
            Tangent bundle TM over the 2-dimensional differentiable manifold M
        """
        return self.tensor_bundle(1, 0, dest_map=dest_map)

    def cotangent_bundle(self, dest_map=None):
        r"""
        Return the cotangent bundle possibly along a destination map with base
        space ``self``.

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.vector_bundle.TensorBundle`
            for complete documentation.

        INPUT:

        - ``dest_map`` -- (default: ``None``) destination map
          `\Phi:\ M \rightarrow N`
          (type: :class:`~sage.manifolds.differentiable.diff_map.DiffMap`) from
          which the cotangent bundle is pulled back; if
          ``None``, it is assumed that `N=M` and `\Phi` is the identity map of
          `M` (case of the standard tangent bundle over `M`)

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: cTM = M.cotangent_bundle(); cTM
            Cotangent bundle T*M over the 2-dimensional differentiable
             manifold M
        """
        return self.tensor_bundle(0, 1, dest_map=dest_map)

    def tensor_bundle(self, k, l, dest_map=None):
        r"""
        Return a tensor bundle of type `(k, l)` defined over ``self``, possibly
        along a destination map.

        INPUT:

        - ``k`` -- the contravariant rank of the tensor bundle
        - ``l`` -- the covariant rank of the tensor bundle
        - ``dest_map`` -- (default: ``None``) destination map
          `\Phi:\ M \rightarrow N`
          (type: :class:`~sage.manifolds.differentiable.diff_map.DiffMap`) from
          which the tensor bundle is pulled back; if
          ``None``, it is assumed that `N=M` and `\Phi` is the identity map of
          `M` (case of the standard tangent bundle over `M`)

        OUTPUT:

        - a
          :class:`~sage.manifolds.differentiable.vector_bundle.TensorBundle`
          representing a tensor bundle of type-`(k,l)` over ``self``

        EXAMPLES:

        A tensor bundle over a parallelizable 2-dimensional differentiable
        manifold::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: M.tensor_bundle(1, 2)
            Tensor bundle T^(1,2)M over the 2-dimensional differentiable
             manifold M

        The special case of the tangent bundle as tensor bundle of type (1,0)::

            sage: M.tensor_bundle(1,0)
            Tangent bundle TM over the 2-dimensional differentiable manifold M

        The result is cached::

            sage: M.tensor_bundle(1, 2) is M.tensor_bundle(1, 2)
            True

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.vector_bundle.TensorBundle`
            for more examples and documentation.
        """
        if dest_map is None:
            dest_map = self.identity_map()
        if dest_map not in self._tensor_bundles:
            from sage.manifolds.differentiable.vector_bundle import TensorBundle
            self._tensor_bundles[dest_map] = {(k, l):
                                              TensorBundle(self, k, l,
                                                           dest_map=dest_map)}
        else:
            if (k, l) not in self._tensor_bundles[dest_map]:
                from sage.manifolds.differentiable.vector_bundle import TensorBundle
                self._tensor_bundles[dest_map][(k, l)] = TensorBundle(self, k,
                                                           l, dest_map=dest_map)
        return self._tensor_bundles[dest_map][(k, l)]

    def vector_field_module(
        self, dest_map: Optional[DiffMap] = None, force_free: bool = False
    ) -> Union[VectorFieldModule, VectorFieldFreeModule]:
        r"""
        Return the set of vector fields defined on ``self``, possibly
        with values in another differentiable manifold, as a module over the
        algebra of scalar fields defined on the manifold.

        See :class:`~sage.manifolds.differentiable.vectorfield_module.VectorFieldModule`
        for a complete documentation.

        INPUT:

        - ``dest_map`` -- (default: ``None``) destination map, i.e. a
          differentiable map `\Phi:\ M \rightarrow N`, where `M` is the
          current manifold and `N` a differentiable manifold;
          if ``None``, it is assumed that `N = M` and that `\Phi` is the
          identity map (case of vector fields *on* `M`), otherwise
          ``dest_map`` must be a
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`
        - ``force_free`` -- boolean (default: ``False``); if set to ``True``, force
          the construction of a *free* module (this implies that `N` is
          parallelizable)

        OUTPUT:

        - a
          :class:`~sage.manifolds.differentiable.vectorfield_module.VectorFieldModule`
          (or if `N` is parallelizable, a
          :class:`~sage.manifolds.differentiable.vectorfield_module.VectorFieldFreeModule`)
          representing the `C^k(M)`-module `\mathfrak{X}(M,\Phi)` of vector
          fields on `M` taking values on `\Phi(M)\subset N`

        EXAMPLES:

        Vector field module `\mathfrak{X}(U) := \mathfrak{X}(U,\mathrm{Id}_U)`
        of the complement `U` of the two poles on the sphere `\mathbb{S}^2`::

            sage: S2 = Manifold(2, 'S^2')
            sage: U = S2.open_subset('U')  # the complement of the two poles
            sage: spher_coord.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi') # spherical coordinates
            sage: XU = U.vector_field_module() ; XU
            Free module X(U) of vector fields on the Open subset U of
             the 2-dimensional differentiable manifold S^2
            sage: XU.category()
            Category of finite dimensional modules over Algebra of
             differentiable scalar fields on the Open subset U of
             the 2-dimensional differentiable manifold S^2
            sage: XU.base_ring()
            Algebra of differentiable scalar fields on the Open subset U of
             the 2-dimensional differentiable manifold S^2
            sage: XU.base_ring() is U.scalar_field_algebra()
            True

        `\mathfrak{X}(U)` is a free module because `U` is parallelizable
        (being a chart domain)::

            sage: U.is_manifestly_parallelizable()
            True

        Its rank is the manifold's dimension::

            sage: XU.rank()
            2

        The elements of `\mathfrak{X}(U)` are vector fields on `U`::

            sage: XU.an_element()
            Vector field on the Open subset U of the 2-dimensional
             differentiable manifold S^2
            sage: XU.an_element().display()
            2 ∂/∂th + 2 ∂/∂ph

        Vector field module `\mathfrak{X}(U,\Phi)` of the
        `\RR^3`-valued vector fields along `U`, associated with the
        embedding `\Phi` of `\mathbb{S}^2` into `\RR^3`::

            sage: R3 = Manifold(3, 'R^3')
            sage: cart_coord.<x, y, z> = R3.chart()
            sage: Phi = U.diff_map(R3,
            ....:      [sin(th)*cos(ph), sin(th)*sin(ph), cos(th)], name='Phi')
            sage: XU_R3 = U.vector_field_module(dest_map=Phi) ; XU_R3
            Free module X(U,Phi) of vector fields along the Open subset U of
             the 2-dimensional differentiable manifold S^2 mapped into the
             3-dimensional differentiable manifold R^3
            sage: XU_R3.base_ring()
            Algebra of differentiable scalar fields on the Open subset U of the
             2-dimensional differentiable manifold S^2

        `\mathfrak{X}(U,\Phi)` is a free module because `\RR^3`
        is parallelizable and its rank is 3::

            sage: XU_R3.rank()
            3

        Without any information on the manifold, the vector field module is
        not free by default::

            sage: M = Manifold(2, 'M')
            sage: XM = M.vector_field_module()
            sage: isinstance(XM, FiniteRankFreeModule)
            False

        In particular, declaring a coordinate chart on ``M`` would yield an
        error::

            sage: X.<x,y> = M.chart()
            Traceback (most recent call last):
            ...
            ValueError: the Module X(M) of vector fields on the 2-dimensional
             differentiable manifold M has already been constructed as a
             non-free module, which implies that the 2-dimensional
             differentiable manifold M is not parallelizable and hence cannot
             be the domain of a coordinate chart

        Similarly, one cannot declare a vector frame on `M`::

            sage: e = M.vector_frame('e')
            Traceback (most recent call last):
            ...
            ValueError: the Module X(M) of vector fields on the 2-dimensional
             differentiable manifold M has already been constructed as a
             non-free module and therefore cannot have a basis

        One shall use the keyword ``force_free=True`` to construct a free
        module before declaring the chart::

            sage: M = Manifold(2, 'M')
            sage: XM = M.vector_field_module(force_free=True)
            sage: X.<x,y> = M.chart()  # OK
            sage: e = M.vector_frame('e')  # OK

        If one declares the chart or the vector frame before asking for the
        vector field module, the latter is initialized as a free module,
        without the need to specify ``force_free=True``. Indeed, the
        information that `M` is the domain of a chart or a vector frame implies
        that `M` is parallelizable and is therefore sufficient to assert that
        `\mathfrak{X}(M)` is a free module over `C^k(M)`::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: XM = M.vector_field_module()
            sage: isinstance(XM, FiniteRankFreeModule)
            True
            sage: M.is_manifestly_parallelizable()
            True
        """
        from sage.manifolds.differentiable.vectorfield_module import (
            VectorFieldFreeModule,
            VectorFieldModule,
        )
        if dest_map is None:
            dest_map = self.identity_map()
        codomain = dest_map._codomain
        if dest_map not in self._vector_field_modules:
            if codomain.is_manifestly_parallelizable() or force_free:
                self._vector_field_modules[dest_map] = \
                                 VectorFieldFreeModule(self, dest_map=dest_map)
            else:
                self._vector_field_modules[dest_map] = \
                                     VectorFieldModule(self, dest_map=dest_map)
        return self._vector_field_modules[dest_map]

    def tensor_field_module(self, tensor_type, dest_map=None):
        r"""
        Return the set of tensor fields of a given type defined on ``self``,
        possibly with values in another manifold, as a module over
        the algebra of scalar fields defined on ``self``.

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.tensorfield_module.TensorFieldModule`
            for a complete documentation.

        INPUT:

        - ``tensor_type`` -- pair `(k,l)` with `k` being the contravariant
          rank and `l` the covariant rank
        - ``dest_map`` -- (default: ``None``) destination map, i.e. a
          differentiable map `\Phi:\ M \rightarrow N`, where `M` is the
          current manifold and `N` a differentiable manifold;
          if ``None``, it is assumed that `N = M` and that `\Phi` is the
          identity map (case of tensor fields *on* `M`), otherwise
          ``dest_map`` must be a
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`

        OUTPUT:

        - a
          :class:`~sage.manifolds.differentiable.tensorfield_module.TensorFieldModule`
          (or if `N` is parallelizable, a
          :class:`~sage.manifolds.differentiable.tensorfield_module.TensorFieldFreeModule`)
          representing the module `\mathcal{T}^{(k,l)}(M,\Phi)` of type-`(k,l)`
          tensor fields on `M` taking values on `\Phi(M)\subset N`

        EXAMPLES:

        Module of type-`(2,1)` tensor fields on a 3-dimensional open subset of
        a differentiable manifold::

            sage: M = Manifold(3, 'M')
            sage: U = M.open_subset('U')
            sage: c_xyz.<x,y,z> = U.chart()
            sage: TU = U.tensor_field_module((2,1)) ; TU
            Free module T^(2,1)(U) of type-(2,1) tensors fields on the Open
             subset U of the 3-dimensional differentiable manifold M
            sage: TU.category()
            Category of tensor products of finite dimensional modules
             over Algebra of differentiable scalar fields
              on the Open subset U of the 3-dimensional differentiable manifold M
            sage: TU.base_ring()
            Algebra of differentiable scalar fields on the Open subset U of
             the 3-dimensional differentiable manifold M
            sage: TU.base_ring() is U.scalar_field_algebra()
            True
            sage: TU.an_element()
            Tensor field of type (2,1) on the Open subset U of the
             3-dimensional differentiable manifold M
            sage: TU.an_element().display()
            2 ∂/∂x⊗∂/∂x⊗dx
        """
        return self.vector_field_module(dest_map=dest_map).tensor_module(*tensor_type)

    def diff_form_module(self, degree, dest_map=None):
        r"""
        Return the set of differential forms of a given degree defined on
        ``self``, possibly with values in another manifold, as a module
        over the algebra of scalar fields defined on ``self``.

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.diff_form_module.DiffFormModule`
            for complete documentation.

        INPUT:

        - ``degree`` -- positive integer; the degree `p` of the
          differential forms
        - ``dest_map`` -- (default: ``None``) destination map, i.e. a
          differentiable map `\Phi:\ M \rightarrow N`, where `M` is the
          current manifold and `N` a differentiable manifold;
          if ``None``, it is assumed that `N = M` and that `\Phi` is the
          identity map (case of differential forms *on* `M`), otherwise
          ``dest_map`` must be a
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`

        OUTPUT:

        - a
          :class:`~sage.manifolds.differentiable.diff_form_module.DiffFormModule`
          (or if `N` is parallelizable, a
          :class:`~sage.manifolds.differentiable.diff_form_module.DiffFormFreeModule`)
          representing the module `\Omega^p(M,\Phi)` of `p`-forms on `M`
          taking values on `\Phi(M)\subset N`

        EXAMPLES:

        Module of 2-forms on a 3-dimensional parallelizable manifold::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: M.diff_form_module(2)
            Free module Omega^2(M) of 2-forms on the 3-dimensional
             differentiable manifold M
            sage: M.diff_form_module(2).category()
            Category of finite dimensional modules over Algebra of
             differentiable scalar fields on the 3-dimensional
             differentiable manifold M
            sage: M.diff_form_module(2).base_ring()
            Algebra of differentiable scalar fields on the 3-dimensional
             differentiable manifold M
            sage: M.diff_form_module(2).rank()
            3

        The outcome is cached::

            sage: M.diff_form_module(2) is M.diff_form_module(2)
            True
        """
        return self.vector_field_module(dest_map=dest_map).dual_exterior_power(degree)

    def mixed_form_algebra(self, dest_map=None):
        r"""
        Return the set of mixed forms defined on ``self``, possibly with values
        in another manifold, as a graded algebra.

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.mixed_form_algebra.MixedFormAlgebra`
            for complete documentation.

        INPUT:

        - ``dest_map`` -- (default: ``None``) destination map, i.e. a
          differentiable map `\Phi:\ M \rightarrow N`, where `M` is the
          current manifold and `N` a differentiable manifold;
          if ``None``, it is assumed that `N = M` and that `\Phi` is the
          identity map (case of mixed forms *on* `M`), otherwise
          ``dest_map`` must be a
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`

        OUTPUT:

        - a
          :class:`~sage.manifolds.differentiable.mixed_form_algebra.MixedFormAlgebra`
          representing the graded algebra `\Omega^*(M,\Phi)` of mixed forms on `M`
          taking values on `\Phi(M)\subset N`

        EXAMPLES:

        Graded algebra of mixed forms on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: M.mixed_form_algebra()
            Graded algebra Omega^*(M) of mixed differential forms on the
             2-dimensional differentiable manifold M
            sage: M.mixed_form_algebra().category()
            Join of Category of graded algebras over Symbolic Ring and Category of chain complexes over Symbolic Ring
            sage: M.mixed_form_algebra().base_ring()
            Symbolic Ring

        The outcome is cached::

            sage: M.mixed_form_algebra() is M.mixed_form_algebra()
            True
        """
        vmodule = self.vector_field_module(dest_map=dest_map)
        return MixedFormAlgebra(vmodule)

    de_rham_complex = mixed_form_algebra

    def multivector_module(self, degree, dest_map=None):
        r"""
        Return the set of multivector fields of a given degree defined
        on ``self``, possibly with values in another manifold, as a
        module over the algebra of scalar fields defined on ``self``.

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.multivector_module.MultivectorModule`
            for complete documentation.

        INPUT:

        - ``degree`` -- positive integer; the degree `p` of the
          multivector fields
        - ``dest_map`` -- (default: ``None``) destination map, i.e. a
          differentiable map `\Phi:\ M \rightarrow N`, where `M` is the
          current manifold and `N` a differentiable manifold;
          if ``None``, it is assumed that `N = M` and that `\Phi` is the
          identity map (case of multivector fields *on* `M`), otherwise
          ``dest_map`` must be a
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`

        OUTPUT:

        - a
          :class:`~sage.manifolds.differentiable.multivector_module.MultivectorModule`
          (or if `N` is parallelizable, a
          :class:`~sage.manifolds.differentiable.multivector_module.MultivectorFreeModule`)
          representing the module `\Omega^p(M,\Phi)` of `p`-forms on `M`
          taking values on `\Phi(M)\subset N`

        EXAMPLES:

        Module of 2-vector fields on a 3-dimensional parallelizable
        manifold::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: M.multivector_module(2)
            Free module A^2(M) of 2-vector fields on the 3-dimensional
             differentiable manifold M
            sage: M.multivector_module(2).category()
            Category of finite dimensional modules over Algebra of
             differentiable scalar fields on the 3-dimensional
             differentiable manifold M
            sage: M.multivector_module(2).base_ring()
            Algebra of differentiable scalar fields on the 3-dimensional
             differentiable manifold M
            sage: M.multivector_module(2).rank()
            3

        The outcome is cached::

            sage: M.multivector_module(2) is M.multivector_module(2)
            True
        """
        return self.vector_field_module(dest_map=dest_map).exterior_power(degree)

    def automorphism_field_group(self, dest_map=None):
        r"""
        Return the group of tangent-space automorphism fields defined on
        ``self``, possibly with values in another manifold, as a module
        over the algebra of scalar fields defined on ``self``.

        If `M` is the current manifold and `\Phi` a differentiable map
        `\Phi: M \rightarrow N`, where `N` is a differentiable manifold,
        this method called with ``dest_map`` being `\Phi` returns the
        general linear group `\mathrm{GL}(\mathfrak{X}(M, \Phi))` of the module
        `\mathfrak{X}(M, \Phi)` of vector fields along `M` with values in
        `\Phi(M) \subset N`.

        INPUT:

        - ``dest_map`` -- (default: ``None``) destination map, i.e. a
          differentiable map `\Phi:\ M \rightarrow N`, where `M` is the
          current manifold and `N` a differentiable manifold;
          if ``None``, it is assumed that `N = M` and that `\Phi` is the
          identity map, otherwise ``dest_map`` must be a
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`

        OUTPUT:

        - a
          :class:`~sage.manifolds.differentiable.automorphismfield_group.AutomorphismFieldParalGroup`
          (if `N` is parallelizable) or a
          :class:`~sage.manifolds.differentiable.automorphismfield_group.AutomorphismFieldGroup`
          (if `N` is not parallelizable) representing
          `\mathrm{GL}(\mathfrak{X}(U, \Phi))`

        EXAMPLES:

        Group of tangent-space automorphism fields of a 2-dimensional
        differentiable manifold::

            sage: M = Manifold(2, 'M')
            sage: M.automorphism_field_group()
            General linear group of the Module X(M) of vector fields on the
             2-dimensional differentiable manifold M
            sage: M.automorphism_field_group().category()
            Category of groups

        .. SEEALSO::

            For more examples, see
            :class:`~sage.manifolds.differentiable.automorphismfield_group.AutomorphismFieldParalGroup`
            and
            :class:`~sage.manifolds.differentiable.automorphismfield_group.AutomorphismFieldGroup`.
        """
        return self.vector_field_module(dest_map=dest_map).general_linear_group()

    def vector_field(self, *comp, **kwargs):
        r"""
        Define a vector field on ``self``.

        Via the argument ``dest_map``, it is possible to let the vector field
        take its values on another manifold. More precisely, if `M` is
        the current manifold, `N` a differentiable manifold and
        `\Phi:\  M \rightarrow N` a differentiable map, a *vector field
        along* `M` *with values on* `N` is a differentiable map

        .. MATH::

            v:\ M  \longrightarrow TN

        (`TN` being the tangent bundle of `N`) such that

        .. MATH::

            \forall p \in M,\ v(p) \in T_{\Phi(p)} N,

        where `T_{\Phi(p)} N` is the tangent space to `N` at the
        point `\Phi(p)`.

        The standard case of vector fields *on* `M` corresponds
        to `N = M` and `\Phi = \mathrm{Id}_M`. Other common cases are `\Phi`
        being an immersion and `\Phi` being a curve in `N` (`M` is then
        an open interval of `\RR`).

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.vectorfield.VectorField`
            and
            :class:`~sage.manifolds.differentiable.vectorfield.VectorFieldParal`
            for a complete documentation.

        INPUT:

        - ``comp`` -- (optional) either the components of the vector field
          with respect to the vector frame specified by the argument ``frame``
          or a dictionary of components, the keys of which are vector frames or
          pairs ``(f, c)`` where ``f`` is a vector frame and ``c`` the chart
          in which the components are expressed
        - ``frame`` -- (default: ``None``; unused if ``comp`` is not given or
          is a dictionary) vector frame in which the components are given; if
          ``None``, the default vector frame of ``self`` is assumed
        - ``chart`` -- (default: ``None``; unused if ``comp`` is not given or
          is a dictionary) coordinate chart in which the components are
          expressed; if ``None``, the default chart on the domain of ``frame``
          is assumed
        - ``name`` -- (default: ``None``) name given to the vector field
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          vector field; if none is provided, the LaTeX symbol is set to
          ``name``
        - ``dest_map`` -- (default: ``None``) the destination map
          `\Phi:\ M \rightarrow N`; if ``None``, it is assumed that `N = M`
          and that `\Phi` is the identity map (case of a vector field
          *on* `M`), otherwise ``dest_map`` must be a
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`

        OUTPUT:

        - a
          :class:`~sage.manifolds.differentiable.vectorfield.VectorField`
          (or if `N` is parallelizable, a
          :class:`~sage.manifolds.differentiable.vectorfield.VectorFieldParal`)
          representing the defined vector field

        EXAMPLES:

        A vector field on a open subset of a 3-dimensional differentiable
        manifold::

            sage: M = Manifold(3, 'M')
            sage: U = M.open_subset('U')
            sage: c_xyz.<x,y,z> = U.chart()
            sage: v = U.vector_field(y, -x*z, 1+y, name='v'); v
            Vector field v on the Open subset U of the 3-dimensional
             differentiable manifold M
            sage: v.display()
            v = y ∂/∂x - x*z ∂/∂y + (y + 1) ∂/∂z

        The vector fields on `U` form the set `\mathfrak{X}(U)`, which is a
        module over the algebra `C^k(U)` of differentiable scalar fields
        on `U`::

            sage: v.parent()
            Free module X(U) of vector fields on the Open subset U of the
             3-dimensional differentiable manifold M
            sage: v in U.vector_field_module()
            True

        For more examples, see
        :class:`~sage.manifolds.differentiable.vectorfield.VectorField` and
        :class:`~sage.manifolds.differentiable.vectorfield.VectorFieldParal`.
        """
        name = kwargs.pop('name', None)
        latex_name = kwargs.pop('latex_name', None)
        dest_map = kwargs.pop('dest_map', None)
        vmodule = self.vector_field_module(dest_map)  # the parent
        resu = vmodule.element_class(vmodule, name=name, latex_name=latex_name)
        if comp:
            # Some components are to be initialized
            resu._init_components(*comp, **kwargs)
        return resu

    def tensor_field(self, *args, **kwargs):
        r"""
        Define a tensor field on ``self``.

        Via the argument ``dest_map``, it is possible to let the tensor field
        take its values on another manifold. More precisely, if `M` is
        the current manifold, `N` a differentiable manifold,
        `\Phi:\  M \rightarrow N` a differentiable map and `(k,l)`
        a pair of nonnegative integers, a *tensor field of type* `(k,l)`
        *along* `M` *with values on* `N` is a differentiable map

        .. MATH::

            t:\ M  \longrightarrow T^{(k,l)} N

        (`T^{(k,l)}N` being the tensor bundle of type `(k,l)` over `N`)
        such that

        .. MATH::

            \forall p \in M,\ t(p) \in T^{(k,l)}(T_{\Phi(p)} N),

        where `T^{(k,l)}(T_{\Phi(p)} N)` is the space of tensors of type
        `(k,l)` on the tangent space `T_{\Phi(p)} N`.

        The standard case of tensor fields *on* `M` corresponds
        to `N=M` and `\Phi = \mathrm{Id}_M`. Other common cases are `\Phi`
        being an immersion and `\Phi` being a curve in `N` (`M` is then
        an open interval of `\RR`).

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.tensorfield.TensorField`
            and
            :class:`~sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal`
            for a complete documentation.

        INPUT:

        - ``k`` -- the contravariant rank `k`, the tensor type being `(k,l)`
        - ``l`` -- the covariant rank `l`, the tensor type being `(k,l)`
        - ``comp`` -- (optional) either the components of the tensor field
          with respect to the vector frame specified by the argument
          ``frame`` or a dictionary of components, the keys of which are vector
          frames or pairs ``(f, c)`` where ``f`` is a vector frame and ``c``
          the chart in which the components are expressed
        - ``frame`` -- (default: ``None``; unused if ``comp`` is not given or
          is a dictionary) vector frame in which the components are given; if
          ``None``, the default vector frame of ``self`` is assumed
        - ``chart`` -- (default: ``None``; unused if ``comp`` is not given or
          is a dictionary) coordinate chart in which the components are
          expressed; if ``None``, the default chart on the domain of ``frame``
          is assumed
        - ``name`` -- (default: ``None``) name given to the tensor field
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          tensor field; if ``None``, the LaTeX symbol is set to ``name``
        - ``sym`` -- (default: ``None``) a symmetry or a list of symmetries
          among the tensor arguments: each symmetry is described by a tuple
          containing the positions of the involved arguments, with the
          convention ``position=0`` for the first argument; for instance:

          * ``sym = (0,1)`` for a symmetry between the 1st and 2nd arguments
          * ``sym = [(0,2), (1,3,4)]`` for a symmetry between the 1st and 3rd
            arguments and a symmetry between the 2nd, 4th and 5th arguments

        - ``antisym`` -- (default: ``None``) antisymmetry or list of
          antisymmetries among the arguments, with the same convention as for
          ``sym``
        - ``dest_map`` -- (default: ``None``) the destination map
          `\Phi:\ M \rightarrow N`; if ``None``, it is assumed that `N = M`
          and that `\Phi` is the identity map (case of a tensor field
          *on* `M`), otherwise ``dest_map`` must be a
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`

        OUTPUT:

        - a :class:`~sage.manifolds.differentiable.tensorfield.TensorField`
          (or if `N` is parallelizable, a
          :class:`~sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal`)
          representing the defined tensor field

        EXAMPLES:

        A tensor field of type `(2,0)` on a 2-dimensional differentiable
        manifold::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: t = M.tensor_field(2, 0, [[1+x, -y], [0, x*y]], name='T'); t
            Tensor field T of type (2,0) on the 2-dimensional differentiable
             manifold M
            sage: t.display()
            T = (x + 1) ∂/∂x⊗∂/∂x - y ∂/∂x⊗∂/∂y + x*y ∂/∂y⊗∂/∂y

        The type `(2,0)` tensor fields on `M` form the set
        `\mathcal{T}^{(2,0)}(M)`, which is a module over the algebra `C^k(M)`
        of differentiable scalar fields on `M`::

            sage: t.parent()
            Free module T^(2,0)(M) of type-(2,0) tensors fields on the
             2-dimensional differentiable manifold M
            sage: t in M.tensor_field_module((2,0))
            True

        For more examples, see
        :class:`~sage.manifolds.differentiable.tensorfield.TensorField` and
        :class:`~sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal`.
        """
        k = args[0]
        l = args[1]
        name = kwargs.pop('name', None)
        latex_name = kwargs.pop('latex_name', None)
        sym = kwargs.pop('sym', None)
        antisym = kwargs.pop('antisym', None)
        dest_map = kwargs.pop('dest_map', None)
        vmodule = self.vector_field_module(dest_map)
        resu = vmodule.tensor((k, l), name=name, latex_name=latex_name,
                              sym=sym, antisym=antisym)
        if len(args) > 2:
            # Some components are to be initialized
            resu._init_components(args[2], **kwargs)
        return resu

    def sym_bilin_form_field(self, *comp, **kwargs):
        r"""
        Define a field of symmetric bilinear forms on ``self``.

        Via the argument ``dest_map``, it is possible to let the field
        take its values on another manifold. More precisely, if `M` is
        the current manifold, `N` a differentiable manifold and
        `\Phi:\  M \rightarrow N` a differentiable map, a *field of
        symmetric bilinear forms along* `M` *with values on* `N` is a
        differentiable map

        .. MATH::

            t:\ M  \longrightarrow T^{(0,2)}N

        (`T^{(0,2)} N` being the tensor bundle of type `(0,2)` over `N`)
        such that

        .. MATH::

            \forall p \in M,\ t(p) \in S(T_{\Phi(p)} N),

        where `S(T_{\Phi(p)} N)` is the space of symmetric bilinear forms on
        the tangent space `T_{\Phi(p)} N`.

        The standard case of fields of symmetric bilinear forms *on* `M`
        corresponds to `N = M` and `\Phi = \mathrm{Id}_M`. Other common
        cases are `\Phi` being an immersion and `\Phi` being a curve in `N`
        (`M` is then an open interval of `\RR`).

        INPUT:

        - ``comp`` -- (optional) either the components of the field of
          symmetric bilinear forms with respect to the vector frame specified
          by the argument ``frame`` or a dictionary of components, the keys of
          which are vector frames or pairs ``(f, c)`` where ``f`` is a vector
          frame and ``c`` the chart in which the components are expressed
        - ``frame`` -- (default: ``None``; unused if ``comp`` is not given or
          is a dictionary) vector frame in which the components are given; if
          ``None``, the default vector frame of ``self`` is assumed
        - ``chart`` -- (default: ``None``; unused if ``comp`` is not given or
          is a dictionary) coordinate chart in which the components are
          expressed; if ``None``, the default chart on the domain of ``frame``
          is assumed
        - ``name`` -- (default: ``None``) name given to the field
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          field; if none is provided, the LaTeX symbol is set to ``name``
        - ``dest_map`` -- (default: ``None``) the destination map
          `\Phi:\ M \rightarrow N`; if ``None``, it is assumed that `N = M`
          and that `\Phi` is the identity map (case of a field *on* `M`),
          otherwise ``dest_map`` must be an instance of instance of
          class :class:`~sage.manifolds.differentiable.diff_map.DiffMap`

        OUTPUT:

        - a :class:`~sage.manifolds.differentiable.tensorfield.TensorField`
          (or if `N` is parallelizable, a
          :class:`~sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal`)
          of tensor type `(0,2)` and symmetric representing the defined
          field of symmetric bilinear forms

        EXAMPLES:

        A field of symmetric bilinear forms on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: t = M.sym_bilin_form_field(name='T'); t
            Field of symmetric bilinear forms T on the 2-dimensional
             differentiable manifold M

        Such a object is a tensor field of rank 2 and type `(0,2)`::

            sage: t.parent()
            Free module T^(0,2)(M) of type-(0,2) tensors fields on the
             2-dimensional differentiable manifold M
            sage: t.tensor_rank()
            2
            sage: t.tensor_type()
            (0, 2)

        The LaTeX symbol is deduced from the name or can be specified when
        creating the object::

            sage: latex(t)
            T
            sage: om = M.sym_bilin_form_field(name='Omega', latex_name=r'\Omega')
            sage: latex(om)
            \Omega

        Setting the components in the manifold's default vector frame::

            sage: t[0,0], t[0,1], t[1,1] = -1, x, x*y

        The unset components are either zero or deduced by symmetry::

            sage: t[1, 0]
            x
            sage: t[:]
            [ -1   x]
            [  x x*y]

        One can also set the components while defining the field of symmetric
        bilinear forms::

            sage: t = M.sym_bilin_form_field([[-1, x], [x, x*y]], name='T')

        A symmetric bilinear form acts on vector pairs::

            sage: v1 = M.vector_field(y, x, name='V_1')
            sage: v2 = M.vector_field(x+y, 2, name='V_2')
            sage: s = t(v1,v2) ; s
            Scalar field T(V_1,V_2) on the 2-dimensional differentiable
             manifold M
            sage: s.expr()
            x^3 + (3*x^2 + x)*y - y^2
            sage: s.expr() - t[0,0]*v1[0]*v2[0] - \
            ....: t[0,1]*(v1[0]*v2[1]+v1[1]*v2[0]) - t[1,1]*v1[1]*v2[1]
            0
            sage: latex(s)
            T\left(V_1,V_2\right)

        Adding two symmetric bilinear forms results in another symmetric
        bilinear form::

            sage: a = M.sym_bilin_form_field([[1, 2], [2, 3]])
            sage: b = M.sym_bilin_form_field([[-1, 4], [4, 5]])
            sage: s = a + b ; s
            Field of symmetric bilinear forms on the 2-dimensional
             differentiable manifold M
            sage: s[:]
            [0 6]
            [6 8]

        But adding a symmetric bilinear from with a non-symmetric bilinear
        form results in a generic type `(0,2)` tensor::

            sage: c = M.tensor_field(0, 2, [[-2, -3], [1,7]])
            sage: s1 = a + c ; s1
            Tensor field of type (0,2) on the 2-dimensional differentiable
             manifold M
            sage: s1[:]
            [-1 -1]
            [ 3 10]
            sage: s2 = c + a ; s2
            Tensor field of type (0,2) on the 2-dimensional differentiable
             manifold M
            sage: s2[:]
            [-1 -1]
            [ 3 10]
        """
        name = kwargs.pop('name', None)
        latex_name = kwargs.pop('latex_name', None)
        dest_map = kwargs.pop('dest_map', None)
        vmodule = self.vector_field_module(dest_map)
        resu = vmodule.tensor((0, 2), name=name, latex_name=latex_name,
                              sym=(0,1))
        if comp:
            # Some components are to be initialized
            resu._init_components(*comp, **kwargs)
        return resu

    def multivector_field(self, *args, **kwargs):
        r"""
        Define a multivector field on ``self``.

        Via the argument ``dest_map``, it is possible to let the
        multivector field take its values on another manifold. More
        precisely, if `M` is the current manifold, `N` a differentiable
        manifold, `\Phi:\  M \rightarrow N` a differentiable map and `p`
        a nonnegative integer, a *multivector field of degree* `p` (or
        `p`-*vector field*) *along* `M` *with values on* `N` is a
        differentiable map

        .. MATH::

            t:\ M  \longrightarrow T^{(p,0)} N

        (`T^{(p,0)} N` being the tensor bundle of type `(p,0)` over `N`)
        such that

        .. MATH::

            \forall x \in M,\quad t(x) \in \Lambda^p(T_{\Phi(x)} N),

        where `\Lambda^p(T_{\Phi(x)} N)` is the `p`-th exterior power
        of the tangent vector space `T_{\Phi(x)} N`.

        The standard case of a `p`-vector field *on* `M` corresponds
        to `N = M` and `\Phi = \mathrm{Id}_M`. Other common cases are
        `\Phi` being an immersion and `\Phi` being a curve in `N` (`M`
        is then an open interval of `\RR`).

        For `p = 1`, one can use the method
        :meth:`~sage.manifolds.differentiable.manifold.DifferentiableManifold.vector_field`
        instead.

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.multivectorfield.MultivectorField`
            and
            :class:`~sage.manifolds.differentiable.multivectorfield.MultivectorFieldParal`
            for a complete documentation.

        INPUT:

        - ``degree`` -- the degree `p` of the multivector field (i.e.
          its tensor rank)
        - ``comp`` -- (optional) either the components of the multivector field
          with respect to the vector frame specified by the argument ``frame``
          or a dictionary of components, the keys of which are vector frames
          or pairs ``(f, c)`` where ``f`` is a vector frame and ``c`` the chart
          in which the components are expressed
        - ``frame`` -- (default: ``None``; unused if ``comp`` is not given or
          is a dictionary) vector frame in which the components are given; if
          ``None``, the default vector frame of ``self`` is assumed
        - ``chart`` -- (default: ``None``; unused if ``comp`` is not given or
          is a dictionary) coordinate chart in which the components are
          expressed; if ``None``, the default chart on the domain of ``frame``
          is assumed
        - ``name`` -- (default: ``None``) name given to the multivector
          field
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote
          the multivector field; if none is provided, the LaTeX symbol
          is set to ``name``
        - ``dest_map`` -- (default: ``None``) the destination map
          `\Phi:\ M \rightarrow N`; if ``None``, it is assumed that
          `N = M` and that `\Phi` is the identity map (case of a
          multivector field *on* `M`), otherwise ``dest_map`` must be a
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`

        OUTPUT:

        - the `p`-vector field as a
          :class:`~sage.manifolds.differentiable.multivectorfield.MultivectorField`
          (or if `N` is parallelizable, a
          :class:`~sage.manifolds.differentiable.multivectorfield.MultivectorFieldParal`)

        EXAMPLES:

        A 2-vector field on a 3-dimensional differentiable manifold::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: h = M.multivector_field(2, name='H'); h
            2-vector field H on the 3-dimensional differentiable manifold M
            sage: h[0,1], h[0,2], h[1,2] = x+y, x*z, -3
            sage: h.display()
            H = (x + y) ∂/∂x∧∂/∂y + x*z ∂/∂x∧∂/∂z - 3 ∂/∂y∧∂/∂z

        For more examples, see
        :class:`~sage.manifolds.differentiable.multivectorfield.MultivectorField`
        and
        :class:`~sage.manifolds.differentiable.multivectorfield.MultivectorFieldParal`.
        """
        degree = args[0]
        name = kwargs.pop('name', None)
        latex_name = kwargs.pop('latex_name', None)
        dest_map = kwargs.pop('dest_map', None)
        vmodule = self.vector_field_module(dest_map)
        resu = vmodule.alternating_contravariant_tensor(degree, name=name,
                                                        latex_name=latex_name)
        if len(args) > 1:
            # Some components are to be initialized
            resu._init_components(args[1], **kwargs)
        return resu

    def diff_form(self, *args, **kwargs) -> DiffForm:
        r"""
        Define a differential form on ``self``.

        Via the argument ``dest_map``, it is possible to let the
        differential form take its values on another manifold. More
        precisely, if `M` is the current manifold, `N` a differentiable
        manifold, `\Phi:\  M \rightarrow N` a differentiable map and `p`
        a nonnegative integer, a *differential form of degree* `p` (or
        `p`-*form*) *along* `M` *with values on* `N` is a differentiable
        map

        .. MATH::

            t:\ M  \longrightarrow T^{(0,p)}N

        (`T^{(0,p)} N` being the tensor bundle of type `(0,p)` over `N`)
        such that

        .. MATH::

            \forall x \in M,\quad t(x) \in \Lambda^p(T^*_{\Phi(x)} N),

        where `\Lambda^p(T^*_{\Phi(x)} N)` is the `p`-th exterior power
        of the dual of the tangent space `T_{\Phi(x)} N`.

        The standard case of a differential form *on* `M` corresponds
        to `N = M` and `\Phi = \mathrm{Id}_M`. Other common cases are
        `\Phi` being an immersion and `\Phi` being a curve in `N` (`M`
        is then an open interval of `\RR`).

        For `p = 1`, one can use the method
        :meth:`~sage.manifolds.differentiable.manifold.DifferentiableManifold.one_form`
        instead.

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.diff_form.DiffForm` and
            :class:`~sage.manifolds.differentiable.diff_form.DiffFormParal`
            for a complete documentation.

        INPUT:

        - ``degree`` -- the degree `p` of the differential form (i.e.
          its tensor rank)
        - ``comp`` -- (optional) either the components of the differential
          form with respect to the vector frame specified by the argument
          ``frame`` or a dictionary of components, the keys of which are vector
          frames or pairs ``(f, c)`` where ``f`` is a vector frame and ``c``
          the chart in which the components are expressed
        - ``frame`` -- (default: ``None``; unused if ``comp`` is not given or
          is a dictionary) vector frame in which the components are given; if
          ``None``, the default vector frame of ``self`` is assumed
        - ``chart`` -- (default: ``None``; unused if ``comp`` is not given or
          is a dictionary) coordinate chart in which the components are
          expressed; if ``None``, the default chart on the domain of ``frame``
          is assumed
        - ``name`` -- (default: ``None``) name given to the differential
          form
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote
          the differential form; if none is provided, the LaTeX symbol
          is set to ``name``
        - ``dest_map`` -- (default: ``None``) the destination map
          `\Phi:\ M \rightarrow N`; if ``None``, it is assumed that
          `N = M` and that `\Phi` is the identity map (case of a
          differential form *on* `M`), otherwise ``dest_map`` must be a
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`

        OUTPUT:

        - the `p`-form as a
          :class:`~sage.manifolds.differentiable.diff_form.DiffForm`
          (or if `N` is parallelizable, a
          :class:`~sage.manifolds.differentiable.diff_form.DiffFormParal`)

        EXAMPLES:

        A 2-form on a 3-dimensional differentiable manifold::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: f = M.diff_form(2, name='F'); f
            2-form F on the 3-dimensional differentiable manifold M
            sage: f[0,1], f[1,2] = x+y, x*z
            sage: f.display()
            F = (x + y) dx∧dy + x*z dy∧dz

        For more examples, see
        :class:`~sage.manifolds.differentiable.diff_form.DiffForm` and
        :class:`~sage.manifolds.differentiable.diff_form.DiffFormParal`.
        """
        degree = args[0]
        name = kwargs.pop('name', None)
        latex_name = kwargs.pop('latex_name', None)
        dest_map = kwargs.pop('dest_map', None)
        vmodule = self.vector_field_module(dest_map)
        resu = vmodule.alternating_form(degree, name=name,
                                        latex_name=latex_name)
        if len(args) > 1:
            # Some components are to be initialized
            resu._init_components(args[1], **kwargs)
        return resu

    def one_form(self, *comp, **kwargs) -> DiffForm:
        r"""
        Define a 1-form on the manifold.

        Via the argument ``dest_map``, it is possible to let the
        1-form take its values on another manifold. More precisely,
        if `M` is the current manifold, `N` a differentiable
        manifold and `\Phi:\ M \rightarrow N` a differentiable map,
        a *1-form along* `M` *with values on* `N` is a differentiable
        map

        .. MATH::

            t:\ M  \longrightarrow T^* N

        (`T^* N` being the cotangent bundle of `N`) such that

        .. MATH::

            \forall p \in M,\quad t(p) \in T^*_{\Phi(p)}N,

        where `T^*_{\Phi(p)}` is the dual of the tangent space
        `T_{\Phi(p)} N`.

        The standard case of a 1-form *on* `M` corresponds to `N = M`
        and `\Phi = \mathrm{Id}_M`. Other common cases are `\Phi`
        being an immersion and `\Phi` being a curve in `N` (`M` is then
        an open interval of `\RR`).

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.diff_form.DiffForm` and
            :class:`~sage.manifolds.differentiable.diff_form.DiffFormParal`
            for a complete documentation.

        INPUT:

        - ``comp`` -- (optional) either the components of 1-form with respect
          to the vector frame specified by the argument ``frame`` or a
          dictionary of components, the keys of which are vector frames or
          pairs ``(f, c)`` where ``f`` is a vector frame and ``c`` the chart
          in which the components are expressed
        - ``frame`` -- (default: ``None``; unused if ``comp`` is not given or
          is a dictionary) vector frame in which the components are given; if
          ``None``, the default vector frame of ``self`` is assumed
        - ``chart`` -- (default: ``None``; unused if ``comp`` is not given or
          is a dictionary) coordinate chart in which the components are
          expressed; if ``None``, the default chart on the domain of ``frame``
          is assumed
        - ``name`` -- (default: ``None``) name given to the 1-form
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote
          the 1-form; if none is provided, the LaTeX symbol is set to
          ``name``
        - ``dest_map`` -- (default: ``None``) the destination map
          `\Phi:\ M \rightarrow N`; if ``None``, it is assumed that
          `N = M` and that `\Phi` is the identity map (case of a 1-form
          *on* `M`), otherwise ``dest_map`` must be a
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`

        OUTPUT:

        - the 1-form as a
          :class:`~sage.manifolds.differentiable.diff_form.DiffForm`
          (or if `N` is parallelizable, a
          :class:`~sage.manifolds.differentiable.diff_form.DiffFormParal`)

        EXAMPLES:

        A 1-form on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: om = M.one_form(-y, 2+x, name='omega', latex_name=r'\omega')
            sage: om
            1-form omega on the 2-dimensional differentiable manifold M
            sage: om.display()
            omega = -y dx + (x + 2) dy
            sage: om.parent()
            Free module Omega^1(M) of 1-forms on the 2-dimensional
             differentiable manifold M

        For more examples, see
        :class:`~sage.manifolds.differentiable.diff_form.DiffForm` and
        :class:`~sage.manifolds.differentiable.diff_form.DiffFormParal`.
        """
        name = kwargs.pop('name', None)
        latex_name = kwargs.pop('latex_name', None)
        dest_map = kwargs.pop('dest_map', None)
        vmodule = self.vector_field_module(dest_map)
        resu = vmodule.linear_form(name=name, latex_name=latex_name)
        if comp:
            # Some components are to be initialized
            resu._init_components(*comp, **kwargs)
        return resu

    def mixed_form(self, comp=None, name=None, latex_name=None, dest_map=None):
        r"""
        Define a mixed form on ``self``.

        Via the argument ``dest_map``, it is possible to let the
        mixed form take its values on another manifold. More
        precisely, if `M` is the current manifold, `N` a differentiable
        manifold, `\Phi:\  M \rightarrow N` a differentiable map, a
        *mixed form along* `\Phi` can be considered as a differentiable map

        .. MATH::

            a: M  \longrightarrow \bigoplus^n_{k=0} T^{(0,k)}N

        (`T^{(0,k)} N` being the tensor bundle of type `(0,k)` over `N`, `\oplus`
        being the Whitney sum and `n` being the dimension of `N`) such that

        .. MATH::

            \forall x \in M,\quad a(x) \in \bigoplus^n_{k=0} \Lambda^k(T^*_{\Phi(x)} N),

        where `\Lambda^k(T^*_{\Phi(x)} N)` is the `k`-th exterior power
        of the dual of the tangent space `T_{\Phi(x)} N`.

        The standard case of a mixed form *on* `M` corresponds
        to `N = M` and `\Phi = \mathrm{Id}_M`.

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.mixed_form.MixedForm`
            for complete documentation.

        INPUT:

        - ``comp`` -- (default: ``None``) homogeneous components of the mixed
          form as a list; if none is provided, the components are set to
          innocent unnamed differential forms
        - ``name`` -- (default: ``None``) name given to the differential form
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote
          the differential form; if none is provided, the LaTeX symbol
          is set to ``name``
        - ``dest_map`` -- (default: ``None``) the destination map
          `\Phi:\ M \rightarrow N`; if ``None``, it is assumed that
          `N = M` and that `\Phi` is the identity map (case of a
          differential form *on* `M`), otherwise ``dest_map`` must be a
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`

        OUTPUT:

        - the mixed form as a
          :class:`~sage.manifolds.differentiable.mixed_form.MixedForm`

        EXAMPLES:

        A mixed form on an open subset of a 3-dimensional differentiable
        manifold::

            sage: M = Manifold(3, 'M')
            sage: U = M.open_subset('U', latex_name=r'\mathcal{U}'); U
            Open subset U of the 3-dimensional differentiable manifold M
            sage: c_xyz.<x,y,z> = U.chart()
            sage: f = U.mixed_form(name='F'); f
            Mixed differential form F on the Open subset U of the 3-dimensional
             differentiable manifold M

        See the documentation of class
        :class:`~sage.manifolds.differentiable.mixed_form.MixedForm` for
        more examples.
        """
        algebra = self.mixed_form_algebra(dest_map=dest_map)
        resu = algebra.element_class(algebra, name=name, latex_name=latex_name)
        if comp is not None:
            resu[:] = comp
        return resu

    def symplectic_form(
        self, name: Optional[str] = None, latex_name: Optional[str] = None
    ):
        r"""
        Construct a symplectic form on the current manifold.

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.symplectic_form.SymplecticForm`

        EXAMPLES:

        Standard symplectic form on `\RR^2`::

            sage: M.<q, p> = EuclideanSpace(2)
            sage: omega = M.symplectic_form('omega', r'\omega')
            sage: omega.set_comp()[1,2] = -1
            sage: omega.display()
            omega = -dq∧dp
        """
        return self.vector_field_module().symplectic_form(name, latex_name)

    def poisson_tensor(
        self, name: Optional[str] = None, latex_name: Optional[str] = None
    ):
        r"""
        Construct a Poisson tensor on the current manifold.

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.poisson_tensor.PoissonTensorField`

        EXAMPLES:

        Standard Poisson tensor on `\RR^2`::

            sage: M.<q, p> = EuclideanSpace(2)
            sage: poisson = M.poisson_tensor('varpi')
            sage: poisson.set_comp()[1,2] = -1
            sage: poisson.display()
            varpi = -e_q∧e_p
        """
        return self.vector_field_module().poisson_tensor(name, latex_name)

    def automorphism_field(self, *comp, **kwargs):
        r"""
        Define a field of automorphisms (invertible endomorphisms in each
        tangent space) on ``self``.

        Via the argument ``dest_map``, it is possible to let the
        field take its values on another manifold. More precisely,
        if `M` is the current manifold, `N` a differentiable
        manifold and `\Phi:\  M \rightarrow N` a differentiable map,
        a *field of automorphisms along* `M` *with values on* `N` is a
        differentiable map

        .. MATH::

            t:\ M  \longrightarrow T^{(1,1)} N

        (`T^{(1,1)} N` being the tensor bundle of type `(1,1)` over `N`)
        such that

        .. MATH::

            \forall p \in M,\ t(p) \in \mathrm{GL}\left(T_{\Phi(p)} N \right),

        where `\mathrm{GL}\left(T_{\Phi(p)} N \right)` is the general linear
        group of the tangent space `T_{\Phi(p)} N`.

        The standard case of a field of automorphisms *on* `M` corresponds
        to `N = M` and `\Phi = \mathrm{Id}_M`. Other common cases are `\Phi`
        being an immersion and `\Phi` being a curve in `N` (`M` is then
        an open interval of `\RR`).

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.automorphismfield.AutomorphismField`
            and
            :class:`~sage.manifolds.differentiable.automorphismfield.AutomorphismFieldParal`
            for a complete documentation.

        INPUT:

        - ``comp`` -- (optional) either the components of the field of
          automorphisms with respect to the vector frame specified by the
          argument ``frame`` or a dictionary of components, the keys of which
          are vector frames or pairs ``(f, c)`` where ``f`` is a vector frame
          and ``c`` the chart in which the components are expressed
        - ``frame`` -- (default: ``None``; unused if ``comp`` is not given or
          is a dictionary) vector frame in which the components are given; if
          ``None``, the default vector frame of ``self`` is assumed
        - ``chart`` -- (default: ``None``; unused if ``comp`` is not given or
          is a dictionary) coordinate chart in which the components are
          expressed; if ``None``, the default chart on the domain of ``frame``
          is assumed
        - ``name`` -- (default: ``None``) name given to the field
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          field; if none is provided, the LaTeX symbol is set to ``name``
        - ``dest_map`` -- (default: ``None``) the destination map
          `\Phi:\ M \rightarrow N`; if ``None``, it is assumed that `N = M`
          and that `\Phi` is the identity map (case of a field of
          automorphisms *on* `M`), otherwise ``dest_map`` must be a
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`

        OUTPUT:

        - a
          :class:`~sage.manifolds.differentiable.automorphismfield.AutomorphismField`
          (or if `N` is parallelizable, a
          :class:`~sage.manifolds.differentiable.automorphismfield.AutomorphismFieldParal`)
          representing the defined field of automorphisms

        EXAMPLES:

        A field of automorphisms on a 2-dimensional manifold::

            sage: M = Manifold(2,'M')
            sage: X.<x,y> = M.chart()
            sage: a = M.automorphism_field([[1+x^2, 0], [0, 1+y^2]], name='A')
            sage: a
            Field of tangent-space automorphisms A on the 2-dimensional
             differentiable manifold M
            sage: a.parent()
            General linear group of the Free module X(M) of vector fields on
             the 2-dimensional differentiable manifold M
            sage: a(X.frame()[0]).display()
            A(∂/∂x) = (x^2 + 1) ∂/∂x
            sage: a(X.frame()[1]).display()
            A(∂/∂y) = (y^2 + 1) ∂/∂y

        For more examples, see
        :class:`~sage.manifolds.differentiable.automorphismfield.AutomorphismField`
        and
        :class:`~sage.manifolds.differentiable.automorphismfield.AutomorphismFieldParal`.
        """
        name = kwargs.pop('name', None)
        latex_name = kwargs.pop('latex_name', None)
        dest_map = kwargs.pop('dest_map', None)
        vmodule = self.vector_field_module(dest_map)
        resu = vmodule.automorphism(name=name, latex_name=latex_name)
        if comp:
            # Some components are to be initialized
            resu._init_components(*comp, **kwargs)
        return resu

    def tangent_identity_field(self, dest_map=None):
        r"""
        Return the field of identity maps in the tangent spaces on ``self``.

        Via the argument ``dest_map``, it is possible to let the
        field take its values on another manifold. More precisely,
        if `M` is the current manifold, `N` a differentiable
        manifold and `\Phi:\  M \rightarrow N` a differentiable map,
        a *field of identity maps along* `M` *with values on* `N` is a
        differentiable map

        .. MATH::

            t:\ M \longrightarrow T^{(1,1)} N

        (`T^{(1,1)} N` being the tensor bundle of type `(1,1)` over `N`) such
        that

        .. MATH::

            \forall p \in M,\ t(p) = \mathrm{Id}_{T_{\Phi(p)} N},

        where `\mathrm{Id}_{T_{\Phi(p)} N}` is the identity map of the
        tangent space `T_{\Phi(p)} N`.

        The standard case of a field of identity maps *on* `M` corresponds
        to `N = M` and `\Phi = \mathrm{Id}_M`. Other common cases are `\Phi`
        being an immersion and `\Phi` being a curve in `N` (`M` is then
        an open interval of `\RR`).

        INPUT:

        - ``name`` -- (string; default: 'Id') name given to the field of
          identity maps
        - ``latex_name`` -- (string; default: ``None``) LaTeX symbol to denote
          the field of identity map; if none is provided, the LaTeX symbol is
          set to '\mathrm{Id}' if ``name`` is 'Id' and to ``name`` otherwise
        - ``dest_map`` -- (default: ``None``) the destination map
          `\Phi:\ M \rightarrow N`; if ``None``, it is assumed that `N = M`
          and that `\Phi` is the identity map (case of a field of identity
          maps *on* `M`), otherwise ``dest_map`` must be a
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`

        OUTPUT:

        - a
          :class:`~sage.manifolds.differentiable.automorphismfield.AutomorphismField`
          (or if `N` is parallelizable, a
          :class:`~sage.manifolds.differentiable.automorphismfield.AutomorphismFieldParal`)
          representing the field of identity maps

        EXAMPLES:

        Field of tangent-space identity maps on a 3-dimensional manifold::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = M.chart()
            sage: a = M.tangent_identity_field(); a
            Field of tangent-space identity maps on the 3-dimensional
             differentiable manifold M
            sage: a.comp()
            Kronecker delta of size 3x3

        For more examples, see
        :class:`~sage.manifolds.differentiable.automorphismfield.AutomorphismField`.
        """
        vmodule = self.vector_field_module(dest_map)
        return vmodule.identity_map()

    def set_orientation(self, orientation):
        r"""
        Set the preferred orientation of ``self``.

        INPUT:

        - ``orientation`` -- either a chart / list of charts, or a vector
          frame / list of vector frames, covering ``self``

        .. WARNING::

            It is the user's responsibility that the orientation set here
            is indeed an orientation. There is no check going on in the
            background. See :meth:`orientation` for the definition of an
            orientation.

        EXAMPLES:

        Set an orientation on a manifold::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart(); c_uv.<u,v> = M.chart()
            sage: M.set_orientation(c_uv)
            sage: M.orientation()
            [Coordinate frame (M, (∂/∂u,∂/∂v))]

        Instead of a chart, a vector frame can be given, too::

            sage: M.set_orientation(c_xy.frame())
            sage: M.orientation()
            [Coordinate frame (M, (∂/∂x,∂/∂y))]

        Set an orientation in the non-trivial case::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: M.declare_union(U, V)
            sage: c_xy.<x,y> = U.chart(); c_uv.<u,v> = V.chart()
            sage: M.set_orientation([c_xy, c_uv])
            sage: M.orientation()
            [Coordinate frame (U, (∂/∂x,∂/∂y)),
             Coordinate frame (V, (∂/∂u,∂/∂v))]

        Again, the vector frame notion can be used instead::

            sage: M.set_orientation([c_xy.frame(), c_uv.frame()])
            sage: M.orientation()
            [Coordinate frame (U, (∂/∂x,∂/∂y)),
             Coordinate frame (V, (∂/∂u,∂/∂v))]
        """
        from sage.manifolds.differentiable.vectorframe import VectorFrame
        chart_type = self._structure.chart
        if isinstance(orientation, chart_type):
            orientation = [orientation.frame()]
        elif isinstance(orientation, VectorFrame):
            orientation = [orientation]
        elif isinstance(orientation, (list, tuple)):
            if isinstance(orientation[0], chart_type):
                orientation = [c.frame() for c in orientation]
            else:
                orientation = list(orientation)
        else:
            raise TypeError("orientation must be a chart/frame or a "
                            "list/tuple of charts/frames")
        dom_union = None
        for frame in orientation:
            if not isinstance(frame, VectorFrame):
                raise ValueError("orientation must consist of vector frames")
            dom = frame._domain
            if not dom.is_subset(self):
                raise ValueError("{} must be defined ".format(frame) +
                                 "on a subset of {}".format(self))
            if dom_union is not None:
                dom_union = dom.union(dom_union)
            else:
                dom_union = dom
        if dom_union != self:
            raise ValueError("frame domains must cover {}".format(self))
        self._orientation = orientation

    def orientation(self):
        r"""
        Get the preferred orientation of ``self`` if available.

        An *orientation* on a differentiable manifold is an atlas of charts
        whose transition maps are pairwise orientation preserving, i.e. whose
        Jacobian determinants are pairwise positive.

        A differentiable manifold with an orientation is called *orientable*.

        A differentiable manifold is orientable if and only if the tangent
        bundle is orientable in terms of a vector bundle,
        see :meth:`~sage.manifolds.vector_bundle.TopologicalVectorBundle.orientation`.

        .. NOTE::

            In contrast to topological manifolds,
            see :meth:`~sage.manifolds.manifold.TopologicalManifold.orientation`,
            differentiable manifolds preferably use the notion of
            orientability in terms of the tangent bundle.

        The trivial case corresponds to the manifold being parallelizable,
        i.e. admitting a frame covering the whole manifold. In that case,
        if no preferred orientation has been manually set before, one of those
        frames (usually the default frame) is set to the preferred
        orientation on ``self`` and returned here.

        EXAMPLES:

        In case one frame already covers the manifold, an orientation
        is readily obtained::

            sage: M = Manifold(3, 'M')
            sage: c.<x,y,z> = M.chart()
            sage: M.orientation()
            [Coordinate frame (M, (∂/∂x,∂/∂y,∂/∂z))]

        However, orientations are usually not easy to obtain::

            sage: M = Manifold(2, 'M')
            sage: U = M.open_subset('U'); V = M.open_subset('V')
            sage: M.declare_union(U, V)
            sage: c_xy.<x,y> = U.chart(); c_uv.<u,v> = V.chart()
            sage: M.orientation()
            []

        In that case, the orientation can be set by the user; either in
        terms of charts or in terms of frames::

            sage: M.set_orientation([c_xy, c_uv])
            sage: M.orientation()
            [Coordinate frame (U, (∂/∂x,∂/∂y)),
             Coordinate frame (V, (∂/∂u,∂/∂v))]
            sage: M.set_orientation([c_xy.frame(), c_uv.frame()])
            sage: M.orientation()
            [Coordinate frame (U, (∂/∂x,∂/∂y)),
             Coordinate frame (V, (∂/∂u,∂/∂v))]

        The orientation on submanifolds are inherited from the ambient
        manifold::

            sage: W = U.intersection(V, name='W')
            sage: W.orientation()
            [Vector frame (W, (∂/∂x,∂/∂y))]
        """
        if not self._orientation:
            # try to get an orientation from super domains:
            for sdom in self.open_supersets():
                sorient = sdom._orientation
                if sorient:
                    rst_orient = [f.restrict(self) for f in sorient]
                    # clear multiple domains:
                    rst_orient = list(self._get_min_covering(rst_orient))
                    self._orientation = rst_orient
                    break
            else:
                # Trivial case:
                if self.is_manifestly_parallelizable():
                    # Try the default frame:
                    def_frame = self._def_frame
                    if def_frame is not None:
                        if def_frame._domain is self:
                            self._orientation = [def_frame]
                    # Still no orientation? Choose arbitrary frame:
                    if not self._orientation:
                        for frame in self._covering_frames:
                            dest_map = frame.destination_map()
                            if dest_map.is_identity():
                                self._orientation = [frame]
                                break
        return list(self._orientation)

    def default_frame(self):
        r"""
        Return the default vector frame defined on ``self``.

        By *vector frame*, it is meant a field on the manifold that provides,
        at each point `p`, a vector basis of the tangent space at `p`.

        Unless changed via :meth:`set_default_frame`, the default frame is
        the first one defined on the manifold, usually implicitly as the
        coordinate basis associated with the first chart defined on the
        manifold.

        OUTPUT:

        - a :class:`~sage.manifolds.differentiable.vectorframe.VectorFrame`
          representing the default vector frame

        EXAMPLES:

        The default vector frame is often the coordinate frame associated
        with the first chart defined on the manifold::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: M.default_frame()
            Coordinate frame (M, (∂/∂x,∂/∂y))
        """
        return self._def_frame

    def set_default_frame(self, frame):
        r"""
        Changing the default vector frame on ``self``.

        INPUT:

        - ``frame`` --
          :class:`~sage.manifolds.differentiable.vectorframe.VectorFrame`
          a vector frame defined on some subset of ``self``

        EXAMPLES:

        Changing the default frame on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: e = M.vector_frame('e')
            sage: M.default_frame()
            Coordinate frame (M, (∂/∂x,∂/∂y))
            sage: M.set_default_frame(e)
            sage: M.default_frame()
            Vector frame (M, (e_0,e_1))
        """
        from sage.manifolds.differentiable.vectorframe import VectorFrame
        if not isinstance(frame, VectorFrame):
            raise TypeError("{} is not a vector frame".format(frame))
        if not frame._domain.is_subset(self):
            raise ValueError("the frame must be defined on the {}".format(self))
        self._def_frame = frame
        frame._fmodule.set_default_basis(frame)

    def change_of_frame(self, frame1, frame2):
        r"""
        Return a change of vector frames defined on ``self``.

        INPUT:

        - ``frame1`` -- vector frame 1
        - ``frame2`` -- vector frame 2

        OUTPUT:

        - a
          :class:`~sage.manifolds.differentiable.automorphismfield.AutomorphismField`
          representing, at each point, the vector space automorphism `P`
          that relates frame 1, `(e_i)` say, to frame 2, `(n_i)` say,
          according to `n_i = P(e_i)`

        EXAMPLES:

        Change of vector frames induced by a change of coordinates::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: c_uv.<u,v> = M.chart()
            sage: c_xy.transition_map(c_uv, (x+y, x-y))
            Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))
            sage: M.change_of_frame(c_xy.frame(), c_uv.frame())
            Field of tangent-space automorphisms on the 2-dimensional
             differentiable manifold M
            sage: M.change_of_frame(c_xy.frame(), c_uv.frame())[:]
            [ 1/2  1/2]
            [ 1/2 -1/2]
            sage: M.change_of_frame(c_uv.frame(), c_xy.frame())
            Field of tangent-space automorphisms on the 2-dimensional
             differentiable manifold M
            sage: M.change_of_frame(c_uv.frame(), c_xy.frame())[:]
            [ 1  1]
            [ 1 -1]
            sage: M.change_of_frame(c_uv.frame(), c_xy.frame()) == \
            ....:       M.change_of_frame(c_xy.frame(), c_uv.frame()).inverse()
            True

        In the present example, the manifold `M` is parallelizable, so
        that the module `X(M)` of vector fields on `M` is free. A change
        of frame on `M` is then identical to a change of basis in `X(M)`::

            sage: XM = M.vector_field_module() ; XM
            Free module X(M) of vector fields on the 2-dimensional
             differentiable manifold M
            sage: XM.print_bases()
            Bases defined on the Free module X(M) of vector fields on the
             2-dimensional differentiable manifold M:
             - (M, (∂/∂x,∂/∂y)) (default basis)
             - (M, (∂/∂u,∂/∂v))
            sage: XM.change_of_basis(c_xy.frame(), c_uv.frame())
            Field of tangent-space automorphisms on the 2-dimensional
             differentiable manifold M
            sage: M.change_of_frame(c_xy.frame(), c_uv.frame()) is \
            ....:  XM.change_of_basis(c_xy.frame(), c_uv.frame())
            True
        """
        if (frame1, frame2) not in self._frame_changes:
            raise ValueError("the change of frame from {} to {}".format(frame1, frame2) +
                             " has not been defined on the {}".format(self))
        return self._frame_changes[(frame1, frame2)]

    def set_change_of_frame(self, frame1, frame2, change_of_frame,
                         compute_inverse=True):
        r"""
        Relate two vector frames by an automorphism.

        This updates the internal dictionary ``self._frame_changes``.

        INPUT:

        - ``frame1`` -- frame 1, denoted `(e_i)` below
        - ``frame2`` -- frame 2, denoted `(f_i)` below
        - ``change_of_frame`` -- instance of class
          :class:`~sage.manifolds.differentiable.automorphismfield.AutomorphismFieldParal`
          describing the automorphism `P` that relates the basis `(e_i)` to
          the basis `(f_i)` according to `f_i = P(e_i)`
        - ``compute_inverse`` -- boolean (default: ``True``); if set to True, the inverse
          automorphism is computed and the change from basis `(f_i)` to `(e_i)`
          is set to it in the internal dictionary ``self._frame_changes``

        EXAMPLES:

        Connecting two vector frames on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: e = M.vector_frame('e')
            sage: f = M.vector_frame('f')
            sage: a = M.automorphism_field()
            sage: a[e,:] = [[1,2],[0,3]]
            sage: M.set_change_of_frame(e, f, a)
            sage: f[0].display(e)
            f_0 = e_0
            sage: f[1].display(e)
            f_1 = 2 e_0 + 3 e_1
            sage: e[0].display(f)
            e_0 = f_0
            sage: e[1].display(f)
            e_1 = -2/3 f_0 + 1/3 f_1
            sage: M.change_of_frame(e,f)[e,:]
            [1 2]
            [0 3]
        """
        from sage.manifolds.differentiable.automorphismfield import (
            AutomorphismFieldParal,
        )
        fmodule = frame1._fmodule
        if frame2._fmodule != fmodule:
            raise ValueError("the two frames are not defined on the same " +
                             "vector field module")
        if not isinstance(change_of_frame, AutomorphismFieldParal):
            raise TypeError("the argument change_of_frame must be some " +
                            "instance of AutomorphismFieldParal")
        fmodule.set_change_of_basis(frame1, frame2, change_of_frame,
                                    compute_inverse=compute_inverse)
        for sdom in self.open_supersets():
            sdom._frame_changes[(frame1, frame2)] = change_of_frame
        if compute_inverse:
            if (frame2, frame1) not in self._frame_changes:
                for sdom in self.open_supersets():
                    sdom._frame_changes[(frame2, frame1)] = change_of_frame.inverse()

    def vector_frame(self, *args, **kwargs) -> VectorFrame:
        r"""
        Define a vector frame on ``self``.

        A *vector frame* is a field on the manifold that provides, at each
        point `p` of the manifold, a vector basis of the tangent space at `p`
        (or at `\Phi(p)` when ``dest_map`` is not ``None``, see below).

        The vector frame can be defined from a set of `n` linearly independent
        vector fields, `n` being the dimension of ``self``.

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.vectorframe.VectorFrame`
            for complete documentation.

        INPUT:

        - ``symbol`` -- either a string, to be used as a
          common base for the symbols of the vector fields constituting the
          vector frame, or a list/tuple of strings, representing the individual
          symbols of the vector fields; can be omitted only if ``from_frame``
          is not ``None`` (see below)
        - ``vector_fields`` -- tuple or list of `n` linearly independent vector
          fields on the manifold ``self`` (`n` being the dimension of ``self``)
          defining the vector frame; can be omitted if the vector frame is
          created from scratch or if ``from_frame`` is not ``None``
        - ``latex_symbol`` -- (default: ``None``) either a string, to be used
          as a common base for the LaTeX symbols of the vector fields
          constituting the vector frame, or a list/tuple of strings,
          representing the individual LaTeX symbols of the vector fields;
          if ``None``, ``symbol`` is used in place of ``latex_symbol``
        - ``dest_map`` -- (default: ``None``)
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`;
          destination map `\Phi:\ U \rightarrow M`, where `U` is ``self`` and
          `M` is a differentiable manifold; for each `p\in U`, the vector
          frame evaluated at `p` is a basis of the tangent space
          `T_{\Phi(p)}M`; if ``dest_map`` is ``None``, the identity map is
          assumed (case of a vector frame *on* `U`)
        - ``from_frame`` -- (default: ``None``) vector frame `\tilde{e}`
          on the codomain `M` of the destination map `\Phi`; the returned
          frame `e` is then such that for all `p \in U`,
          we have `e(p) = \tilde{e}(\Phi(p))`
        - ``indices`` -- (default: ``None``; used only if ``symbol`` is a
          single string) tuple of strings representing the indices labelling
          the vector fields of the frame; if ``None``, the indices will be
          generated as integers within the range declared on ``self``
        - ``latex_indices`` -- (default: ``None``) tuple of strings
          representing the indices for the LaTeX symbols of the vector fields;
          if ``None``, ``indices`` is used instead
        - ``symbol_dual`` -- (default: ``None``) same as ``symbol`` but for the
          dual coframe; if ``None``, ``symbol`` must be a string and is used
          for the common base of the symbols of the elements of the dual
          coframe
        - ``latex_symbol_dual`` -- (default: ``None``) same as ``latex_symbol``
          but for the dual coframe

        OUTPUT:

        - a :class:`~sage.manifolds.differentiable.vectorframe.VectorFrame`
          representing the defined vector frame

        EXAMPLES:

        Defining a vector frame from two linearly independent vector
        fields on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: e0 = M.vector_field(1+x^2, 1+y^2)
            sage: e1 = M.vector_field(2, -x*y)
            sage: e = M.vector_frame('e', (e0, e1)); e
            Vector frame (M, (e_0,e_1))
            sage: e[0].display()
            e_0 = (x^2 + 1) ∂/∂x + (y^2 + 1) ∂/∂y
            sage: e[1].display()
            e_1 = 2 ∂/∂x - x*y ∂/∂y
            sage: (e[0], e[1]) == (e0, e1)
            True

        If the vector fields are not linearly independent, an error is
        raised::

            sage: z = M.vector_frame('z', (e0, -e0))
            Traceback (most recent call last):
            ...
            ValueError: the provided vector fields are not linearly
             independent

        Another example, involving a pair vector fields along a curve::

            sage: R.<t> = manifolds.RealLine()
            sage: c = M.curve([sin(t), sin(2*t)/2], (t, 0, 2*pi), name='c')
            sage: I = c.domain(); I
            Real interval (0, 2*pi)
            sage: v = c.tangent_vector_field()
            sage: v.display()
            c' = cos(t) ∂/∂x + (2*cos(t)^2 - 1) ∂/∂y
            sage: w = I.vector_field(1-2*cos(t)^2, cos(t), dest_map=c)
            sage: u = I.vector_frame('u', (v, w))
            sage: u[0].display()
            u_0 = cos(t) ∂/∂x + (2*cos(t)^2 - 1) ∂/∂y
            sage: u[1].display()
            u_1 = (-2*cos(t)^2 + 1) ∂/∂x + cos(t) ∂/∂y
            sage: (u[0], u[1]) == (v, w)
            True

        It is also possible to create a vector frame from scratch, without
        connecting it to previously defined vector frames or vector fields
        (this can still be performed later via the method
        :meth:`~sage.manifolds.differentiable.manifold.DifferentiableManifold.set_change_of_frame`)::

            sage: f = M.vector_frame('f'); f
            Vector frame (M, (f_0,f_1))
            sage: f[0]
            Vector field f_0 on the 2-dimensional differentiable manifold M

        Thanks to the keywords ``dest_map`` and ``from_frame``, one can also
        define a vector frame from one preexisting on another manifold, via a
        differentiable map (here provided by the curve ``c``)::

            sage: fc = I.vector_frame(dest_map=c, from_frame=f); fc
            Vector frame ((0, 2*pi), (f_0,f_1)) with values on the
             2-dimensional differentiable manifold M
            sage: fc[0]
            Vector field f_0 along the Real interval (0, 2*pi) with values on
             the 2-dimensional differentiable manifold M

        Note that the symbol for ``fc``, namely `f`, is inherited from ``f``,
        the original vector frame.

        .. SEEALSO::

            For more options, in particular for the choice of symbols and
            indices, see
            :class:`~sage.manifolds.differentiable.vectorframe.VectorFrame`.
        """
        from sage.manifolds.differentiable.vectorframe import VectorFrame
        # Input processing
        symbol = None
        vector_fields = None
        n_args = len(args)
        if n_args >= 1:
            symbol = args[0]
        if n_args == 2:
            vector_fields = args[1]
        elif n_args > 2:
            raise TypeError("vector_frame() takes at most two positional "
                            "arguments")
        latex_symbol = kwargs.pop('latex_symbol', None)
        dest_map = kwargs.pop('dest_map', None)
        from_frame = kwargs.pop('from_frame', None)
        indices = kwargs.pop('indices', None)
        latex_indices = kwargs.pop('latex_indices', None)
        symbol_dual = kwargs.pop('symbol_dual', None)
        latex_symbol_dual = kwargs.pop('latex_symbol_dual', None)

        if vector_fields:
            dest_map0 = vector_fields[0].parent().destination_map()
            if dest_map and dest_map is not dest_map0:
                raise ValueError("incompatible values of destination maps")
            dest_map = dest_map0
        resu = VectorFrame(self.vector_field_module(dest_map=dest_map,
                                                    force_free=True),
                           symbol=symbol, latex_symbol=latex_symbol,
                           from_frame=from_frame, indices=indices,
                           latex_indices=latex_indices, symbol_dual=symbol_dual,
                           latex_symbol_dual=latex_symbol_dual)
        if vector_fields:
            linked = False
            try:
                resu._init_from_family(vector_fields)
            except ArithmeticError as err:
                linked = str(err) in ["non-invertible matrix",
                                      "input matrix must be nonsingular"]
            if linked:
                raise ValueError("the provided vector fields are not "
                                 "linearly independent")
            # Adding the newly generated changes of frame to the
            # dictionary _frame_changes of self and its supersets:
            for frame_pair, chge in resu._fmodule._basis_changes.items():
                if resu in frame_pair:
                    for sdom in self.open_supersets():
                        sdom._frame_changes[frame_pair] = chge
        return resu

    def _set_covering_frame(self, frame):
        r"""
        Declare a frame covering ``self``.

        This helper method is invoked by the frame constructor.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: M._covering_frames
            []
            sage: e = M.vector_frame('e')
            sage: M._covering_frames
            [Vector frame (M, (e_0,e_1))]
            sage: M._covering_frames = []
            sage: M._set_covering_frame(e)
            sage: M._covering_frames
            [Vector frame (M, (e_0,e_1))]
        """
        self._covering_frames.append(frame)
        self._parallelizable_parts = set([self])
        # if self contained smaller parallelizable parts, they are forgotten
        for sd in self.open_supersets():
            if not sd.is_manifestly_parallelizable():
                sd._parallelizable_parts.add(self)

    def frames(self):
        r"""
        Return the list of vector frames defined on open subsets of ``self``.

        OUTPUT: list of vector frames defined on open subsets of ``self``

        EXAMPLES:

        Vector frames on subsets of `\RR^2`::

            sage: M = Manifold(2, 'R^2')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: M.frames()
            [Coordinate frame (R^2, (∂/∂x,∂/∂y))]
            sage: e = M.vector_frame('e')
            sage: M.frames()
            [Coordinate frame (R^2, (∂/∂x,∂/∂y)),
             Vector frame (R^2, (e_0,e_1))]
            sage: U = M.open_subset('U', coord_def={c_cart: x^2+y^2<1}) # unit disk
            sage: U.frames()
            [Coordinate frame (U, (∂/∂x,∂/∂y))]
            sage: M.frames()
            [Coordinate frame (R^2, (∂/∂x,∂/∂y)),
             Vector frame (R^2, (e_0,e_1)),
             Coordinate frame (U, (∂/∂x,∂/∂y))]
        """
        return list(self._frames)

    def coframes(self):
        r"""
        Return the list of coframes defined on open subsets of ``self``.

        OUTPUT: list of coframes defined on open subsets of ``self``

        EXAMPLES:

        Coframes on subsets of `\RR^2`::

            sage: M = Manifold(2, 'R^2')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: M.coframes()
            [Coordinate coframe (R^2, (dx,dy))]
            sage: e = M.vector_frame('e')
            sage: M.coframes()
            [Coordinate coframe (R^2, (dx,dy)), Coframe (R^2, (e^0,e^1))]
            sage: U = M.open_subset('U', coord_def={c_cart: x^2+y^2<1}) # unit disk
            sage: U.coframes()
            [Coordinate coframe (U, (dx,dy))]
            sage: e.restrict(U)
            Vector frame (U, (e_0,e_1))
            sage: U.coframes()
            [Coordinate coframe (U, (dx,dy)), Coframe (U, (e^0,e^1))]
            sage: M.coframes()
            [Coordinate coframe (R^2, (dx,dy)),
             Coframe (R^2, (e^0,e^1)),
             Coordinate coframe (U, (dx,dy)),
             Coframe (U, (e^0,e^1))]
        """
        return list(self._coframes)

    def changes_of_frame(self):
        r"""
        Return all the changes of vector frames defined on ``self``.

        OUTPUT:

        - dictionary of fields of tangent-space automorphisms representing
          the changes of frames, the keys being the pair of frames

        EXAMPLES:

        Let us consider a first vector frame on a 2-dimensional
        differentiable manifold::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: e = X.frame(); e
            Coordinate frame (M, (∂/∂x,∂/∂y))

        At this stage, the dictionary of changes of frame is empty::

            sage: M.changes_of_frame()
            {}

        We introduce a second frame on the manifold, relating it to
        frame ``e`` by a field of tangent space automorphisms::

            sage: a = M.automorphism_field(name='a')
            sage: a[:] = [[-y, x], [1, 2]]
            sage: f = e.new_frame(a, 'f'); f
            Vector frame (M, (f_0,f_1))

        Then we have::

            sage: M.changes_of_frame()  # random (dictionary output)
            {(Coordinate frame (M, (∂/∂x,∂/∂y)),
              Vector frame (M, (f_0,f_1))): Field of tangent-space
               automorphisms on the 2-dimensional differentiable manifold M,
             (Vector frame (M, (f_0,f_1)),
              Coordinate frame (M, (∂/∂x,∂/∂y))): Field of tangent-space
               automorphisms on the 2-dimensional differentiable manifold M}

        Some checks::

            sage: M.changes_of_frame()[(e,f)] == a
            True
            sage: M.changes_of_frame()[(f,e)] == a^(-1)
            True
        """
        return self._frame_changes.copy()

    def is_manifestly_parallelizable(self):
        r"""
        Return ``True`` if ``self`` is known to be a parallelizable
        and ``False`` otherwise.

        If ``False`` is returned, either the manifold is not parallelizable
        or no vector frame has been defined on it yet.

        EXAMPLES:

        A just created manifold is a priori not manifestly parallelizable::

            sage: M = Manifold(2, 'M')
            sage: M.is_manifestly_parallelizable()
            False

        Defining a vector frame on it makes it parallelizable::

            sage: e = M.vector_frame('e')
            sage: M.is_manifestly_parallelizable()
            True

        Defining a coordinate chart on the whole manifold also makes it
        parallelizable::

            sage: N = Manifold(4, 'N')
            sage: X.<t,x,y,z> = N.chart()
            sage: N.is_manifestly_parallelizable()
            True
        """
        return bool(self._covering_frames)

    def tangent_space(self, point, base_ring=None):
        r"""
        Tangent space to ``self`` at a given point.

        INPUT:

        - ``point`` -- :class:`~sage.manifolds.point.ManifoldPoint`;
          point `p` on the manifold

        - ``base_ring`` -- (default: the symbolic ring) the base ring

        OUTPUT:

        - :class:`~sage.manifolds.differentiable.tangent_space.TangentSpace`
          representing the tangent vector space `T_{p} M`, where `M` is the
          current manifold

        EXAMPLES:

        A tangent space to a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((2, -3), name='p')
            sage: Tp = M.tangent_space(p); Tp
            Tangent space at Point p on the 2-dimensional differentiable
             manifold M
            sage: Tp.category()
            Category of finite dimensional vector spaces over Symbolic Ring
            sage: dim(Tp)
            2

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.tangent_space.TangentSpace`
            for more examples.
        """
        from sage.manifolds.differentiable.tangent_space import TangentSpace
        from sage.manifolds.point import ManifoldPoint
        if not isinstance(point, ManifoldPoint):
            raise TypeError("{} is not a manifold point".format(point))
        if point not in self:
            raise ValueError("{} is not a point on the {}".format(point, self))
        return TangentSpace(point, base_ring=base_ring)

    def curve(self, coord_expression, param, chart=None,
              name=None, latex_name=None):
        r"""
        Define a differentiable curve in the manifold.

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.curve.DifferentiableCurve`
            for details.

        INPUT:

        - ``coord_expression`` -- either

          - (i) a dictionary whose keys are charts on the manifold and values
            the coordinate expressions (as lists or tuples) of the curve in
            the given chart
          - (ii) a single coordinate expression in a given chart on the
            manifold, the latter being provided by the argument ``chart``

          in both cases, if the dimension of the manifold is 1, a single
          coordinate expression can be passed instead of a tuple with
          a single element
        - ``param`` -- tuple of the type ``(t, t_min, t_max)``, where

          * ``t`` is the curve parameter used in ``coord_expression``;
          * ``t_min`` is its minimal value;
          * ``t_max`` its maximal value;

          if ``t_min=-Infinity`` and ``t_max=+Infinity``, they can be
          omitted and ``t`` can be passed for ``param`` instead of the
          tuple ``(t, t_min, t_max)``
        - ``chart`` -- (default: ``None``) chart on the manifold used for
          case (ii) above; if ``None`` the default chart of the manifold is
          assumed
        - ``name`` -- (default: ``None``) string; symbol given to the curve
        - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote
          the curve; if none is provided, ``name`` will be used

        OUTPUT: :class:`~sage.manifolds.differentiable.curve.DifferentiableCurve`

        EXAMPLES:

        The lemniscate of Gerono in the 2-dimensional Euclidean plane::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: R.<t> = manifolds.RealLine()
            sage: c = M.curve([sin(t), sin(2*t)/2], (t, 0, 2*pi), name='c') ; c
            Curve c in the 2-dimensional differentiable manifold M

        The same definition with the coordinate expression passed as a
        dictionary::

            sage: c = M.curve({X: [sin(t), sin(2*t)/2]}, (t, 0, 2*pi), name='c') ; c
            Curve c in the 2-dimensional differentiable manifold M

        An example of definition with ``t_min`` and ``t_max`` omitted: a helix
        in `\RR^3`::

            sage: R3 = Manifold(3, 'R^3')
            sage: X.<x,y,z> = R3.chart()
            sage: c = R3.curve([cos(t), sin(t), t], t, name='c') ; c
            Curve c in the 3-dimensional differentiable manifold R^3
            sage: c.domain() # check that t is unbounded
            Real number line ℝ

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.curve.DifferentiableCurve`
            for more examples, including plots.
        """
        from sage.manifolds.differentiable.examples.real_line import RealLine
        if not isinstance(param, (tuple, list)):
            param = (param, minus_infinity, infinity)
        elif len(param) != 3:
            raise ValueError("the argument 'param' must be of the form " +
                             "(t, t_min, t_max)")
        t = param[0]
        t_min = param[1]
        t_max = param[2]
        real_field = RealLine(names=(repr(t),))
        interval = real_field.open_interval(t_min, t_max)
        curve_set = Hom(interval, self)
        if not isinstance(coord_expression, dict):
            # Turn coord_expression into a dictionary:
            if chart is None:
                chart = self._def_chart
            elif chart not in self._atlas:
                raise ValueError("the {} has not been ".format(chart) +
                                 "defined on the {}".format(self))
            if isinstance(coord_expression, (tuple, list)):
                coord_expression = {chart: coord_expression}
            else:
                # case self.dim()=1
                coord_expression = {chart: (coord_expression,)}
        return curve_set(coord_expression, name=name, latex_name=latex_name)

    def integrated_curve(self, equations_rhs, velocities, curve_param,
                         initial_tangent_vector, chart=None, name=None,
                         latex_name=None, verbose=False, across_charts=False):
        r"""
        Construct a curve defined by a system of second order
        differential equations in the coordinate functions.

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.integrated_curve.IntegratedCurve`
            for details.

        INPUT:

        - ``equations_rhs`` -- list of the right-hand sides of the
          equations on the velocities only
        - ``velocities`` -- list of the symbolic expressions used in
          ``equations_rhs`` to denote the velocities
        - ``curve_param`` -- tuple of the type ``(t, t_min, t_max)``,
          where

          * ``t`` is the symbolic variable used in ``equations_rhs`` to
            denote the parameter of the curve;
          * ``t_min`` is its minimal (finite) value;
          * ``t_max`` its maximal (finite) value.

        - ``initial_tangent_vector`` --
          :class:`~sage.manifolds.differentiable.tangent_vector.TangentVector`;
          initial tangent vector of the curve
        - ``chart`` -- (default: ``None``) chart on the manifold in
          which the equations are given; if ``None`` the default chart
          of the manifold is assumed
        - ``name`` -- (default: ``None``) string; symbol given to the curve
        - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote
          the curve; if none is provided, ``name`` will be used

        OUTPUT: :class:`~sage.manifolds.differentiable.integrated_curve.IntegratedCurve`

        EXAMPLES:

        Trajectory of a particle of unit mass and unit charge in a
        unit, uniform, stationary magnetic field::

            sage: M = Manifold(3, 'M')
            sage: X.<x1,x2,x3> = M.chart()
            sage: t = var('t')
            sage: D = X.symbolic_velocities()
            sage: eqns = [D[1], -D[0], SR(0)]
            sage: p = M.point((0,0,0), name='p')
            sage: Tp = M.tangent_space(p)
            sage: v = Tp((1,0,1))
            sage: c = M.integrated_curve(eqns, D, (t,0,6), v, name='c'); c
            Integrated curve c in the 3-dimensional differentiable
             manifold M
            sage: sys = c.system(verbose=True)
            Curve c in the 3-dimensional differentiable manifold M
             integrated over the Real interval (0, 6) as a solution to
             the following system, written with respect to
             Chart (M, (x1, x2, x3)):
            <BLANKLINE>
            Initial point: Point p on the 3-dimensional differentiable
             manifold M with coordinates [0, 0, 0] with respect to
             Chart (M, (x1, x2, x3))
            Initial tangent vector: Tangent vector at Point p on the
             3-dimensional differentiable manifold M with
             components [1, 0, 1] with respect to Chart (M, (x1, x2, x3))
            <BLANKLINE>
            d(x1)/dt = Dx1
            d(x2)/dt = Dx2
            d(x3)/dt = Dx3
            d(Dx1)/dt = Dx2
            d(Dx2)/dt = -Dx1
            d(Dx3)/dt = 0
            <BLANKLINE>
            sage: sol = c.solve()
            sage: interp = c.interpolate()
            sage: p = c(1.3, verbose=True)
            Evaluating point coordinates from the interpolation
             associated with the key 'cubic spline-interp-odeint'
             by default...
            sage: p
            Point on the 3-dimensional differentiable manifold M
            sage: p.coordinates()     # abs tol 1e-12
            (0.9635581599167499, -0.7325011788437327, 1.3)
            sage: tgt_vec = c.tangent_vector_eval_at(3.7, verbose=True)
            Evaluating tangent vector components from the interpolation
             associated with the key 'cubic spline-interp-odeint'
             by default...
            sage: tgt_vec[:]    # abs tol 1e-12
            [-0.8481007454066425, 0.5298350137284363, 1.0]
        """

        from sage.manifolds.differentiable.examples.real_line import RealLine
        from sage.manifolds.differentiable.manifold_homset import IntegratedCurveSet

        if len(curve_param) != 3:
            raise ValueError("the argument 'curve_param' must be of the form " +
                             "(t, t_min, t_max)")
        t = curve_param[0]
        t_min = curve_param[1]
        t_max = curve_param[2]
        real_field = RealLine(names=(repr(t),))
        interval = real_field.open_interval(t_min, t_max)
        integrated_curve_set = IntegratedCurveSet(interval, self) # not
        # possible to use Hom(interval, self)
        return integrated_curve_set(equations_rhs, velocities, t,
                                    initial_tangent_vector, chart=chart,
                                    name=name, latex_name=latex_name,
                                    verbose=verbose, across_charts=across_charts)

    def integrated_autoparallel_curve(self, affine_connection,
                        curve_param, initial_tangent_vector, chart=None,
                        name=None, latex_name=None, verbose=False,
                        across_charts=False):
        r"""
        Construct an autoparallel curve on the manifold with respect to
        a given affine connection.

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.integrated_curve.IntegratedAutoparallelCurve`
            for details.

        INPUT:

        - ``affine_connection`` --
          :class:`~sage.manifolds.differentiable.affine_connection.AffineConnection`;
          affine connection with respect to which the curve is autoparallel
        - ``curve_param`` -- tuple of the type ``(t, t_min, t_max)``,
          where

          * ``t`` is the symbolic variable to be used as the parameter
            of the curve (the equations defining an instance of
            :class:`~sage.manifolds.differentiable.integrated_curve.IntegratedAutoparallelCurve`
            are such that ``t`` will actually be an affine parameter
            of the curve);
          * ``t_min`` is its minimal (finite) value;
          * ``t_max`` its maximal (finite) value.

        - ``initial_tangent_vector`` --
          :class:`~sage.manifolds.differentiable.tangent_vector.TangentVector`;
          initial tangent vector of the curve
        - ``chart`` -- (default: ``None``) chart on the manifold in
          which the equations are given ; if ``None`` the default chart
          of the manifold is assumed
        - ``name`` -- (default: ``None``) string; symbol given to the curve
        - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote
          the curve; if none is provided, ``name`` will be used

        OUTPUT: :class:`~sage.manifolds.differentiable.integrated_curve.IntegratedAutoparallelCurve`

        EXAMPLES:

        Autoparallel curves associated with the Mercator projection of
        the 2-sphere `\mathbb{S}^{2}`::

            sage: S2 = Manifold(2, 'S^2', start_index=1)
            sage: polar.<th,ph> = S2.chart('th ph')
            sage: epolar = polar.frame()
            sage: ch_basis = S2.automorphism_field()
            sage: ch_basis[1,1], ch_basis[2,2] = 1, 1/sin(th)
            sage: epolar_ON=S2.default_frame().new_frame(ch_basis,'epolar_ON')

        Set the affine connection associated with Mercator projection;
        it is metric compatible but it has non-vanishing torsion::

            sage: nab = S2.affine_connection('nab')
            sage: nab.set_coef(epolar_ON)[:]
            [[[0, 0], [0, 0]], [[0, 0], [0, 0]]]
            sage: g = S2.metric('g')
            sage: g[1,1], g[2,2] = 1, (sin(th))^2
            sage: nab(g)[:]
            [[[0, 0], [0, 0]], [[0, 0], [0, 0]]]
            sage: nab.torsion()[:]
            [[[0, 0], [0, 0]], [[0, cos(th)/sin(th)], [-cos(th)/sin(th), 0]]]

        Declare an integrated autoparallel curve with respect to this
        connection::

            sage: p = S2.point((pi/4, 0), name='p')
            sage: Tp = S2.tangent_space(p)
            sage: v = Tp((1,1), basis=epolar_ON.at(p))
            sage: t = var('t')
            sage: c = S2.integrated_autoparallel_curve(nab, (t, 0, 2.3),
            ....:                              v, chart=polar, name='c')
            sage: sys = c.system(verbose=True)
            Autoparallel curve c in the 2-dimensional differentiable
             manifold S^2 equipped with Affine connection nab on the
             2-dimensional differentiable manifold S^2, and integrated
             over the Real interval (0, 2.30000000000000) as a solution to the
             following equations, written with respect to
             Chart (S^2, (th, ph)):
            <BLANKLINE>
            Initial point: Point p on the 2-dimensional differentiable
             manifold S^2 with coordinates [1/4*pi, 0] with respect to
             Chart (S^2, (th, ph))
            Initial tangent vector: Tangent vector at Point p on the
             2-dimensional differentiable manifold S^2 with
             components [1, sqrt(2)] with respect to
             Chart (S^2, (th, ph))
            <BLANKLINE>
            d(th)/dt = Dth
            d(ph)/dt = Dph
            d(Dth)/dt = 0
            d(Dph)/dt = -Dph*Dth*cos(th)/sin(th)
            <BLANKLINE>
            sage: sol = c.solve()
            sage: interp = c.interpolate()
            sage: p = c(1.3, verbose=True)
            Evaluating point coordinates from the interpolation
             associated with the key 'cubic spline-interp-odeint'
             by default...
            sage: p
            Point on the 2-dimensional differentiable manifold S^2
            sage: polar(p)     # abs tol 1e-12
            (2.0853981633974477, 1.4203177070475606)
            sage: tgt_vec = c.tangent_vector_eval_at(1.3, verbose=True)
            Evaluating tangent vector components from the interpolation
             associated with the key 'cubic spline-interp-odeint'
             by default...
            sage: tgt_vec[:]    # abs tol 1e-12
            [1.000000000000011, 1.148779968412235]
        """

        from sage.manifolds.differentiable.examples.real_line import RealLine
        from sage.manifolds.differentiable.manifold_homset import (
            IntegratedAutoparallelCurveSet,
        )

        if len(curve_param) != 3:
            raise ValueError("the argument 'curve_param' must be " +
                             "of the form (t, t_min, t_max)")
        t = curve_param[0]
        t_min = curve_param[1]
        t_max = curve_param[2]
        real_field = RealLine(names=(repr(t),))
        interval = real_field.open_interval(t_min, t_max)
        autoparallel_curve_set = IntegratedAutoparallelCurveSet(interval,
                                                                self)
        # not possible to use Hom(interval, self)
        return autoparallel_curve_set(affine_connection, t,
                                      initial_tangent_vector,
                                      chart=chart, name=name,
                                      latex_name=latex_name,
                                      verbose=verbose,
                                      across_charts=across_charts)

    def integrated_geodesic(self, metric, curve_param,
                            initial_tangent_vector, chart=None,
                            name=None, latex_name=None, verbose=False,
                            across_charts=False):
        r"""
        Construct a geodesic on the manifold with respect to a given metric.

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.integrated_curve.IntegratedGeodesic`
            for details.

        INPUT:

        - ``metric`` --
          :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`
          metric with respect to which the curve is a geodesic
        - ``curve_param`` -- tuple of the type ``(t, t_min, t_max)``,
          where

          * ``t`` is the symbolic variable to be used as the parameter
            of the curve (the equations defining an instance of
            :class:`~sage.manifolds.differentiable.integrated_curve.IntegratedGeodesic`
            are such that ``t`` will actually be an affine parameter
            of the curve);
          * ``t_min`` is its minimal (finite) value;
          * ``t_max`` its maximal (finite) value.

        - ``initial_tangent_vector`` --
          :class:`~sage.manifolds.differentiable.tangent_vector.TangentVector`;
          initial tangent vector of the curve
        - ``chart`` -- (default: ``None``) chart on the manifold in
          which the equations are given; if ``None`` the default chart
          of the manifold is assumed
        - ``name`` -- (default: ``None``) string; symbol given to the curve
        - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote
          the curve; if none is provided, ``name`` will be used

        OUTPUT: :class:`~sage.manifolds.differentiable.integrated_curve.IntegratedGeodesic`

        EXAMPLES:

        Geodesics of the unit 2-sphere `\mathbb{S}^{2}`::

            sage: S2 = Manifold(2, 'S^2', start_index=1)
            sage: polar.<th,ph> = S2.chart('th ph')
            sage: epolar = polar.frame()

        Set the standard metric tensor `g` on `\mathbb{S}^{2}`::

            sage: g = S2.metric('g')
            sage: g[1,1], g[2,2] = 1, (sin(th))^2

        Declare an integrated geodesic with respect to this metric::

            sage: p = S2.point((pi/4, 0), name='p')
            sage: Tp = S2.tangent_space(p)
            sage: v = Tp((1, 1), basis=epolar.at(p))
            sage: t = var('t')
            sage: c = S2.integrated_geodesic(g, (t, 0, 6), v,
            ....:                                 chart=polar, name='c')
            sage: sys = c.system(verbose=True)
            Geodesic c in the 2-dimensional differentiable manifold S^2
             equipped with Riemannian metric g on the 2-dimensional
             differentiable manifold S^2, and integrated over the Real
             interval (0, 6) as a solution to the following geodesic
             equations, written with respect to Chart (S^2, (th, ph)):
            <BLANKLINE>
            Initial point: Point p on the 2-dimensional differentiable
            manifold S^2 with coordinates [1/4*pi, 0] with respect to
            Chart (S^2, (th, ph))
            Initial tangent vector: Tangent vector at Point p on the
            2-dimensional differentiable manifold S^2 with
            components [1, 1] with respect to Chart (S^2, (th, ph))
            <BLANKLINE>
            d(th)/dt = Dth
            d(ph)/dt = Dph
            d(Dth)/dt = Dph^2*cos(th)*sin(th)
            d(Dph)/dt = -2*Dph*Dth*cos(th)/sin(th)
            <BLANKLINE>
            sage: sol = c.solve()
            sage: interp = c.interpolate()
            sage: p = c(1.3, verbose=True)
            Evaluating point coordinates from the interpolation
             associated with the key 'cubic spline-interp-odeint'
             by default...
            sage: p
            Point on the 2-dimensional differentiable manifold S^2
            sage: p.coordinates()     # abs tol 1e-12
            (2.2047435672397526, 0.7986602654406825)
            sage: tgt_vec = c.tangent_vector_eval_at(3.7, verbose=True)
            Evaluating tangent vector components from the interpolation
             associated with the key 'cubic spline-interp-odeint'
             by default...
            sage: tgt_vec[:]    # abs tol 1e-12
            [-1.0907409234671228, 0.6205670379855032]
        """
        from sage.manifolds.differentiable.examples.real_line import RealLine
        from sage.manifolds.differentiable.manifold_homset import IntegratedGeodesicSet

        if len(curve_param) != 3:
            raise ValueError("the argument 'curve_param' must be of " +
                             "the form (t, t_min, t_max)")
        t = curve_param[0]
        t_min = curve_param[1]
        t_max = curve_param[2]
        real_field = RealLine(names=(repr(t),))
        interval = real_field.open_interval(t_min, t_max)
        integrated_geodesic_set = IntegratedGeodesicSet(interval, self)
        return integrated_geodesic_set(metric, t, initial_tangent_vector,
                                       chart=chart, name=name,
                                       latex_name=latex_name,
                                       verbose=verbose,
                                       across_charts=across_charts)

    def affine_connection(self, name, latex_name=None):
        r"""
        Define an affine connection on the manifold.

        See :class:`~sage.manifolds.differentiable.affine_connection.AffineConnection`
        for a complete documentation.

        INPUT:

        - ``name`` -- name given to the affine connection
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          affine connection

        OUTPUT:

        - the affine connection, as an instance of
          :class:`~sage.manifolds.differentiable.affine_connection.AffineConnection`

        EXAMPLES:

        Affine connection on an open subset of a 3-dimensional smooth manifold::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: A = M.open_subset('A', latex_name=r'\mathcal{A}')
            sage: nab = A.affine_connection('nabla', r'\nabla') ; nab
            Affine connection nabla on the Open subset A of the 3-dimensional
             differentiable manifold M

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.affine_connection.AffineConnection`
            for more examples.
        """
        from sage.manifolds.differentiable.affine_connection import AffineConnection
        return AffineConnection(self, name, latex_name)

    def metric(self, name: str, signature: Optional[int] = None,
               latex_name: Optional[str] = None,
               dest_map: Optional[DiffMap] = None) -> PseudoRiemannianMetric:
        r"""
        Define a pseudo-Riemannian metric on the manifold.

        A *pseudo-Riemannian metric* is a field of nondegenerate symmetric
        bilinear forms acting in the tangent spaces. See
        :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`
        for a complete documentation.

        INPUT:

        - ``name`` -- name given to the metric
        - ``signature`` -- (default: ``None``) signature `S` of the metric as a
          single integer: `S = n_+ - n_-`, where `n_+` (resp. `n_-`) is the
          number of positive terms (resp. number of negative terms) in any
          diagonal writing of the metric components; if ``signature`` is not
          provided, `S` is set to the manifold's dimension (Riemannian
          signature)
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          metric; if ``None``, it is formed from ``name``
        - ``dest_map`` -- (default: ``None``) instance of
          class :class:`~sage.manifolds.differentiable.diff_map.DiffMap`
          representing the destination map `\Phi:\ U \rightarrow M`, where `U`
          is the current manifold; if ``None``, the identity map is assumed
          (case of a metric tensor field *on* `U`)

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`
          representing the defined pseudo-Riemannian metric.

        EXAMPLES:

        Metric on a 3-dimensional manifold::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = M.chart()
            sage: g = M.metric('g'); g
            Riemannian metric g on the 3-dimensional differentiable manifold M

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`
            for more examples.
        """
        vmodule = self.vector_field_module(dest_map)
        return vmodule.metric(name, signature=signature, latex_name=latex_name)

    def degenerate_metric(self, name, latex_name=None, dest_map=None):
        r"""
        Define a degenerate (or null or lightlike) metric on the manifold.

        A *degenerate metric* is a field of degenerate symmetric
        bilinear forms acting in the tangent spaces.

        See
        :class:`~sage.manifolds.differentiable.metric.DegenerateMetric`
        for a complete documentation.

        INPUT:

        - ``name`` -- name given to the metric
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          metric; if ``None``, it is formed from ``name``
        - ``dest_map`` -- (default: ``None``) instance of
          class :class:`~sage.manifolds.differentiable.diff_map.DiffMap`
          representing the destination map `\Phi:\ U \rightarrow M`, where `U`
          is the current manifold; if ``None``, the identity map is assumed
          (case of a metric tensor field *on* `U`)

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.metric.DegenerateMetric`
          representing the defined degenerate metric.

        EXAMPLES:

        Lightlike cone::

            sage: M = Manifold(3, 'M'); X.<x,y,z> = M.chart()
            sage: g = M.degenerate_metric('g'); g
            degenerate metric g on the 3-dimensional differentiable manifold M
            sage: det(g)
            Scalar field zero on the 3-dimensional differentiable manifold M
            sage: g.parent()
            Free module T^(0,2)(M) of type-(0,2) tensors fields on the
            3-dimensional differentiable manifold M
            sage: g[0,0], g[0,1], g[0,2] = (y^2 + z^2)/(x^2 + y^2 + z^2), \
            ....: - x*y/(x^2 + y^2 + z^2), - x*z/(x^2 + y^2 + z^2)
            sage: g[1,1], g[1,2], g[2,2] = (x^2 + z^2)/(x^2 + y^2 + z^2), \
            ....: - y*z/(x^2 + y^2 + z^2), (x^2 + y^2)/(x^2 + y^2 + z^2)
            sage: g.disp()
            g = (y^2 + z^2)/(x^2 + y^2 + z^2) dx⊗dx - x*y/(x^2 + y^2 + z^2) dx⊗dy
            - x*z/(x^2 + y^2 + z^2) dx⊗dz - x*y/(x^2 + y^2 + z^2) dy⊗dx
            + (x^2 + z^2)/(x^2 + y^2 + z^2) dy⊗dy - y*z/(x^2 + y^2 + z^2) dy⊗dz
            - x*z/(x^2 + y^2 + z^2) dz⊗dx - y*z/(x^2 + y^2 + z^2) dz⊗dy
            + (x^2 + y^2)/(x^2 + y^2 + z^2) dz⊗dz

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.metric.DegenerateMetric`
            for more examples.
        """
        vmodule = self.vector_field_module(dest_map)
        dim = vmodule.ambient_domain().dimension()
        return vmodule.metric(name, signature=(0,dim-1,1), latex_name=latex_name)

    def riemannian_metric(self, name, latex_name=None, dest_map=None):
        r"""
        Define a Riemannian metric on the manifold.

        A *Riemannian metric* is a field of positive definite symmetric
        bilinear forms acting in the tangent spaces.

        See
        :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`
        for a complete documentation.

        INPUT:

        - ``name`` -- name given to the metric
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          metric; if ``None``, it is formed from ``name``
        - ``dest_map`` -- (default: ``None``) instance of
          class :class:`~sage.manifolds.differentiable.diff_map.DiffMap`
          representing the destination map `\Phi:\ U \rightarrow M`, where `U`
          is the current manifold; if ``None``, the identity map is assumed
          (case of a metric tensor field *on* `U`)

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`
          representing the defined Riemannian metric.

        EXAMPLES:

        Metric of the hyperbolic plane `H^2`::

            sage: H2 = Manifold(2, 'H^2', start_index=1)
            sage: X.<x,y> = H2.chart('x y:(0,+oo)')  # Poincaré half-plane coord.
            sage: g = H2.riemannian_metric('g')
            sage: g[1,1], g[2,2] = 1/y^2, 1/y^2
            sage: g
            Riemannian metric g on the 2-dimensional differentiable manifold H^2
            sage: g.display()
            g = y^(-2) dx⊗dx + y^(-2) dy⊗dy
            sage: g.signature()
            2

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`
            for more examples.
        """
        vmodule = self.vector_field_module(dest_map)
        dim = vmodule.ambient_domain().dimension()
        return vmodule.metric(name, signature=dim, latex_name=latex_name)

    def lorentzian_metric(self, name, signature='positive', latex_name=None,
                          dest_map=None):
        r"""
        Define a Lorentzian metric on the manifold.

        A *Lorentzian metric* is a field of nondegenerate symmetric bilinear
        forms acting in the tangent spaces, with signature `(-,+,\cdots,+)` or
        `(+,-,\cdots,-)`.

        See
        :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`
        for a complete documentation.

        INPUT:

        - ``name`` -- name given to the metric
        - ``signature`` -- (default: ``'positive'``) sign of the metric
          signature:

          * if set to 'positive', the signature is n-2, where n is the
            manifold's dimension, i.e. `(-,+,\cdots,+)`
          * if set to 'negative', the signature is -n+2, i.e. `(+,-,\cdots,-)`

        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          metric; if ``None``, it is formed from ``name``
        - ``dest_map`` -- (default: ``None``) instance of
          class :class:`~sage.manifolds.differentiable.diff_map.DiffMap`
          representing the destination map `\Phi:\ U \rightarrow M`, where `U`
          is the current manifold; if ``None``, the identity map is assumed
          (case of a metric tensor field *on* `U`)

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`
          representing the defined Lorentzian metric.

        EXAMPLES:

        Metric of Minkowski spacetime::

            sage: M = Manifold(4, 'M')
            sage: X.<t,x,y,z> = M.chart()
            sage: g = M.lorentzian_metric('g'); g
            Lorentzian metric g on the 4-dimensional differentiable manifold M
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1, 1, 1, 1
            sage: g.display()
            g = -dt⊗dt + dx⊗dx + dy⊗dy + dz⊗dz
            sage: g.signature()
            2

        Choice of a negative signature::

            sage: g = M.lorentzian_metric('g', signature='negative'); g
            Lorentzian metric g on the 4-dimensional differentiable manifold M
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = 1, -1, -1, -1
            sage: g.display()
            g = dt⊗dt - dx⊗dx - dy⊗dy - dz⊗dz
            sage: g.signature()
            -2
        """
        vmodule = self.vector_field_module(dest_map)
        dim = vmodule.ambient_domain().dimension()
        if signature == 'positive':
            signat = dim - 2
        else:
            signat = 2 - dim
        return vmodule.metric(name, signature=signat, latex_name=latex_name)

    def tangent_vector(self, *args, **kwargs):
        r"""
        Define a tangent vector at a given point of ``self``.

        INPUT:

        - ``point`` -- :class:`~sage.manifolds.point.ManifoldPoint`;
          point `p` on ``self``
        - ``comp`` -- components of the vector with respect to the basis
          specified by the argument ``basis``, either as an iterable or as a
          sequence of `n` components, `n` being the dimension of ``self`` (see
          examples below)
        - ``basis`` -- (default: ``None``)
          :class:`~sage.tensor.modules.free_module_basis.FreeModuleBasis`;
          basis of the tangent space at `p` with respect to which the
          components are defined; if ``None``, the default basis of the tangent
          space is used
        - ``name`` -- (default: ``None``) string; symbol given to the vector
        - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote
          the vector; if ``None``, ``name`` will be used

        OUTPUT:

        - :class:`~sage.manifolds.differentiable.tangent_vector.TangentVector`
          representing the tangent vector at point `p`

        EXAMPLES:

        Vector at a point `p` of the Euclidean plane::

            sage: E.<x,y>= EuclideanSpace()
            sage: p = E((1, 2), name='p')
            sage: v = E.tangent_vector(p, -1, 3, name='v'); v
            Vector v at Point p on the Euclidean plane E^2
            sage: v.display()
            v = -e_x + 3 e_y
            sage: v.parent()
            Tangent space at Point p on the Euclidean plane E^2
            sage: v in E.tangent_space(p)
            True

        An alias of ``tangent_vector`` is ``vector``::

            sage: v = E.vector(p, -1, 3, name='v'); v
            Vector v at Point p on the Euclidean plane E^2

        The components can be passed as a tuple or a list::

            sage: v1 = E.vector(p, (-1, 3)); v1
            Vector at Point p on the Euclidean plane E^2
            sage: v1 == v
            True

        or as an object created by the ``vector`` function::

            sage: v2 = E.vector(p, vector([-1, 3])); v2
            Vector at Point p on the Euclidean plane E^2
            sage: v2 == v
            True

        Example of use with the options ``basis`` and ``latex_name``::

            sage: polar_basis = E.polar_frame().at(p)
            sage: polar_basis
            Basis (e_r,e_ph) on the Tangent space at Point p on the Euclidean plane E^2
            sage: v = E.vector(p, 2, -1, basis=polar_basis, name='v',
            ....:              latex_name=r'\vec{v}')
            sage: v
            Vector v at Point p on the Euclidean plane E^2
            sage: v.display(polar_basis)
            v = 2 e_r - e_ph
            sage: v.display()
            v = 4/5*sqrt(5) e_x + 3/5*sqrt(5) e_y
            sage: latex(v)
            \vec{v}

        TESTS::

            sage: E.vector(-1, 3)
            Traceback (most recent call last):
            ...
            TypeError: -1 is not a manifold point
            sage: E.vector([-1, 3])
            Traceback (most recent call last):
            ...
            TypeError: a point and a set of components must be provided
            sage: E.vector(p, 4, 2, 1)
            Traceback (most recent call last):
            ...
            ValueError: 2 components must be provided
        """
        basis = kwargs.pop('basis', None)
        name = kwargs.pop('name', None)
        latex_name = kwargs.pop('latex_name', None)
        if len(args) < 2:
            raise TypeError("a point and a set of components must be provided")
        point = args[0]
        tspace = self.tangent_space(point) # checks on point are performed here
        comp0 = args[1]
        if hasattr(comp0, '__len__') and hasattr(comp0, '__getitem__'):
            # comp0 is a list/vector of components
            comp = comp0
        else:
            # the components are provided as args[1], args[2], ..., args[dim]
            dim = self._dim
            if len(args) != dim + 1:
                raise ValueError(f"{dim} components must be provided")
            comp = args[1:dim + 1]
        return tspace._element_constructor_(comp=comp, basis=basis, name=name,
                                            latex_name=latex_name)

    vector = tangent_vector
