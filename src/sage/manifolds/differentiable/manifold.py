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
    <type 'sage.symbolic.expression.Expression'>

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

    sage: M.list_of_subsets()
    [2-dimensional differentiable manifold S^2,
     Open subset U of the 2-dimensional differentiable manifold S^2,
     Open subset V of the 2-dimensional differentiable manifold S^2,
     Open subset W of the 2-dimensional differentiable manifold S^2]

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
    f: S^2 --> R
    on U: (x, y) |--> arctan(x^2 + y^2)
    on V: (u, v) |--> 1/2*pi - arctan(u^2 + v^2)
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
    Category of commutative algebras over Symbolic Ring

A differentiable manifold has a default vector frame, which, unless otherwise
specified, is the coordinate frame associated with the first defined chart::

    sage: M.default_frame()
    Coordinate frame (U, (d/dx,d/dy))
    sage: latex(M.default_frame())
    \left(U, \left(\frac{\partial}{\partial x },\frac{\partial}{\partial y }\right)\right)
    sage: M.default_frame() is stereoN.frame()
    True

A vector field on the sphere::

    sage: w = M.vector_field('w')
    sage: w[stereoN.frame(), :] = [x, y]
    sage: w.add_comp_by_continuation(stereoS.frame(), W, stereoS)
    sage: w.display() # display in the default frame (stereoN.frame())
    w = x d/dx + y d/dy
    sage: w.display(stereoS.frame())
    w = -u d/du - v d/dv
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
    w(f): S^2 --> R
    on U: (x, y) |--> 2*(x^2 + y^2)/(x^4 + 2*x^2*y^2 + y^4 + 1)
    on V: (u, v) |--> 2*(u^2 + v^2)/(u^4 + 2*u^2*v^2 + v^4 + 1)
    sage: w(f) == f.differential()(w)
    True

The value of the vector field at point `p` is a vector tangent to the sphere::

    sage: w.at(p)
    Tangent vector w at Point p on the 2-dimensional differentiable manifold S^2
    sage: w.at(p).display()
    w = d/dx + 2 d/dy
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
    Module /\^1(S^2) of 1-forms on the 2-dimensional differentiable manifold S^2
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

    sage: M = Manifold(1, 'C*', field='complex'); M
    1-dimensional complex manifold C*

We introduce a first open subset, which is actually
`\CC = \CC^*\setminus\{\infty\}` if we interpret `\CC^*` as the Alexandroff
one-point compactification of `\CC`::

    sage: U = M.open_subset('U')

A natural chart on `U` is then nothing but the identity map of `\CC`, hence
we denote the associated coordinate by `z`::

    sage: Z.<z> = U.chart()

The origin of the complex plane is the point of coordinate `z=0`::

    sage: O = U.point((0,), chart=Z, name='O'); O
    Point O on the 1-dimensional complex manifold C*

Another open subset of `\CC^*` is `V = \CC^*\setminus\{O\}`::

    sage: V = M.open_subset('V')

We define a chart on `V` such that the point at infinity is the point of
coordinate 0 in this chart::

    sage: W.<w> = V.chart(); W
    Chart (V, (w,))
    sage: inf = M.point((0,), chart=W, name='inf', latex_name=r'\infty')
    sage: inf
    Point inf on the 1-dimensional complex manifold C*

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
    Point i on the 1-dimensional complex manifold C*

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

    sage: M.list_of_subsets()
    [Open subset A of the 1-dimensional complex manifold C*,
     1-dimensional complex manifold C*,
     Open subset U of the 1-dimensional complex manifold C*,
     Open subset V of the 1-dimensional complex manifold C*]
    sage: M.atlas()
    [Chart (U, (z,)), Chart (V, (w,)), Chart (A, (z,)), Chart (A, (w,))]

A constant map `\CC^* \rightarrow \CC`::

    sage: f = M.constant_scalar_field(3+2*I, name='f'); f
    Scalar field f on the 1-dimensional complex manifold C*
    sage: f.display()
    f: C* --> C
    on U: z |--> 2*I + 3
    on V: w |--> 2*I + 3
    sage: f(O)
    2*I + 3
    sage: f(i)
    2*I + 3
    sage: f(inf)
    2*I + 3
    sage: f.parent()
    Algebra of differentiable scalar fields on the 1-dimensional complex
     manifold C*
    sage: f.parent().category()
    Category of commutative algebras over Symbolic Ring

A vector field on the Riemann sphere::

    sage: v = M.vector_field(name='v')
    sage: v[Z.frame(), 0] = z^2
    sage: v.add_comp_by_continuation(W.frame(), U.intersection(V), W)
    sage: v.display(Z.frame())
    v = z^2 d/dz
    sage: v.display(W.frame())
    v = -d/dw
    sage: v.parent()
    Module X(C*) of vector fields on the 1-dimensional complex manifold C*

The vector field `v` acting on the scalar field `f`::

    sage: v(f)
    Scalar field v(f) on the 1-dimensional complex manifold C*

Since `f` is constant, `v(f)` is vanishing::

    sage: v(f).display()
    v(f): C* --> C
    on U: z |--> 0
    on V: w |--> 0

The value of the vector field `v` at the point `\infty` is a vector tangent to
the Riemann sphere::

    sage: v.at(inf)
    Tangent vector v at Point inf on the 1-dimensional complex manifold C*
    sage: v.at(inf).display()
    v = -d/dw
    sage: v.at(inf).parent()
    Tangent space at Point inf on the 1-dimensional complex manifold C*

AUTHORS:

- Eric Gourgoulhon (2015): initial version
- Travis Scrimshaw (2016): review tweaks

REFERENCES:

- [Lee2013]_
- [KN1963]_
- [Huy2005]_
- [Ser1992]_
- [Ber2008]_
- [BG1988]_

"""

#*****************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#       Copyright (C) 2016 Travis Scrimshaw <tscrimsh@umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.categories.manifolds import Manifolds
from sage.categories.homset import Hom
from sage.rings.all import CC
from sage.rings.real_mpfr import RR
from sage.rings.infinity import infinity, minus_infinity
from sage.rings.integer import Integer
from sage.misc.latex import latex
from sage.manifolds.manifold import TopologicalManifold

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
    - ``ambient`` -- (default: ``None``) if not ``None``, must be a
      differentiable manifold; the created object is then an open subset of
      ``ambient``
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
    def __init__(self, n, name, field, structure, ambient=None,
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
        if ambient is None:
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
        elif not isinstance(ambient, DifferentiableManifold):
            raise TypeError("the argument 'ambient' must be a " +
                            "differentiable manifold")
        TopologicalManifold.__init__(self, n, name, field, structure,
                                     ambient=ambient,
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
        self._vector_field_modules = {}

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

    def open_subset(self, name, latex_name=None, coord_def={}):
        r"""
        Create an open subset of the manifold.

        An open subset is a set that is (i) included in the manifold and (ii)
        open with respect to the manifold's topology. It is a differentiable
        manifold by itself. Hence the returned object is an instance of
        :class:`DifferentiableManifold`.

        INPUT:

        - ``name`` -- name given to the open subset
        - ``latex_name`` --  (default: ``None``) LaTeX symbol to denote the
          subset; if none is provided, it is set to ``name``
        - ``coord_def`` -- (default: {}) definition of the subset in
          terms of coordinates; ``coord_def`` must a be dictionary with keys
          charts in the manifold's atlas and values the symbolic expressions
          formed by the coordinates to define the subset.

        OUTPUT:

        - the open subset, as an instance of :class:`DifferentiableManifold`

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

            sage: A.list_of_subsets()
            [Open subset A of the 2-dimensional differentiable manifold M,
             Open subset B of the 2-dimensional differentiable manifold M]
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

        An point in ``U``::

            sage: p = U.an_element(); p
            Point on the 2-dimensional differentiable manifold M
            sage: X(p)  # the coordinates (x,y) of p
            (0, 0)
            sage: p in U
            True

        Checking whether various points, defined by their coordinates w.r.t.
        chart ``X``,  are in ``U``::

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
                                   self._structure, ambient=self._manifold,
                                   diff_degree=self._diff_degree,
                                   latex_name=latex_name,
                                   start_index=self._sindex)
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
        #!# update vector frames and change of frames
        return resu

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

    def diff_mapping(self, codomain, coord_functions=None, chart1=None,
                     chart2=None, name=None, latex_name=None):
        r"""
        Deprecated.

        Use :meth:`diff_map` instead.

        EXAMPLE::

            sage: M = Manifold(2, 'M'); X.<x,y> = M.chart()
            sage: N = Manifold(2, 'N'); Y.<u,v> = N.chart()
            sage: Phi = M.diff_mapping(N, {(X,Y): [x+y, x-y]}, name='Phi')
            doctest:...: DeprecationWarning: Use diff_map() instead.
            See http://trac.sagemath.org/18783 for details.
            sage: Phi
            Differentiable map Phi from the 2-dimensional differentiable
             manifold M to the 2-dimensional differentiable manifold N

        """
        from sage.misc.superseded import deprecation
        deprecation(18783, 'Use diff_map() instead.')
        return self.diff_map(codomain, coord_functions=coord_functions,
                             chart1=chart1, chart2=chart2, name=name,
                             latex_name=latex_name)

    def diffeomorphism(self, codomain, coord_functions=None, chart1=None,
                       chart2=None, name=None, latex_name=None):
        r"""
        Define a diffeomorphism between the current manifold and another one.

        See :class:`~sage.manifolds.differentiable.diff_map.DiffMap` for a
        complete documentation.

        INPUT:

        - ``codomain`` -- codomain of the diffeomorphism (the arrival manifold
          or some subset of it)
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

        EXAMPLE:

        Diffeomorphism between the open unit disk in `\RR^2` and `\RR^2`::

            sage: M = Manifold(2, 'M')  # the open unit disk
            sage: forget()  # for doctests only
            sage: c_xy.<x,y> = M.chart('x:(-1,1) y:(-1,1)')  # Cartesian coord on M
            sage: c_xy.add_restrictions(x^2+y^2<1)
            sage: N = Manifold(2, 'N')  # R^2
            sage: c_XY.<X,Y> = N.chart()  # canonical coordinates on R^2
            sage: Phi = M.diffeomorphism(N, [x/sqrt(1-x^2-y^2), y/sqrt(1-x^2-y^2)],
            ....:                        name='Phi', latex_name=r'\Phi')
            sage: Phi
            Diffeomorphism Phi from the 2-dimensional differentiable manifold M
             to the 2-dimensional differentiable manifold N
            sage: Phi.display()
            Phi: M --> N
               (x, y) |--> (X, Y) = (x/sqrt(-x^2 - y^2 + 1), y/sqrt(-x^2 - y^2 + 1))

        The inverse diffeomorphism::

            sage: Phi^(-1)
            Diffeomorphism Phi^(-1) from the 2-dimensional differentiable
             manifold N to the 2-dimensional differentiable manifold M
            sage: (Phi^(-1)).display()
            Phi^(-1): N --> M
               (X, Y) |--> (x, y) = (X/sqrt(X^2 + Y^2 + 1), Y/sqrt(X^2 + Y^2 + 1))

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
        return homset(coord_functions, name=name, latex_name=latex_name,
                      is_isomorphism=True)

    def vector_field_module(self, dest_map=None, force_free=False):
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
        - ``force_free`` -- (default: ``False``) if set to ``True``, force
          the construction of a *free* module (this implies that `N` is
          parallelizable)

        OUTPUT:

        - a
          :class:`~sage.manifolds.differentiable.vectorfield_module.VectorFieldModule`
          (or if `N` is parallelizable, a
          :class:`~sage.manifolds.differentiable.vectorfield_module.VectorFieldFreeModule`)
          representing the module `\mathcal{X}(M,\Phi)` of vector fields on
          `M` taking values on `\Phi(M)\subset N`

        EXAMPLES:

        Vector field module `\mathcal{X}(U) := \mathcal{X}(U,\mathrm{Id}_U)`
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

        `\mathcal{X}(U)` is a free module because `U` is parallelizable
        (being a chart domain)::

            sage: U.is_manifestly_parallelizable()
            True

        Its rank is the manifold's dimension::

            sage: XU.rank()
            2

        The elements of `\mathcal{X}(U)` are vector fields on `U`::

            sage: XU.an_element()
            Vector field on the Open subset U of the 2-dimensional
             differentiable manifold S^2
            sage: XU.an_element().display()
            2 d/dth + 2 d/dph

        Vector field module `\mathcal{X}(U,\Phi)` of the
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

        `\mathcal{X}(U,\Phi)` is a free module because `\RR^3`
        is parallelizable and its rank is 3::

            sage: XU_R3.rank()
            3

        """
        from sage.manifolds.differentiable.vectorfield_module import \
                                       VectorFieldModule, VectorFieldFreeModule
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
            for complete documentation.

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
          tensor fields on `M` taking values on `\Phi(M)\subset M`

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
            Category of finite dimensional modules over Algebra of
             differentiable scalar fields on the Open subset U of the
             3-dimensional differentiable manifold M
            sage: TU.base_ring()
            Algebra of differentiable scalar fields on the Open subset U of
             the 3-dimensional differentiable manifold M
            sage: TU.base_ring() is U.scalar_field_algebra()
            True
            sage: TU.an_element()
            Tensor field of type (2,1) on the Open subset U of the
             3-dimensional differentiable manifold M
            sage: TU.an_element().display()
            2 d/dx*d/dx*dx

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
          representing the module `\Lambda^p(M,\Phi)` of `p`-forms on `M`
          taking values on `\Phi(M)\subset N`

        EXAMPLES:

        Module of 2-forms on a 3-dimensional parallelizable manifold::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: M.diff_form_module(2)
            Free module /\^2(M) of 2-forms on the 3-dimensional differentiable
             manifold M
            sage: M.diff_form_module(2).category()
            Category of finite dimensional modules over Algebra of
             differentiable scalar fields on the 3-dimensional differentiable
             manifold M
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

    def automorphism_field_group(self, dest_map=None):
        r"""
        Return the group of tangent-space automorphism fields defined on
        ``self``, possibly with values in another manifold, as a module
        over the algebra of scalar fields defined on ``self``.

        If `M` is the current manifold and `\Phi` a differentiable map
        `\Phi: M \rightarrow N`, where `N` is a differentiable manifold,
        this method called with ``dest_map`` being `\Phi` returns the
        general linear group `\mathrm{GL}(\mathcal{X}(M, \Phi))` of the module
        `\mathcal{X}(M, \Phi)` of vector fields along `M` with values in
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
          `\mathrm{GL}(\mathcal{X}(U, \Phi))`

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

    def vector_field(self, name=None, latex_name=None, dest_map=None):
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

        See :class:`~sage.manifolds.differentiable.vectorfield.VectorField`
        for a complete documentation.

        INPUT:

        - ``name`` -- (default: ``None``) name given to the vector field
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          vector field; if none is provided, the LaTeX symbol is set to
          ``name``
        - ``dest_map`` -- (default: ``None``) the destination map
          `\Phi:\ M \rightarrow N`; if ``None``, it is assumed that `N = M`
          andthat `\Phi` is the identity map (case of a vector field
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
            sage: v = U.vector_field('v'); v
            Vector field v on the Open subset U of the 3-dimensional
             differentiable manifold M

        The vector fields on `U` form the set `\mathcal{X}(U)`, which is a
        module over the algebra `C^k(U)` of differentiable scalar fields
        on `U`::

            sage: v.parent()
            Free module X(U) of vector fields on the Open subset U of the
             3-dimensional differentiable manifold M
            sage: v in U.vector_field_module()
            True

        .. SEEALSO::

            For more examples, see
            :class:`~sage.manifolds.differentiable.vectorfield.VectorField`.

        """
        vmodule = self.vector_field_module(dest_map)  # the parent
        return vmodule.element_class(vmodule, name=name, latex_name=latex_name)

    def tensor_field(self, k, l, name=None, latex_name=None, sym=None,
                     antisym=None, dest_map=None):
        r"""
        Define a tensor field on ``self``.

        Via the argument ``dest_map``, it is possible to let the tensor field
        take its values on another manifold. More precisely, if `M` is
        the current manifold, `N` a differentiable manifold,
        `\Phi:\  M \rightarrow N` a differentiable map and `(k,l)`
        a pair of non-negative integers, a *tensor field of type* `(k,l)`
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

        See :class:`~sage.manifolds.differentiable.tensorfield.TensorField`
        for a complete documentation.

        INPUT:

        - ``k`` -- the contravariant rank `k`, the tensor type being `(k,l)`
        - ``l`` -- the covariant rank `l`, the tensor type being `(k,l)`
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

        A tensor field of type `(2,0)` on an open subset of a 3-dimensional
        differentiable manifold::

            sage: M = Manifold(3, 'M')
            sage: U = M.open_subset('U')
            sage: c_xyz.<x,y,z> = U.chart()
            sage: t = U.tensor_field(2, 0, 'T'); t
            Tensor field T of type (2,0) on the Open subset U of the
             3-dimensional differentiable manifold M

        The type `(2,0)` tensor fields on `U` form the set
        `\mathcal{T}^{(2,0)}(U)`, which is a module over the algebra `C^k(U)`
        of differentiable scalar fields on `U`::

            sage: t.parent()
            Free module T^(2,0)(U) of type-(2,0) tensors fields on the Open
             subset U of the 3-dimensional differentiable manifold M
            sage: t in U.tensor_field_module((2,0))
            True

        .. SEEALSO:

            For more examples, see
            :class:`~sage.manifolds.differentiable.tensorfield.TensorField`.

        """
        vmodule = self.vector_field_module(dest_map)
        return vmodule.tensor((k,l), name=name, latex_name=latex_name,
                              sym=sym, antisym=antisym)

    def sym_bilin_form_field(self, name=None, latex_name=None, dest_map=None):
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

        A field of symmetric bilinear forms on a 3-dimensional manifold::

            sage: M = Manifold(3, 'M')
            sage: c_xyz.<x,y,z> = M.chart()
            sage: t = M.sym_bilin_form_field('T'); t
            Field of symmetric bilinear forms T on the 3-dimensional
             differentiable manifold M

        Such a object is a tensor field of rank 2 and type `(0,2)`::

            sage: t.parent()
            Free module T^(0,2)(M) of type-(0,2) tensors fields on the
             3-dimensional differentiable manifold M
            sage: t.tensor_rank()
            2
            sage: t.tensor_type()
            (0, 2)

        The LaTeX symbol is deduced from the name or can be specified when
        creating the object::

            sage: latex(t)
            T
            sage: om = M.sym_bilin_form_field('Omega', r'\Omega')
            sage: latex(om)
            \Omega

        Components with respect to some vector frame::

            sage: e = M.vector_frame('e') ; M.set_default_frame(e)
            sage: t.set_comp()
            Fully symmetric 2-indices components w.r.t. Vector frame
             (M, (e_0,e_1,e_2))
            sage: type(t.comp())
            <class 'sage.tensor.modules.comp.CompFullySym'>

        For the default frame, the components are accessed with the
        square brackets::

            sage: t[0,0], t[0,1], t[0,2] = (1, 2, 3)
            sage: t[1,1], t[1,2] = (4, 5)
            sage: t[2,2] = 6

        The other components are deduced by symmetry::

            sage: t[1,0], t[2,0], t[2,1]
            (2, 3, 5)
            sage: t[:]
            [1 2 3]
            [2 4 5]
            [3 5 6]

        A symmetric bilinear form acts on vector pairs::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: t = M.sym_bilin_form_field('T')
            sage: t[0,0], t[0,1], t[1,1] = (-1, x, y*x)
            sage: v1 = M.vector_field('V_1')
            sage: v1[:] = (y,x)
            sage: v2 = M.vector_field('V_2')
            sage: v2[:] = (x+y,2)
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

            sage: a = M.sym_bilin_form_field()
            sage: a[0,0], a[0,1], a[1,1] = (1,2,3)
            sage: b = M.sym_bilin_form_field()
            sage: b[0,0], b[0,1], b[1,1] = (-1,4,5)
            sage: s = a + b ; s
            Field of symmetric bilinear forms on the 2-dimensional
             differentiable manifold M
            sage: s[:]
            [0 6]
            [6 8]

        But adding a symmetric bilinear from with a non-symmetric bilinear
        form results in a generic type `(0,2)` tensor::

            sage: c = M.tensor_field(0,2)
            sage: c[:] = [[-2, -3], [1,7]]
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
        return self.tensor_field(0, 2, name=name, latex_name=latex_name,
                                 sym=(0,1))

    def diff_form(self, degree, name=None, latex_name=None,
                  dest_map=None):
        r"""
        Define a differential form on ``self``.

        Via the argument ``dest_map``, it is possible to let the
        differential form take its values on another manifold. More
        precisely, if `M` is the current manifold, `N` a differentiable
        manifold, `\Phi:\  M \rightarrow N` a differentiable map and `p`
        a non-negative integer, a *differential form of degree* `p` (or
        `p`-form) *along* `M` *with values on* `N` is a differentiable map

        .. MATH::

            t:\ M  \longrightarrow T^{(0,p)}N

        (`T^{(0,p)} N` being the tensor bundle of type `(0,p)` over `N`)
        such that

        .. MATH::

            \forall p \in M,\ t(p) \in \Lambda^p(T^*_{\Phi(p)} N),

        where `\Lambda^p(T^*_{\Phi(p)} N)` is the `p`-th exterior power
        of the dual of the tangent space `T_{\Phi(p)} N`.

        The standard case of a differential form *on* `M` corresponds
        to `N = M` and `\Phi = \mathrm{Id}_M`. Other common cases are `\Phi`
        being an immersion and `\Phi` being a curve in `N` (`M` is then
        an open interval of `\RR`).

        For `p = 1`, one can use the method
        :meth:`~sage.manifolds.differentiable.manifold.DifferentiableManifold.one_form`
        instead.

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.diff_form.DiffForm` for
            complete documentation.

        INPUT:

        - ``degree`` -- the degree `p` of the differential form (i.e. its
          tensor rank)
        - ``name`` -- (default: ``None``) name given to the differential form
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          differential form; if none is provided, the LaTeX symbol is set to
          ``name``
        - ``dest_map`` -- (default: ``None``) the destination map
          `\Phi:\ M \rightarrow N`; if ``None``, it is assumed that `N = M`
          and that `\Phi` is the identity map (case of a differential
          form *on* `M`), otherwise ``dest_map`` must be a
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`

        OUTPUT:

        - the `p`-form as a
          :class:`~sage.manifolds.differentiable.diff_form.DiffForm`
          (or if `N` is parallelizable, a
          :class:`~sage.manifolds.differentiable.diff_form.DiffFormParal`)

        EXAMPLES:

        A 2-form on a open subset of a 4-dimensional differentiable manifold::

            sage: M = Manifold(4, 'M')
            sage: A = M.open_subset('A', latex_name=r'\mathcal{A}'); A
            Open subset A of the 4-dimensional differentiable manifold M
            sage: c_xyzt.<x,y,z,t> = A.chart()
            sage: f = A.diff_form(2, 'F'); f
            2-form F on the Open subset A of the 4-dimensional differentiable
             manifold M

        See the documentation of class
        :class:`~sage.manifolds.differentiable.diff_form.DiffForm` for more
        examples.

        """
        vmodule = self.vector_field_module(dest_map)
        return vmodule.alternating_form(degree, name=name,
                                        latex_name=latex_name)

    def one_form(self, name=None, latex_name=None, dest_map=None):
        r"""
        Define a 1-form on the manifold.

        Via the argument ``dest_map``, it is possible to let the
        1-form take its values on another manifold. More precisely,
        if `M` is the current manifold, `N` a differentiable
        manifold and `\Phi:\ M \rightarrow N` a differentiable map,
        a *1-form along* `M` *with values on* `N` is a differentiable map

        .. MATH::

            t:\ M  \longrightarrow T^* N

        (`T^* N` being the cotangent bundle of `N`) such that

        .. MATH::

            \forall p \in M,\ t(p) \in T^*_{\Phi(p)}N,

        where `T^*_{\Phi(p)}` is the dual of the tangent space `T_{\Phi(p)} N`.

        The standard case of a 1-form *on* `M` corresponds to `N = M`
        and `\Phi = \mathrm{Id}_M`. Other common cases are `\Phi`
        being an immersion and `\Phi` being a curve in `N` (`M` is then
        an open interval of `\RR`).

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.diff_form.DiffForm` for
            complete documentation.

        INPUT:

        - ``name`` -- (default: ``None``) name given to the 1-form
        - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          1-form; if none is provided, the LaTeX symbol is set to ``name``
        - ``dest_map`` -- (default: ``None``) the destination map
          `\Phi:\ M \rightarrow N`; if ``None``, it is assumed that `N = M`
          and that `\Phi` is the identity map (case of a 1-form *on* `M`),
          otherwise ``dest_map`` must be a
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`

        OUTPUT:

        - the 1-form as a
          :class:`~sage.manifolds.differentiable.diff_form.DiffForm`
          (or if `N` is parallelizable, a
          :class:`~sage.manifolds.differentiable.diff_form.DiffFormParal`)

        EXAMPLES:

        A 1-form on a 3-dimensional open subset::

            sage: M = Manifold(3, 'M')
            sage: A = M.open_subset('A', latex_name=r'\mathcal{A}')
            sage: X.<x,y,z> = A.chart()
            sage: om = A.one_form('omega', r'\omega') ; om
            1-form omega on the Open subset A of the 3-dimensional
             differentiable manifold M
            sage: om.parent()
            Free module /\^1(A) of 1-forms on the Open subset A of the
             3-dimensional differentiable manifold M

        .. SEEALSO::

            For more examples, see
            :class:`~sage.manifolds.differentiable.diff_form.DiffForm`.

        """
        vmodule = self.vector_field_module(dest_map)
        return vmodule.linear_form(name=name, latex_name=latex_name)

    def automorphism_field(self, name=None, latex_name=None, dest_map=None):
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
            for complete documentation.

        INPUT:

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

        A field of automorphisms on a 3-dimensional manifold::

            sage: M = Manifold(3,'M')
            sage: c_xyz.<x,y,z> = M.chart()
            sage: a = M.automorphism_field('A') ; a
            Field of tangent-space automorphisms A on the 3-dimensional
             differentiable manifold M
            sage: a.parent()
            General linear group of the Free module X(M) of vector fields on
             the 3-dimensional differentiable manifold M

        .. SEEALSO::

            For more examples, see
            :class:`~sage.manifolds.differentiable.automorphismfield.AutomorphismField`.

        """
        vmodule = self.vector_field_module(dest_map)
        return vmodule.automorphism(name=name, latex_name=latex_name)

    def tangent_identity_field(self, name='Id', latex_name=None, dest_map=None):
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

        .. SEEALSO::

            For more examples, see
            :class:`~sage.manifolds.differentiable.automorphismfield.AutomorphismField`.

        """
        vmodule = self.vector_field_module(dest_map)
        return vmodule.identity_map(name=name, latex_name=latex_name)

    def default_frame(self):
        r"""
        Return the default vector frame defined on ``self``.

        By *vector frame*, it is meant a field on the manifold that provides,
        at each point `p`, a vector basis of the tangent space at `p`.

        Unless changed via :meth:`set_default_frame`, the default frame is
        the first one defined on the manifold, usually implicitely as the
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
            Coordinate frame (M, (d/dx,d/dy))

        """
        return self._def_frame

    def set_default_frame(self, frame):
        r"""
        Changing the default vector frame on ``self``.

        INPUT:

        - ``frame`` --
          :class:`~sage.manifolds.differentiable.vectorframe.VectorFrame`
          a vector frame defined on the manifold

        EXAMPLES:

        Changing the default frame on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: e = M.vector_frame('e')
            sage: M.default_frame()
            Coordinate frame (M, (d/dx,d/dy))
            sage: M.set_default_frame(e)
            sage: M.default_frame()
            Vector frame (M, (e_0,e_1))

        """
        from sage.manifolds.differentiable.vectorframe import VectorFrame
        if not isinstance(frame, VectorFrame):
            raise TypeError("{} is not a vector frame".format(frame))
        if frame._domain is not self:
            if self.is_manifestly_parallelizable():
                raise ValueError("the frame domain must coincide with " +
                                 "the {}".format(self))
            if not frame._domain.is_subset(self):
                raise ValueError("the frame must be defined on " +
                                 "the {}".format(self))
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
             - (M, (d/dx,d/dy)) (default basis)
             - (M, (d/du,d/dv))
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
        - ``compute_inverse`` (default: True) -- if set to True, the inverse
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
        from sage.manifolds.differentiable.automorphismfield import AutomorphismFieldParal
        fmodule = frame1._fmodule
        if frame2._fmodule != fmodule:
            raise ValueError("the two frames are not defined on the same " +
                             "vector field module")
        if not isinstance(change_of_frame, AutomorphismFieldParal):
            raise TypeError("the argument change_of_frame must be some " +
                            "instance of AutomorphismFieldParal")
        fmodule.set_change_of_basis(frame1, frame2, change_of_frame,
                                    compute_inverse=compute_inverse)
        for sdom in self._supersets:
            sdom._frame_changes[(frame1, frame2)] = change_of_frame
        if compute_inverse:
            if (frame2, frame1) not in self._frame_changes:
                for sdom in self._supersets:
                    sdom._frame_changes[(frame2, frame1)] = change_of_frame.inverse()

    def vector_frame(self, symbol=None, latex_symbol=None, dest_map=None,
                     from_frame=None):
        r"""
        Define a vector frame on ``self``.

        A *vector frame* is a field on the manifold that provides, at each
        point `p` of the manifold, a vector basis of the tangent space at `p`.

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.vectorframe.VectorFrame`
            for complete documentation.

        INPUT:

        - ``symbol`` -- (default: ``None``) a letter (of a few letters) to
          denote a generic vector of the frame; can be set to ``None`` if the
          parameter ``from_frame`` is filled
        - ``latex_symbol`` -- (default: ``None``) symbol to denote a generic
          vector of the frame; if None, the value of ``symbol`` is used.
        - ``dest_map`` -- (default: ``None``)
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`;
          destination map `\Phi:\ U \rightarrow V`; if ``None``, the
          identity is assumed (case of a vector frame *on* `U`)
        - ``from_frame`` -- (default: ``None``) vector frame `\tilde{e}`
          on the codomain `V` of the destination map `\Phi`; the returned
          frame `e` is then such that for all `p \in U`,
          we have `e(p) = \tilde{e}(\Phi(p))`

        OUTPUT:

        - a :class:`~sage.manifolds.differentiable.vectorframe.VectorFrame`
          representing the defined vector frame

        EXAMPLES:

        Setting a vector frame on a 3-dimensional open subset::

            sage: M = Manifold(3, 'M')
            sage: A = M.open_subset('A', latex_name=r'\mathcal{A}'); A
            Open subset A of the 3-dimensional differentiable manifold M
            sage: c_xyz.<x,y,z> = A.chart()
            sage: e = A.vector_frame('e'); e
            Vector frame (A, (e_0,e_1,e_2))
            sage: e[0]
            Vector field e_0 on the Open subset A of the 3-dimensional
             differentiable manifold M

        .. SEEALSO::

            For more examples, see
            :class:`~sage.manifolds.differentiable.vectorframe.VectorFrame`.

        """
        from sage.manifolds.differentiable.vectorframe import VectorFrame
        return VectorFrame(self.vector_field_module(dest_map=dest_map,
                                                    force_free=True),
                           symbol=symbol, latex_symbol=latex_symbol,
                           from_frame=from_frame)

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
        for sd in self._supersets:
            if not sd.is_manifestly_parallelizable():
                sd._parallelizable_parts.add(self)

    def frames(self):
        r"""
        Return the list of vector frames defined on open subsets of ``self``.

        OUTPUT:

        - list of vector frames defined on open subsets of ``self``

        EXAMPLES:

        Vector frames on subsets of `\RR^2`::

            sage: M = Manifold(2, 'R^2')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: M.frames()
            [Coordinate frame (R^2, (d/dx,d/dy))]
            sage: e = M.vector_frame('e')
            sage: M.frames()
            [Coordinate frame (R^2, (d/dx,d/dy)),
             Vector frame (R^2, (e_0,e_1))]
            sage: U = M.open_subset('U', coord_def={c_cart: x^2+y^2<1}) # unit disk
            sage: U.frames()
            [Coordinate frame (U, (d/dx,d/dy))]
            sage: M.frames()
            [Coordinate frame (R^2, (d/dx,d/dy)),
             Vector frame (R^2, (e_0,e_1)),
             Coordinate frame (U, (d/dx,d/dy))]

        """
        return self._frames

    def coframes(self):
        r"""
        Return the list of coframes defined on open subsets of ``self``.

        OUTPUT:

        - list of coframes defined on open subsets of ``self``

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
        return self._coframes

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
            Coordinate frame (M, (d/dx,d/dy))

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
            {(Coordinate frame (M, (d/dx,d/dy)),
              Vector frame (M, (f_0,f_1))): Field of tangent-space
               automorphisms on the 2-dimensional differentiable manifold M,
             (Vector frame (M, (f_0,f_1)),
              Coordinate frame (M, (d/dx,d/dy))): Field of tangent-space
               automorphisms on the 2-dimensional differentiable manifold M}

        Some checks::

            sage: M.changes_of_frame()[(e,f)] == a
            True
            sage: M.changes_of_frame()[(f,e)] == a^(-1)
            True

        """
        return self._frame_changes

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

    def tangent_space(self, point):
        r"""
        Tangent space to ``self`` at a given point.

        INPUT:

        - ``point`` -- :class:`~sage.manifolds.point.ManifoldPoint`;
          point `p` on the manifold

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
        from sage.manifolds.point import ManifoldPoint
        from sage.manifolds.differentiable.tangent_space import TangentSpace
        if not isinstance(point, ManifoldPoint):
            raise TypeError("{} is not a manifold point".format(point))
        if point not in self:
            raise ValueError("{} is not a point on the {}".format(point, self))
        return TangentSpace(point)

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
        - ``param`` -- a tuple of the type ``(t, t_min, t_max)``, where

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

        OUTPUT:

        - :class:`~sage.manifolds.differentiable.curve.DifferentiableCurve`

        EXAMPLES:

        The lemniscate of Gerono in the 2-dimensional Euclidean plane::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: R.<t> = RealLine()
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
            Real number line R

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.curve.DifferentiableCurve`
            for more examples, including plots.

        """
        from sage.manifolds.differentiable.real_line import RealLine
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

        EXAMPLE:

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
        from sage.manifolds.differentiable.affine_connection import \
                                                               AffineConnection
        return AffineConnection(self, name, latex_name)

    def metric(self, name, signature=None, latex_name=None, dest_map=None):
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
          (case of a metric field *on* `U`)

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`
          representing the defined pseudo-Riemannian metric.

        EXAMPLE:

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
          (case of a metric field *on* `U`)

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`
          representing the defined Riemannian metric.

        EXAMPLE:

        Metric of the hyperbolic plane `H^2`::

            sage: H2 = Manifold(2, 'H^2', start_index=1)
            sage: X.<x,y> = H2.chart('x y:(0,+oo)')  # Poincare half-plane coord.
            sage: g = H2.riemannian_metric('g')
            sage: g[1,1], g[2,2] = 1/y^2, 1/y^2
            sage: g
            Riemannian metric g on the 2-dimensional differentiable manifold H^2
            sage: g.display()
            g = y^(-2) dx*dx + y^(-2) dy*dy
            sage: g.signature()
            2

        .. SEEALSO::

            :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`
            for more examples.

        """
        vmodule = self.vector_field_module(dest_map)
        dim = vmodule.ambient_domain().dimension()
        return vmodule.metric(name, signature=dim, latex_name=latex_name)

    def riemann_metric(self, name, latex_name=None, dest_map=None):
        r"""
        Deprecated.

        Use :meth:`riemannian_metric` instead.

        EXAMPLE::

            sage: M = Manifold(3, 'M')
            sage: g = M.riemann_metric('g')
            doctest:...: DeprecationWarning: Use riemannian_metric() instead.
            See http://trac.sagemath.org/19209 for details.
            sage: g
            Riemannian metric g on the 3-dimensional differentiable manifold M

        """
        from sage.misc.superseded import deprecation
        deprecation(19209, 'Use riemannian_metric() instead.')
        return self.riemannian_metric(name, latex_name=latex_name,
                                      dest_map=dest_map)

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
        - ``signature`` -- (default: 'positive') sign of the metric
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
          (case of a metric field *on* `U`)

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric`
          representing the defined Lorentzian metric.

        EXAMPLE:

        Metric of Minkowski spacetime::

            sage: M = Manifold(4, 'M')
            sage: X.<t,x,y,z> = M.chart()
            sage: g = M.lorentzian_metric('g'); g
            Lorentzian metric g on the 4-dimensional differentiable manifold M
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = -1, 1, 1, 1
            sage: g.display()
            g = -dt*dt + dx*dx + dy*dy + dz*dz
            sage: g.signature()
            2

        Choice of a negative signature::

            sage: g = M.lorentzian_metric('g', signature='negative'); g
            Lorentzian metric g on the 4-dimensional differentiable manifold M
            sage: g[0,0], g[1,1], g[2,2], g[3,3] = 1, -1, -1, -1
            sage: g.display()
            g = dt*dt - dx*dx - dy*dy - dz*dz
            sage: g.signature()
            -2

        """
        vmodule = self.vector_field_module(dest_map)
        dim = vmodule.ambient_domain().dimension()
        if signature=='positive':
            signat = dim - 2
        else:
            signat = 2 - dim
        return vmodule.metric(name, signature=signat, latex_name=latex_name)

    def lorentz_metric(self, name, signature='positive', latex_name=None,
                       dest_map=None):
        r"""
        Deprecated.

        Use :meth:`lorentzian_metric` instead.

        EXAMPLE::

            sage: M = Manifold(4, 'M')
            sage: g = M.lorentz_metric('g')
            doctest:...: DeprecationWarning: Use lorentzian_metric() instead.
            See http://trac.sagemath.org/19209 for details.
            sage: g
            Lorentzian metric g on the 4-dimensional differentiable manifold M

        """
        from sage.misc.superseded import deprecation
        deprecation(19209, 'Use lorentzian_metric() instead.')
        return self.lorentzian_metric(name, signature=signature,
                                      latex_name=latex_name, dest_map=dest_map)
