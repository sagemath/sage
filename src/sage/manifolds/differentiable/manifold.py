r"""
Differentiable Manifolds

Given a non-discrete topological field `K` (in most applications, `K = \RR` or
`K = \CC`; see however [4]_ for `K = \QQ_p` and [5]_ for other fields),
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

AUTHORS:

- Eric Gourgoulhon (2015): initial version

REFERENCES:

.. [1] \J.M. Lee : *Introduction to Smooth Manifolds*, 2nd ed., Springer
   (New York) (2012); :doi:`10.1007/978-1-4419-9982-5`
.. [2] \S. Kobayashi & K. Nomizu : *Foundations of Differential Geometry*,
   vol. 1, Interscience Publishers (New York) (1963)
.. [3] \D. Huybrechts : *Complex Geometry*, Springer (Berlin) (2005);
   :doi:`10.1007/b137952`
.. [4] \J.-P. Serre : *Lie Algebras and Lie Groups*, 2nd ed., Springer
   (Berlin) (1992); :doi:`10.1007/978-3-540-70634-2`
.. [5] \W. Bertram : *Differential Geometry, Lie Groups and Symmetric Spaces
   over General Base Fields and Rings*, Memoirs of the American Mathematical
   Society, vol. 192 (2008); :doi:`10.1090/memo/0900`; :arxiv:`math/0502168`
.. [6] \M. Berger & B. Gostiaux : *Differential Geometry: Manifolds, Curves and
   Surfaces*, Springer (New York) (1988); :doi:`10.1007/978-1-4612-1033-7`

"""

#*****************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
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
from sage.rings.infinity import infinity
from sage.rings.integer import Integer
from sage.misc.latex import latex
from sage.manifolds.manifold import TopologicalManifold

###############################################################################

class DifferentiableManifold(TopologicalManifold):
    r"""
    Differentiable manifold over a topological field `K`.

    Given a non-discrete topological field `K` (in most applications,
    `K = \RR` or `K = \CC`; see however [4]_ for `K = \QQ_p` and [5]_ for
    other fields), a *differentiable manifold over* `K` is a topological
    manifold `M` over `K` equipped with an atlas whose transitions maps are of
    class `C^k` (i.e. `k`-times  continuously differentiable) for a fixed
    positive integer `k` (possibly `k=\infty`). `M` is then called a
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
        for chart, restrictions in coord_def.iteritems():
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

