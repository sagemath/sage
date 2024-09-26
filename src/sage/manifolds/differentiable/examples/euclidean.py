r"""
Euclidean Spaces

A *Euclidean space of dimension* `n` is an affine space `E`, whose associated
vector space is a `n`-dimensional vector space over `\RR` and is equipped with
a positive definite symmetric bilinear form, called the *scalar product* or
*dot product* [Ber1987]_. A Euclidean space of dimension `n` can also be
viewed as a Riemannian manifold that is diffeomorphic to `\RR^n` and that
has a flat metric `g`. The Euclidean scalar product is then that defined
by the Riemannian metric `g`.

The current implementation of Euclidean spaces is based on the second point of
view. This allows for the introduction of various coordinate systems in
addition to the usual the Cartesian systems. Standard curvilinear systems
(planar, spherical and cylindrical coordinates) are predefined for
2-dimensional and 3-dimensional Euclidean spaces, along with the corresponding
transition maps between them. Another benefit of such an implementation is
the direct use of methods for vector calculus already implemented at the
level of Riemannian manifolds (see, e.g., the methods
:meth:`~sage.manifolds.differentiable.vectorfield.VectorField.cross_product`
and
:meth:`~sage.manifolds.differentiable.vectorfield.VectorField.curl`,
as well as the module :mod:`~sage.manifolds.operators`).

Euclidean spaces are implemented via the following classes:

- :class:`EuclideanSpace` for generic values `n`,
- :class:`EuclideanPlane` for `n = 2`,
- :class:`Euclidean3dimSpace` for `n = 3`.

The user interface is provided by :class:`EuclideanSpace`.

.. _EuclideanSpace_example1:

.. RUBRIC:: Example 1: the Euclidean plane

We start by declaring the Euclidean plane ``E``, with ``(x, y)`` as
Cartesian coordinates::

    sage: E.<x,y> = EuclideanSpace()
    sage: E
    Euclidean plane E^2
    sage: dim(E)
    2

``E`` is automatically endowed with the chart of Cartesian coordinates::

    sage: E.atlas()
    [Chart (E^2, (x, y))]
    sage: cartesian = E.default_chart(); cartesian
    Chart (E^2, (x, y))

Thanks to the use of ``<x,y>`` when declaring ``E``, the coordinates `(x,y)`
have been injected in the global namespace, i.e. the Python variables ``x``
and ``y`` have been created and are available to form symbolic expressions::

    sage: y
    y
    sage: type(y)
    <class 'sage.symbolic.expression.Expression'>
    sage: assumptions()
    [x is real, y is real]

The metric tensor of ``E`` is predefined::

    sage: g = E.metric(); g
    Riemannian metric g on the Euclidean plane E^2
    sage: g.display()
    g = dx⊗dx + dy⊗dy
    sage: g[:]
    [1 0]
    [0 1]

It is a *flat* metric, i.e. it has a vanishing Riemann tensor::

    sage: g.riemann()
    Tensor field Riem(g) of type (1,3) on the Euclidean plane E^2
    sage: g.riemann().display()
    Riem(g) = 0

Polar coordinates `(r,\phi)` are introduced by::

    sage: polar.<r,ph> = E.polar_coordinates()
    sage: polar
    Chart (E^2, (r, ph))

``E`` is now endowed with two coordinate charts::

    sage: E.atlas()
    [Chart (E^2, (x, y)), Chart (E^2, (r, ph))]

The ranges of the coordinates introduced so far are::

    sage: cartesian.coord_range()
    x: (-oo, +oo); y: (-oo, +oo)
    sage: polar.coord_range()
    r: (0, +oo); ph: [0, 2*pi] (periodic)

The transition map from polar coordinates to Cartesian ones is::

    sage: E.coord_change(polar, cartesian).display()
    x = r*cos(ph)
    y = r*sin(ph)

while the reverse one is::

    sage: E.coord_change(cartesian, polar).display()
    r = sqrt(x^2 + y^2)
    ph = arctan2(y, x)

A point of ``E`` is constructed from its coordinates (by default in the
Cartesian chart)::

    sage: p = E((-1,1), name='p'); p
    Point p on the Euclidean plane E^2
    sage: p.parent()
    Euclidean plane E^2

The coordinates of a point are obtained by letting the corresponding chart
act on it::

    sage: cartesian(p)
    (-1, 1)
    sage: polar(p)
    (sqrt(2), 3/4*pi)

At this stage, ``E`` is endowed with three vector frames::

    sage: E.frames()
    [Coordinate frame (E^2, (e_x,e_y)),
     Coordinate frame (E^2, (∂/∂r,∂/∂ph)),
     Vector frame (E^2, (e_r,e_ph))]

The third one is the standard orthonormal frame associated with polar
coordinates, as we can check from the metric components in it::

    sage: polar_frame = E.polar_frame(); polar_frame
    Vector frame (E^2, (e_r,e_ph))
    sage: g[polar_frame,:]
    [1 0]
    [0 1]

The expression of the metric tensor in terms of polar coordinates is::

    sage: g.display(polar)
    g = dr⊗dr + r^2 dph⊗dph

A vector field on ``E``::

    sage: v = E.vector_field(-y, x, name='v'); v
    Vector field v on the Euclidean plane E^2
    sage: v.display()
    v = -y e_x + x e_y
    sage: v[:]
    [-y, x]

By default, the components of ``v``, as returned by ``display`` or the bracket
operator, refer to the Cartesian frame on ``E``; to get the components with
respect to the orthonormal polar frame, one has to specify it explicitly,
generally along with the polar chart for the coordinate expression of the
components::

    sage: v.display(polar_frame, polar)
    v = r e_ph
    sage: v[polar_frame,:,polar]
    [0, r]

Note that the default frame for the display of vector fields can be changed
thanks to the method
:meth:`~sage.manifolds.differentiable.manifold.DifferentiableManifold.set_default_frame`;
in the same vein, the default coordinates can be changed via the method
:meth:`~sage.manifolds.manifold.TopologicalManifold.set_default_chart`::

    sage: E.set_default_frame(polar_frame)
    sage: E.set_default_chart(polar)
    sage: v.display()
    v = r e_ph
    sage: v[:]
    [0, r]
    sage: E.set_default_frame(E.cartesian_frame())  # revert to Cartesian frame
    sage: E.set_default_chart(cartesian)            # and chart

When defining a vector field from components relative to a vector frame
different from the default one, the vector frame has to be specified
explicitly::

    sage: v = E.vector_field(1, 0, frame=polar_frame)
    sage: v.display(polar_frame)
    e_r
    sage: v.display()
    x/sqrt(x^2 + y^2) e_x + y/sqrt(x^2 + y^2) e_y

The argument ``chart`` must be used to specify in which coordinate
chart the components are expressed::

    sage: v = E.vector_field(0, r, frame=polar_frame, chart=polar)
    sage: v.display(polar_frame, polar)
    r e_ph
    sage: v.display()
    -y e_x + x e_y

It is also possible to pass the components as a dictionary, with
a pair (vector frame, chart) as a key::

    sage: v = E.vector_field({(polar_frame, polar): (0, r)})
    sage: v.display(polar_frame, polar)
    r e_ph

The key can be reduced to the vector frame if the chart is the default
one::

    sage: v = E.vector_field({polar_frame: (0, 1)})
    sage: v.display(polar_frame)
    e_ph

Finally, it is possible to construct the vector field without
initializing any component::

    sage: v = E.vector_field(); v
    Vector field on the Euclidean plane E^2

The components can then by set in a second stage, via the square
bracket operator, the unset components being assumed to be zero::

    sage: v[1] = -y
    sage: v.display()  # v[2] is zero
    -y e_x
    sage: v[2] = x
    sage: v.display()
    -y e_x + x e_y

The above is equivalent to::

    sage: v[:] = -y, x
    sage: v.display()
    -y e_x + x e_y

The square bracket operator can also be used to set components in a
vector frame that is not the default one::

    sage: v = E.vector_field(name='v')
    sage: v[polar_frame, 2, polar] = r
    sage: v.display(polar_frame, polar)
    v = r e_ph
    sage: v.display()
    v = -y e_x + x e_y

The value of the vector field ``v`` at point ``p``::

    sage: vp = v.at(p); vp
    Vector v at Point p on the Euclidean plane E^2
    sage: vp.display()
    v = -e_x - e_y
    sage: vp.display(polar_frame.at(p))
    v = sqrt(2) e_ph

A scalar field on ``E``::

    sage: f = E.scalar_field(x*y, name='f'); f
    Scalar field f on the Euclidean plane E^2
    sage: f.display()
    f: E^2 → ℝ
       (x, y) ↦ x*y
       (r, ph) ↦ r^2*cos(ph)*sin(ph)

The value of ``f`` at point ``p``::

    sage: f(p)
    -1

The gradient of ``f``::

    sage: from sage.manifolds.operators import * # to get grad, div, etc.
    sage: w = grad(f); w
    Vector field grad(f) on the Euclidean plane E^2
    sage: w.display()
    grad(f) = y e_x + x e_y
    sage: w.display(polar_frame, polar)
    grad(f) = 2*r*cos(ph)*sin(ph) e_r + (2*cos(ph)^2 - 1)*r e_ph

The dot product of two vector fields::

    sage: s = v.dot(w); s
    Scalar field v.grad(f) on the Euclidean plane E^2
    sage: s.display()
    v.grad(f): E^2 → ℝ
       (x, y) ↦ x^2 - y^2
       (r, ph) ↦ (2*cos(ph)^2 - 1)*r^2
    sage: s.expr()
    x^2 - y^2

The norm is related to the dot product by the standard formula::

    sage: norm(v)^2 == v.dot(v)
    True

The divergence of the vector field ``v``::

    sage: s = div(v); s
    Scalar field div(v) on the Euclidean plane E^2
    sage: s.display()
    div(v): E^2 → ℝ
       (x, y) ↦ 0
       (r, ph) ↦ 0

.. _EuclideanSpace_example2:

.. RUBRIC:: Example 2: Vector calculus in the Euclidean 3-space

We start by declaring the 3-dimensional Euclidean space ``E``, with
``(x,y,z)`` as Cartesian coordinates::

    sage: E.<x,y,z> = EuclideanSpace()
    sage: E
    Euclidean space E^3

A simple vector field on ``E``::

    sage: v = E.vector_field(-y, x, 0, name='v')
    sage: v.display()
    v = -y e_x + x e_y
    sage: v[:]
    [-y, x, 0]

The Euclidean norm of ``v``::

    sage: s = norm(v); s
    Scalar field |v| on the Euclidean space E^3
    sage: s.display()
    |v|: E^3 → ℝ
       (x, y, z) ↦ sqrt(x^2 + y^2)
    sage: s.expr()
    sqrt(x^2 + y^2)

The divergence of ``v`` is zero::

    sage: from sage.manifolds.operators import *
    sage: div(v)
    Scalar field div(v) on the Euclidean space E^3
    sage: div(v).display()
    div(v): E^3 → ℝ
       (x, y, z) ↦ 0

while its curl is a constant vector field along `e_z`::

    sage: w = curl(v); w
    Vector field curl(v) on the Euclidean space E^3
    sage: w.display()
    curl(v) = 2 e_z

The gradient of a scalar field::

    sage: f = E.scalar_field(sin(x*y*z), name='f')
    sage: u = grad(f); u
    Vector field grad(f) on the Euclidean space E^3
    sage: u.display()
    grad(f) = y*z*cos(x*y*z) e_x + x*z*cos(x*y*z) e_y + x*y*cos(x*y*z) e_z

The curl of a gradient is zero::

    sage: curl(u).display()
    curl(grad(f)) = 0

The dot product of two vector fields::

    sage: s = u.dot(v); s
    Scalar field grad(f).v on the Euclidean space E^3
    sage: s.expr()
    (x^2 - y^2)*z*cos(x*y*z)

The cross product of two vector fields::

    sage: a = u.cross(v); a
    Vector field grad(f) x v on the Euclidean space E^3
    sage: a.display()
    grad(f) x v = -x^2*y*cos(x*y*z) e_x - x*y^2*cos(x*y*z) e_y
     + 2*x*y*z*cos(x*y*z) e_z

The scalar triple product of three vector fields::

    sage: triple_product = E.scalar_triple_product()
    sage: s = triple_product(u, v, w); s
    Scalar field epsilon(grad(f),v,curl(v)) on the Euclidean space E^3
    sage: s.expr()
    4*x*y*z*cos(x*y*z)

Let us check that the scalar triple product of `u`, `v` and `w` is
`u\cdot(v\times w)`::

    sage: s == u.dot(v.cross(w))
    True

AUTHORS:

- Eric Gourgoulhon (2018): initial version

REFERENCES:

- \M. Berger: *Geometry I* [Ber1987]_
"""

#*****************************************************************************
#       Copyright (C) 2018 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from sage.functions.trig import cos, sin, atan2
from sage.misc.functional import sqrt
from sage.misc.latex import latex
from sage.rings.real_mpfr import RR
from sage.categories.manifolds import Manifolds
from sage.categories.metric_spaces import MetricSpaces
from sage.manifolds.differentiable.pseudo_riemannian import \
                                                       PseudoRiemannianManifold

###############################################################################

class EuclideanSpace(PseudoRiemannianManifold):
    r"""
    Euclidean space.

    A *Euclidean space of dimension* `n` is an affine space `E`, whose
    associated vector space is a `n`-dimensional vector space over `\RR` and
    is equipped with a positive definite symmetric bilinear form, called
    the *scalar product* or *dot product*.

    Euclidean space of dimension `n` can be viewed as a Riemannian manifold
    that is diffeomorphic to `\RR^n` and that has a flat metric `g`. The
    Euclidean scalar product is the one defined by the Riemannian metric `g`.

    INPUT:

    - ``n`` -- positive integer; dimension of the space over the real field
    - ``name`` -- (default: ``None``) string; name (symbol) given to the
      Euclidean space; if ``None``, the name will be set to ``'E^n'``
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the space; if ``None``, it is set to ``'\mathbb{E}^{n}'`` if
      ``name`` is  ``None`` and to ``name`` otherwise
    - ``coordinates`` -- (default: ``'Cartesian'``) string describing the
      type of coordinates to be initialized at the Euclidean space
      creation; allowed values are

      - ``'Cartesian'`` (canonical coordinates on `\RR^n`)
      - ``'polar'`` for ``n=2`` only (see
        :meth:`~sage.manifolds.differentiable.examples.euclidean.EuclideanPlane.polar_coordinates`)
      - ``'spherical'`` for ``n=3`` only (see
        :meth:`~sage.manifolds.differentiable.examples.euclidean.Euclidean3dimSpace.spherical_coordinates`)
      - ``'cylindrical'`` for ``n=3`` only (see
        :meth:`~sage.manifolds.differentiable.examples.euclidean.Euclidean3dimSpace.cylindrical_coordinates`)

    - ``symbols`` -- (default: ``None``) string defining the coordinate
      text symbols and LaTeX symbols, with the same conventions as the
      argument ``coordinates`` in
      :class:`~sage.manifolds.differentiable.chart.RealDiffChart`, namely
      ``symbols`` is a string of coordinate fields separated by a blank
      space, where each field contains the coordinate's text symbol and
      possibly the coordinate's LaTeX symbol (when the latter is different
      from the text symbol), both symbols being separated by a colon
      (``:``); if ``None``, the symbols will be automatically generated
      according to the value of ``coordinates``
    - ``metric_name`` -- (default: ``'g'``) string; name (symbol) given to the
      Euclidean metric tensor
    - ``metric_latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the Euclidean metric tensor; if none is provided, it is set to
      ``metric_name``
    - ``start_index`` -- (default: 1) integer; lower value of the range of
      indices used for "indexed objects" in the Euclidean space, e.g.
      coordinates of a chart
    - ``names`` -- (default: ``None``) unused argument, except if
      ``symbols`` is not provided; it must then be a tuple containing
      the coordinate symbols (this is guaranteed if the shortcut operator
      ``<,>`` is used)

    If ``names`` is specified, then ``n`` does not have to be specified.

    EXAMPLES:

    Constructing a 2-dimensional Euclidean space::

        sage: E = EuclideanSpace(2); E
        Euclidean plane E^2

    Each call to :class:`EuclideanSpace` creates a different object::

        sage: E1 = EuclideanSpace(2)
        sage: E1 is E
        False
        sage: E1 == E
        False

    The LaTeX symbol of the Euclidean space is by default `\mathbb{E}^n`,
    where `n` is the dimension::

        sage: latex(E)
        \mathbb{E}^{2}

    But both the name and LaTeX names of the Euclidean space can be
    customized::

        sage: F = EuclideanSpace(2, name='F', latex_name=r'\mathcal{F}'); F
        Euclidean plane F
        sage: latex(F)
        \mathcal{F}

    By default, a Euclidean space is created with a single coordinate chart:
    that of Cartesian coordinates::

        sage: E.atlas()
        [Chart (E^2, (x, y))]
        sage: E.cartesian_coordinates()
        Chart (E^2, (x, y))
        sage: E.default_chart() is E.cartesian_coordinates()
        True

    The coordinate variables can be initialized, as the Python variables
    ``x`` and ``y``, by::

        sage: x, y = E.cartesian_coordinates()[:]

    However, it is possible to both construct the Euclidean space and
    initialize the coordinate variables in a single stage, thanks to
    SageMath operator ``<,>``::

        sage: E.<x,y> = EuclideanSpace()

    Note that providing the dimension as an argument of ``EuclideanSpace`` is
    not necessary in that case, since it can be deduced from the number of
    coordinates within ``<,>``. Besides, the coordinate symbols can be
    customized::

        sage: E.<X,Y> = EuclideanSpace()
        sage: E.cartesian_coordinates()
        Chart (E^2, (X, Y))

    By default, the LaTeX symbols of the coordinates coincide with the text
    ones::

        sage: latex(X+Y)
        X + Y

    However, it is possible to customize them, via the argument ``symbols``,
    which must be a string, usually prefixed by ``r`` (for *raw* string, in
    order to allow for the backslash character of LaTeX expressions). This
    string contains the coordinate fields separated by a blank space; each
    field contains the coordinate's text symbol and possibly the coordinate's
    LaTeX symbol (when the latter is different from the text symbol), both
    symbols being separated by a colon (``:``)::

        sage: E.<xi,ze> = EuclideanSpace(symbols=r"xi:\xi ze:\zeta")
        sage: E.cartesian_coordinates()
        Chart (E^2, (xi, ze))
        sage: latex(xi+ze)
        {\xi} + {\zeta}

    Thanks to the argument ``coordinates``, a Euclidean space can be
    constructed with curvilinear coordinates initialized instead of the
    Cartesian ones::

        sage: E.<r,ph> = EuclideanSpace(coordinates='polar')
        sage: E.atlas()   # no Cartesian coordinates have been constructed
        [Chart (E^2, (r, ph))]
        sage: polar = E.polar_coordinates(); polar
        Chart (E^2, (r, ph))
        sage: E.default_chart() is polar
        True
        sage: latex(r+ph)
        {\phi} + r

    The Cartesian coordinates, along with the transition maps to and from
    the curvilinear coordinates, can be constructed at any time by::

        sage: cartesian.<x,y> = E.cartesian_coordinates()
        sage: E.atlas()  # both polar and Cartesian coordinates now exist
        [Chart (E^2, (r, ph)), Chart (E^2, (x, y))]

    The transition maps have been initialized by the command
    ``E.cartesian_coordinates()``::

        sage: E.coord_change(polar, cartesian).display()
        x = r*cos(ph)
        y = r*sin(ph)
        sage: E.coord_change(cartesian, polar).display()
        r = sqrt(x^2 + y^2)
        ph = arctan2(y, x)

    The default name of the Euclidean metric tensor is `g`::

        sage: E.metric()
        Riemannian metric g on the Euclidean plane E^2
        sage: latex(_)
        g

    But this can be customized::

        sage: E = EuclideanSpace(2, metric_name='h')
        sage: E.metric()
        Riemannian metric h on the Euclidean plane E^2
        sage: latex(_)
        h
        sage: E = EuclideanSpace(2, metric_latex_name=r'\mathbf{g}')
        sage: E.metric()
        Riemannian metric g on the Euclidean plane E^2
        sage: latex(_)
        \mathbf{g}

    A 4-dimensional Euclidean space::

        sage: E = EuclideanSpace(4); E
        4-dimensional Euclidean space E^4
        sage: latex(E)
        \mathbb{E}^{4}

    ``E`` is both a real smooth manifold of dimension `4` and a complete metric
    space::

        sage: E.category()
        Join of Category of smooth manifolds over Real Field with 53 bits of
         precision and Category of connected manifolds over Real Field with
         53 bits of precision and Category of complete metric spaces
        sage: dim(E)
        4

    It is endowed with a default coordinate chart, which is that of
    Cartesian coordinates `(x_1,x_2,x_3,x_4)`::

        sage: E.atlas()
        [Chart (E^4, (x1, x2, x3, x4))]
        sage: E.default_chart()
        Chart (E^4, (x1, x2, x3, x4))
        sage: E.default_chart() is E.cartesian_coordinates()
        True

    ``E`` is also endowed with a default metric tensor, which defines the
    Euclidean scalar product::

        sage: g = E.metric(); g
        Riemannian metric g on the 4-dimensional Euclidean space E^4
        sage: g.display()
        g = dx1⊗dx1 + dx2⊗dx2 + dx3⊗dx3 + dx4⊗dx4
    """
    @staticmethod
    def __classcall_private__(cls, n=None, name=None, latex_name=None,
                              coordinates='Cartesian', symbols=None,
                              metric_name='g', metric_latex_name=None,
                              start_index=1, names=None, unique_tag=None):
        r"""
        Determine the correct class to return based upon the input.

        TESTS::

            sage: E2.<x,y> = EuclideanSpace(); E2
            Euclidean plane E^2
            sage: type(E2)
            <class 'sage.manifolds.differentiable.examples.euclidean.EuclideanPlane_with_category'>

            sage: E3.<r,t,p> = EuclideanSpace(coordinates='spherical'); E3
            Euclidean space E^3
            sage: type(E3)
            <class 'sage.manifolds.differentiable.examples.euclidean.Euclidean3dimSpace_with_category'>
            sage: E3.default_frame()._latex_indices
            (r, {\theta}, {\phi})
        """
        if n is None:
            if names is None:
                raise ValueError("either n or names must be specified")
            n = len(names)

        # Parse names into symbols
        if names is not None and symbols is None:
            if n == 2:
                names = list(names)
                if coordinates == 'polar':
                    if names[1] in ['p', 'ph', 'phi']:
                        names[1] += ':\\phi'
                    elif names[1] in ['t', 'th', 'theta']:
                        names[1] += ':\\theta'
            elif n == 3:
                names = list(names)
                # We add the LaTeX symbols when relevant:
                if coordinates == 'spherical':
                    if names[1] in ['t', 'th', 'theta']:
                        names[1] = names[1] + ':\\theta'
                    if names[2] in ['p', 'ph', 'phi']:
                        names[2] = names[2] + ':\\phi'
                elif coordinates == 'cylindrical':
                    if names[0] in ['rh', 'rho']:
                        names[0] = names[0] + ':\\rho'
                    if names[1] in ['p', 'ph', 'phi']:
                        names[1] = names[1] + ':\\phi'

            symbols = ' '.join(names)

        # Technical bit for UniqueRepresentation
        from sage.misc.prandom import getrandbits
        from time import time
        if unique_tag is None:
            unique_tag = getrandbits(128) * time()

        if n == 2:
            return EuclideanPlane(name=name, latex_name=latex_name,
                                  coordinates=coordinates, symbols=symbols,
                                  metric_name=metric_name,
                                  metric_latex_name=metric_latex_name,
                                  start_index=start_index,
                                  unique_tag=unique_tag)
        if n == 3:
            return Euclidean3dimSpace(name=name, latex_name=latex_name,
                                      coordinates=coordinates, symbols=symbols,
                                      metric_name=metric_name,
                                      metric_latex_name=metric_latex_name,
                                      start_index=start_index,
                                      unique_tag=unique_tag)

        return super().__classcall__(cls,
                                     n, name=name, latex_name=latex_name,
                                     coordinates=coordinates, symbols=symbols,
                                     metric_name=metric_name,
                                     metric_latex_name=metric_latex_name,
                                     start_index=start_index,
                                     unique_tag=unique_tag)

    def __init__(self, n, name=None, latex_name=None,
                 coordinates='Cartesian', symbols=None, metric_name='g',
                 metric_latex_name=None, start_index=1, base_manifold=None,
                 category=None, init_coord_methods=None,
                 unique_tag=None):
        r"""
        Construct a Euclidean space.

        INPUT:

        This class also takes the following input:

        - ``base_manifold`` -- (default: ``None``) if not ``None``, must be
          a Euclidean space; the created object is then an open subset
          of ``base_manifold``
        - ``category`` -- (default: ``None``) to specify the category;
          if ``None``,
          ``Manifolds(RR).Smooth() & MetricSpaces().Complete()`` is assumed
        - ``init_coord_methods`` -- (default: ``None``) dictionary of
          methods to initialize the various type of coordinates, with each
          key being a string describing the type of coordinates; to be
          used by derived classes only
        - ``unique_tag`` -- (default: ``None``) tag used to force the
          construction of a new object when all the other arguments have
          been used previously (without ``unique_tag``, the
          :class:`~sage.structure.unique_representation.UniqueRepresentation`
          behavior inherited from
          :class:`~sage.manifolds.differentiable.pseudo_riemannian.PseudoRiemannianManifold`
          would return the previously constructed object corresponding
          to these arguments)

        TESTS::

            sage: E = EuclideanSpace(4); E
            4-dimensional Euclidean space E^4
            sage: E.metric()
            Riemannian metric g on the 4-dimensional Euclidean space E^4
            sage: TestSuite(E).run()
        """
        if name is None:
            name = 'E^{}'.format(n)
            if latex_name is None:
                latex_name = r'\mathbb{E}^{' + str(n) + '}'
        if category is None:
            category = Manifolds(RR).Smooth().Connected() & MetricSpaces().Complete()
            # NB: RR is a proxy for the field of real numbers, until
            #     Issue #24456 is ready
        PseudoRiemannianManifold.__init__(self, n, name, metric_name=metric_name,
                                          signature=n, base_manifold=base_manifold,
                                          latex_name=latex_name,
                                          metric_latex_name=metric_latex_name,
                                          start_index=start_index,
                                          category=category)
        if symbols is None:
            if n == 1:
                if coordinates == 'Cartesian':
                    symbols = 'x'
                else:
                    raise TypeError("unknown coordinate type")
            elif n > 3:
                if coordinates == 'Cartesian':
                    symbols = ''
                    for i in self.irange():
                        symbols += "x{}".format(i) + r":x_{" + str(i) + r"} "
                    symbols = symbols[:-1]
                else:
                    raise TypeError("unknown coordinate type")
            else:
                raise NotImplementedError("dimension not implemented yet")
        self._cartesian_chart = None  # to be constructed later if necessary
        if init_coord_methods is None:
            self._init_coordinates = {'Cartesian': self._init_cartesian}
        else:
            self._init_coordinates = init_coord_methods
        self._init_coordinates[coordinates](symbols)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: E = EuclideanSpace(4)
            sage: E._repr_()
            '4-dimensional Euclidean space E^4'
            sage: E  # indirect doctest
            4-dimensional Euclidean space E^4
        """
        return "{}-dimensional Euclidean space {}".format(self._dim, self._name)

    def _first_ngens(self, n):
        r"""
        Return the list of coordinates of the default chart.

        This is useful only for the use of Sage preparser::

            sage: preparse("E.<x,y,z> = EuclideanSpace()")
            "E = EuclideanSpace(names=('x', 'y', 'z',));
             (x, y, z,) = E._first_ngens(3)"

        TESTS::

            sage: E = EuclideanSpace(2)
            sage: E._first_ngens(2)
            (x, y)
            sage: E.<u,v> = EuclideanSpace()
            sage: E._first_ngens(2)
            (u, v)
        """
        return self._def_chart[:]

    def _init_cartesian(self, symbols):
        r"""
        Construct the chart of Cartesian coordinates and initialize the
        components of the metric tensor in it.

        TESTS::

            sage: E = EuclideanSpace(2)
            sage: E._init_cartesian('x y')
        """
        chart = self.chart(coordinates=symbols)
        self._cartesian_chart = chart
        frame = chart.frame()
        # Renaming (∂/∂x, ∂/∂y, ...) to (e_x, e_y, ...):
        coords = chart[:]
        frame.set_name('e',
                       indices=tuple(str(x) for x in coords),
                       latex_indices=tuple(latex(x) for x in coords))
        g = self.metric()
        gc = g.add_comp(frame)
        for i in self.irange():
            gc[[i, i]] = self.one_scalar_field()
        nabla = g.connection(init_coef=False)  # False to avoid any computation
        nabla.add_coef(frame)  # initialize a zero set of coefficients

    def cartesian_coordinates(self, symbols=None, names=None):
        r"""
        Return the chart of Cartesian coordinates, possibly creating it if it
        does not already exist.

        INPUT:

        - ``symbols`` -- (default: ``None``) string defining the coordinate
          text symbols and LaTeX symbols, with the same conventions as
          the argument ``coordinates`` in
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`; this is
          used only if the Cartesian chart has not been already defined; if
          ``None`` the symbols are generated as `(x_1,\ldots,x_n)`.
        - ``names`` -- (default: ``None``) unused argument, except if
          ``symbols`` is not provided; it must be a tuple containing
          the coordinate symbols (this is guaranteed if the shortcut operator
          ``<,>`` is used)

        OUTPUT:

        - the chart of Cartesian coordinates, as an instance of
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`

        EXAMPLES::

            sage: E = EuclideanSpace(4)
            sage: X = E.cartesian_coordinates(); X
            Chart (E^4, (x1, x2, x3, x4))
            sage: X.coord_range()
            x1: (-oo, +oo); x2: (-oo, +oo); x3: (-oo, +oo); x4: (-oo, +oo)
            sage: X[2]
            x2
            sage: X[:]
            (x1, x2, x3, x4)
            sage: latex(X[:])
            \left({x_{1}}, {x_{2}}, {x_{3}}, {x_{4}}\right)
        """
        if self._cartesian_chart is None:
            if symbols is None:
                symbols = ''
                if names is None:
                    for i in self.irange():
                        symbols += "x{}".format(i) + r":x_{" + str(i) + r"} "
                else:
                    for x in names:
                        symbols += x + ' '
                symbols = symbols[:-1]
            self._init_cartesian(symbols)
        return self._cartesian_chart

    def cartesian_frame(self):
        r"""
        Return the orthonormal vector frame associated with Cartesian
        coordinates.

        OUTPUT: :class:`~sage.manifolds.differentiable.vectorframe.CoordFrame`

        EXAMPLES::

            sage: E = EuclideanSpace(2)
            sage: E.cartesian_frame()
            Coordinate frame (E^2, (e_x,e_y))
            sage: E.cartesian_frame()[1]
            Vector field e_x on the Euclidean plane E^2
            sage: E.cartesian_frame()[:]
            (Vector field e_x on the Euclidean plane E^2,
             Vector field e_y on the Euclidean plane E^2)

        For Cartesian coordinates, the orthonormal frame coincides with the
        coordinate frame::

            sage: E.cartesian_frame() is E.cartesian_coordinates().frame()
            True
        """
        if self._cartesian_chart is None:
            self.cartesian_coordinates()  # creates the Cartesian chart
        # Since the coordinate frame of Cartesian coordinates is orthonormal,
        # we simply return this frame:
        return self._cartesian_chart.frame()

    def dist(self, p, q):
        r"""
        Euclidean distance between two points.

        INPUT:

        - ``p`` -- an element of ``self``
        - ``q`` -- an element of ``self``

        OUTPUT:

        - the Euclidean distance `d(p, q)`

        EXAMPLES::

            sage: E.<x,y> = EuclideanSpace()
            sage: p = E((1,0))
            sage: q = E((0,2))
            sage: E.dist(p, q)
            sqrt(5)
            sage: p.dist(q)  # indirect doctest
            sqrt(5)
        """
        chart = self.cartesian_coordinates()
        coords_p = chart(p)
        coords_q = chart(q)
        d2 = 0
        for xp, xq in zip(coords_p, coords_q):
            dx = xp - xq
            d2 += dx*dx
        return sqrt(d2)

    def sphere(self, radius=1, center=None, name=None, latex_name=None,
               coordinates='spherical', names=None):
        r"""
        Return an `(n-1)`-sphere smoothly embedded in ``self``.

        INPUT:

        - ``radius`` -- (default: ``1``) the radius greater than 1 of the sphere
        - ``center`` -- (default: ``None``) point on ``self`` representing the
          barycenter of the sphere
        - ``name`` -- (default: ``None``) string; name (symbol) given to the
          sphere; if ``None``, the name will be generated according to the input
        - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote
          the sphere; if ``None``, the symbol will be generated according to
          the input
        - ``coordinates`` -- (default: ``'spherical'``) string describing the
          type of coordinates to be initialized at the sphere's creation;
          allowed values are

          - ``'spherical'`` spherical coordinates (see
            :meth:`~sage.manifolds.differentiable.examples.sphere.Sphere.spherical_coordinates`))
          - ``'stereographic'`` stereographic coordinates given by the
            stereographic projection (see
            :meth:`~sage.manifolds.differentiable.examples.sphere.Sphere.stereographic_coordinates`)

        - ``names`` -- (default: ``None``) must be a tuple containing
          the coordinate symbols (this guarantees the shortcut operator
          ``<,>`` to function); if ``None``, the usual conventions are used (see
          examples in
          :class:`~sage.manifolds.differentiable.examples.sphere.Sphere`
          for details)

        EXAMPLES:

        Define a 2-sphere with radius 2 centered at `(1,2,3)` in Cartesian
        coordinates::

            sage: E3 = EuclideanSpace(3)
            sage: c = E3.point((1,2,3), name='c'); c
            Point c on the Euclidean space E^3
            sage: S2_2 = E3.sphere(radius=2, center=c); S2_2
            2-sphere S^2_2(c) of radius 2 smoothly embedded in the Euclidean
             space E^3 centered at the Point c

        The ambient space is precisely our previously defined Euclidean space::

            sage: S2_2.ambient() is E3
            True

        The embedding into Euclidean space::

            sage: S2_2.embedding().display()
            iota: S^2_2(c) → E^3
            on A: (theta, phi) ↦ (x, y, z) = (2*cos(phi)*sin(theta) + 1,
                                                 2*sin(phi)*sin(theta) + 2,
                                                 2*cos(theta) + 3)

        See :class:`~sage.manifolds.differentiable.examples.sphere.Sphere`
        for more examples.
        """
        n = self._dim
        if n == 1:
            raise ValueError('Euclidean space must have dimension of at least 2')
        from .sphere import Sphere
        return Sphere(n-1, radius=radius, ambient_space=self,
                      center=center, name=name, latex_name=latex_name,
                      coordinates=coordinates, names=names)

###############################################################################

class EuclideanPlane(EuclideanSpace):
    r"""
    Euclidean plane.

    A *Euclidean plane* is an affine space `E`, whose associated vector space
    is a 2-dimensional vector space over `\RR` and is equipped with a
    positive definite symmetric bilinear form, called the *scalar product* or
    *dot product*.

    The class :class:`EuclideanPlane` inherits from
    :class:`~sage.manifolds.differentiable.pseudo_riemannian.PseudoRiemannianManifold`
    (via :class:`EuclideanSpace`) since a Euclidean plane can be viewed
    as a Riemannian manifold that is diffeomorphic to `\RR^2` and that has a
    flat metric `g`. The Euclidean scalar product is the one defined by the
    Riemannian metric `g`.

    INPUT:

    - ``name`` -- (default: ``None``) string; name (symbol) given to the
      Euclidean plane; if ``None``, the name will be set to ``'E^2'``
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
      Euclidean plane; if ``None``, it is set to ``'\mathbb{E}^{2}'`` if
      ``name`` is  ``None`` and to ``name`` otherwise
    - ``coordinates`` -- (default: ``'Cartesian'``) string describing the type
      of coordinates to be initialized at the Euclidean plane creation;
      allowed values are ``'Cartesian'`` (see :meth:`cartesian_coordinates`)
      and ``'polar'`` (see :meth:`polar_coordinates`)
    - ``symbols`` -- (default: ``None``) string defining the coordinate text
      symbols and LaTeX symbols, with the same conventions as the argument
      ``coordinates`` in
      :class:`~sage.manifolds.differentiable.chart.RealDiffChart`, namely
      ``symbols`` is a string of coordinate fields separated by a blank space,
      where each field contains the coordinate's text symbol and possibly the
      coordinate's LaTeX symbol (when the latter is different from the text
      symbol), both symbols being separated by a colon (``:``); if ``None``,
      the symbols will be automatically generated according to the value of
      ``coordinates``
    - ``metric_name`` -- (default: ``'g'``) string; name (symbol) given to the
      Euclidean metric tensor
    - ``metric_latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the Euclidean metric tensor; if none is provided, it is set to
      ``metric_name``
    - ``start_index`` -- (default: 1) integer; lower value of the range of
      indices used for "indexed objects" in the Euclidean plane, e.g.
      coordinates of a chart
    - ``base_manifold`` -- (default: ``None``) if not ``None``, must be an
      Euclidean plane; the created object is then an open subset of ``base_manifold``
    - ``category`` -- (default: ``None``) to specify the category; if ``None``,
      ``Manifolds(RR).Smooth() & MetricSpaces().Complete()`` is assumed
    - ``names`` -- (default: ``None``) unused argument, except if
      ``symbols`` is not provided; it must then be a tuple containing
      the coordinate symbols (this is guaranteed if the shortcut operator
      ``<,>`` is used)
    - ``init_coord_methods`` -- (default: ``None``) dictionary of methods
      to initialize the various type of coordinates, with each key being a
      string describing the type of coordinates; to be used by derived classes
      only
    - ``unique_tag`` -- (default: ``None``) tag used to force the construction
      of a new object when all the other arguments have been used previously
      (without ``unique_tag``, the
      :class:`~sage.structure.unique_representation.UniqueRepresentation`
      behavior inherited from
      :class:`~sage.manifolds.differentiable.pseudo_riemannian.PseudoRiemannianManifold`
      would return the previously constructed object corresponding to these
      arguments)

    EXAMPLES:

    One creates a Euclidean plane ``E`` with::

        sage: E.<x,y> = EuclideanSpace(); E
        Euclidean plane E^2

    ``E`` is both a real smooth manifold of dimension `2` and a complete metric
    space::

        sage: E.category()
        Join of Category of smooth manifolds over Real Field with 53 bits of
         precision and Category of connected manifolds over Real Field with
         53 bits of precision and Category of complete metric spaces
        sage: dim(E)
        2

    It is endowed with a default coordinate chart, which is that
    of Cartesian coordinates `(x,y)`::

        sage: E.atlas()
        [Chart (E^2, (x, y))]
        sage: E.default_chart()
        Chart (E^2, (x, y))
        sage: cartesian = E.cartesian_coordinates()
        sage: cartesian is E.default_chart()
        True

    A point of ``E``::

        sage: p = E((3,-2)); p
        Point on the Euclidean plane E^2
        sage: cartesian(p)
        (3, -2)
        sage: p in E
        True
        sage: p.parent() is E
        True

    ``E`` is endowed with a default metric tensor, which defines the
    Euclidean scalar product::

        sage: g = E.metric(); g
        Riemannian metric g on the Euclidean plane E^2
        sage: g.display()
        g = dx⊗dx + dy⊗dy

    Curvilinear coordinates can be introduced on ``E``: see
    :meth:`polar_coordinates`.

    .. SEEALSO::

        :ref:`EuclideanSpace_example1`
    """
    def __init__(self, name=None, latex_name=None, coordinates='Cartesian',
                 symbols=None, metric_name='g', metric_latex_name=None,
                 start_index=1, base_manifold=None, category=None, unique_tag=None):
        r"""
        Construct a Euclidean plane.

        TESTS::

            sage: E = EuclideanSpace(2); E
            Euclidean plane E^2
            sage: E.metric()
            Riemannian metric g on the Euclidean plane E^2
            sage: TestSuite(E).run()
        """
        if coordinates not in ['Cartesian', 'polar']:
            raise TypeError("unknown coordinate type")
        if symbols is None:
            if coordinates == 'Cartesian':
                symbols = 'x y'
            elif coordinates == 'polar':
                symbols = 'r ph:\\phi'
        self._polar_chart = None  # to be constructed later if necessary
        self._polar_frame = None  # orthonormal frame associated to polar coord
        init_coord_methods = {'Cartesian': self._init_cartesian,
                              'polar': self._init_polar}
        EuclideanSpace.__init__(self, 2, name=name,
                                latex_name=latex_name,
                                coordinates=coordinates,
                                symbols=symbols,
                                metric_name=metric_name,
                                metric_latex_name=metric_latex_name,
                                start_index=start_index,
                                base_manifold=base_manifold, category=category,
                                init_coord_methods=init_coord_methods)
        if coordinates == 'polar':
            # The default frame is the polar coordinate frame; we change it
            # to the orthonormal polar frame
            self.set_default_frame(self.polar_frame())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: E = EuclideanSpace(2)
            sage: E._repr_()
            'Euclidean plane E^2'
            sage: E  # indirect doctest
            Euclidean plane E^2
            sage: E = EuclideanSpace(2, name='E')
            sage: E._repr_()
            'Euclidean plane E'
        """
        return "Euclidean plane {}".format(self._name)

    def _init_polar(self, symbols):
        r"""
        Construct the chart of polar coordinates and initialize the
        components of the metric tensor in it.

        TESTS::

            sage: E = EuclideanSpace(2)
            sage: E.atlas()
            [Chart (E^2, (x, y))]
            sage: E._init_polar(r"R Phi:\Phi")
            sage: E.atlas()
            [Chart (E^2, (x, y)), Chart (E^2, (R, Phi))]
        """
        coords = symbols.split()  # list of strings, one per coordinate
        # Adding the coordinate ranges:
        coordinates = (coords[0] + ':(0,+oo) ' + coords[1]
                       + ':(0,2*pi):periodic')
        chart = self.chart(coordinates=coordinates)
        self._polar_chart = chart
        frame = chart.frame()
        # Initialization of the metric components and the associated
        # Christoffel symbols
        g = self.metric()
        gc = g.add_comp(frame)
        i1 = self._sindex
        i2 = i1 + 1
        r, ph = chart[:]
        gc[i1, i1, chart] = 1
        gc[i2, i2, chart] = r**2
        # Orthonormal frame associated with polar coordinates:
        to_orthonormal = self.automorphism_field()
        to_orthonormal[frame, i1, i1, chart] = 1
        to_orthonormal[frame, i2, i2, chart] = 1 / r
        oframe = frame.new_frame(to_orthonormal, 'e',
                                 indices=(str(r), str(ph)),
                                 latex_indices=(latex(r), latex(ph)))
        self._polar_frame = oframe
        g.comp(oframe)

    def _transition_polar_cartesian(self):
        r"""
        Transitions between polar and Cartesian coordinates.

        TESTS::

            sage: E = EuclideanSpace(2)
            sage: E._init_polar(r"r ph:\phi")
            sage: E.atlas()
            [Chart (E^2, (x, y)), Chart (E^2, (r, ph))]
            sage: E.coord_changes() # no transition map has been defined yet
            {}
            sage: E._transition_polar_cartesian()
            sage: len(E.coord_changes())  # 2 transition maps have been created
            2

        Tests of the change-of-frame formulas::

            sage: polar = E.polar_coordinates()
            sage: polar_f = E.polar_frame()
            sage: cart_f = E.cartesian_frame()
            sage: E.change_of_frame(cart_f, polar_f)[:, polar]
            [ cos(ph) -sin(ph)]
            [ sin(ph)  cos(ph)]
            sage: E.change_of_frame(polar_f, cart_f)[:, polar]
            [ cos(ph)  sin(ph)]
            [-sin(ph)  cos(ph)]
        """
        # Transition maps polar chart <-> Cartesian chart
        chart_cart = self._cartesian_chart
        chart_pol = self._polar_chart
        x, y = chart_cart[:]
        r, ph = chart_pol[:]
        pol_to_cart = chart_pol.transition_map(chart_cart,
                                               [r*cos(ph), r*sin(ph)])
        pol_to_cart.set_inverse(sqrt(x**2+y**2), atan2(y,x), check=False)
        # Automorphism Cartesian frame → orthonormal polar frame:
        oframe = self._polar_frame
        cframe = chart_cart.frame()
        sframe = chart_pol.frame()
        chg = self._frame_changes
        cframe_to_oframe = chg[(sframe, oframe)] * chg[(cframe, sframe)]
        # cframe_to_oframe has been computed only in sframe;
        # its components in oframe are computed by
        cmp_of = cframe_to_oframe.comp(oframe)
        # while its components in cframe are obtained by identifying the
        # matrices in cframe and oframe:
        cmp_cf = cframe_to_oframe.add_comp(cframe)
        for i in self.irange():
            for j in self.irange():
                cmp_cf[[i,j]] = cmp_of[[i,j]]
        # Automorphism orthonormal polar frame → Cartesian frame:
        oframe_to_cframe = chg[(sframe, cframe)] * chg[(oframe, sframe)]
        # oframe_to_cframe has been computed only in sframe;
        # its components in oframe are computed by
        cmp_of = oframe_to_cframe.comp(oframe)
        # while its components in cframe are obtained by identifying the
        # matrices in cframe and oframe:
        cmp_cf = oframe_to_cframe.add_comp(cframe)
        for i in self.irange():
            for j in self.irange():
                cmp_cf[[i,j]] = cmp_of[[i,j]]
        # Storage of the results:
        chg[(cframe, oframe)] = cframe_to_oframe
        chg[(oframe, cframe)] = oframe_to_cframe
        vmodule = self.vector_field_module()
        vmodule.set_change_of_basis(cframe, oframe, cframe_to_oframe,
                                    compute_inverse=False)
        vmodule.set_change_of_basis(oframe, cframe, oframe_to_cframe,
                                    compute_inverse=False)

    def cartesian_coordinates(self, symbols=None, names=None):
        r"""
        Return the chart of Cartesian coordinates, possibly creating it if it
        does not already exist.

        INPUT:

        - ``symbols`` -- (default: ``None``) string defining the coordinate
          text symbols and LaTeX symbols, with the same conventions as the
          argument ``coordinates`` in
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`; this is
          used only if the Cartesian chart has not been already defined; if
          ``None`` the symbols are generated as `(x,y)`.
        - ``names`` -- (default: ``None``) unused argument, except if
          ``symbols`` is not provided; it must be a tuple containing
          the coordinate symbols (this is guaranteed if the shortcut operator
          ``<,>`` is used)

        OUTPUT:

        - the chart of Cartesian coordinates, as an instance of
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`

        EXAMPLES::

            sage: E = EuclideanSpace(2)
            sage: E.cartesian_coordinates()
            Chart (E^2, (x, y))
            sage: E.cartesian_coordinates().coord_range()
            x: (-oo, +oo); y: (-oo, +oo)

        An example where the Cartesian coordinates have not been previously
        created::

            sage: E = EuclideanSpace(2, coordinates='polar')
            sage: E.atlas()  # only polar coordinates have been initialized
            [Chart (E^2, (r, ph))]
            sage: E.cartesian_coordinates(symbols='X Y')
            Chart (E^2, (X, Y))
            sage: E.atlas()  # the Cartesian chart has been added to the atlas
            [Chart (E^2, (r, ph)), Chart (E^2, (X, Y))]

        Note that if the Cartesian coordinates have been already initialized,
        the argument ``symbols`` has no effect::

            sage: E.cartesian_coordinates(symbols='x y')
            Chart (E^2, (X, Y))

        The coordinate variables are returned by the square bracket operator::

            sage: E.cartesian_coordinates()[1]
            X
            sage: E.cartesian_coordinates()[2]
            Y
            sage: E.cartesian_coordinates()[:]
            (X, Y)

        It is also possible to use the operator ``<,>`` to set symbolic
        variable containing the coordinates::

            sage: E = EuclideanSpace(2, coordinates='polar')
            sage: cartesian.<u,v> = E.cartesian_coordinates()
            sage: cartesian
            Chart (E^2, (u, v))
            sage: u,v
            (u, v)

        The command ``cartesian.<u,v> = E.cartesian_coordinates()`` is
        actually a shortcut for::

            sage: cartesian = E.cartesian_coordinates(symbols='u v')
            sage: u, v = cartesian[:]
        """
        if self._cartesian_chart is None:
            if symbols is None:
                if names is None:
                    symbols = 'x y'
                else:
                    symbols = names[0] + ' ' + names[1]
            self._init_cartesian(symbols)
            if self._polar_chart:
                self._transition_polar_cartesian()
        return self._cartesian_chart

    def polar_coordinates(self, symbols=None, names=None):
        r"""
        Return the chart of polar coordinates, possibly creating it if it
        does not already exist.

        INPUT:

        - ``symbols`` -- (default: ``None``) string defining the coordinate
          text symbols and LaTeX symbols, with the same conventions as the
          argument ``coordinates`` in
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`; this is
          used only if the polar chart has not been already defined; if
          ``None`` the symbols are generated as `(r,\phi)`.
        - ``names`` -- (default: ``None``) unused argument, except if
          ``symbols`` is not provided; it must be a tuple containing
          the coordinate symbols (this is guaranteed if the shortcut operator
          ``<,>`` is used)

        OUTPUT:

        - the chart of polar coordinates, as an instance of
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`

        EXAMPLES::

            sage: E = EuclideanSpace(2)
            sage: E.polar_coordinates()
            Chart (E^2, (r, ph))
            sage: latex(_)
            \left(\mathbb{E}^{2},(r, {\phi})\right)
            sage: E.polar_coordinates().coord_range()
            r: (0, +oo); ph: [0, 2*pi] (periodic)

        The relation to Cartesian coordinates is::

            sage: E.coord_change(E.polar_coordinates(),
            ....:                E.cartesian_coordinates()).display()
            x = r*cos(ph)
            y = r*sin(ph)
            sage: E.coord_change(E.cartesian_coordinates(),
            ....:                E.polar_coordinates()).display()
            r = sqrt(x^2 + y^2)
            ph = arctan2(y, x)

        The coordinate variables are returned by the square bracket operator::

            sage: E.polar_coordinates()[1]
            r
            sage: E.polar_coordinates()[2]
            ph
            sage: E.polar_coordinates()[:]
            (r, ph)

        They can also be obtained via the operator ``<,>``::

            sage: polar.<r,ph> = E.polar_coordinates(); polar
            Chart (E^2, (r, ph))
            sage: r, ph
            (r, ph)

        Actually, ``polar.<r,ph> = E.polar_coordinates()`` is a shortcut for::

            sage: polar = E.polar_coordinates()
            sage: r, ph = polar[:]

        The coordinate symbols can be customized::

            sage: E = EuclideanSpace(2)
            sage: E.polar_coordinates(symbols=r"r th:\theta")
            Chart (E^2, (r, th))
            sage: latex(E.polar_coordinates())
            \left(\mathbb{E}^{2},(r, {\theta})\right)

        Note that if the polar coordinates have been already initialized, the
        argument ``symbols`` has no effect::

            sage: E.polar_coordinates(symbols=r"R Th:\Theta")
            Chart (E^2, (r, th))
        """
        if self._polar_chart is None:
            if symbols is None:
                if names is None:
                    symbols = 'r ph:\\phi'
                else:
                    symbols = names[0] + ' ' + names[1]
                    if names[1] in ['p', 'ph', 'phi']:
                        symbols += ':\\phi'
                    elif names[1] in ['t', 'th', 'theta']:
                        symbols += ':\\theta'
            self._init_polar(symbols)
            if self._cartesian_chart:
                self._transition_polar_cartesian()
        return self._polar_chart

    def polar_frame(self):
        r"""
        Return the orthonormal vector frame associated with polar
        coordinates.

        OUTPUT:

        - instance of
          :class:`~sage.manifolds.differentiable.vectorframe.VectorFrame`

        EXAMPLES::

            sage: E = EuclideanSpace(2)
            sage: E.polar_frame()
            Vector frame (E^2, (e_r,e_ph))
            sage: E.polar_frame()[1]
            Vector field e_r on the Euclidean plane E^2
            sage: E.polar_frame()[:]
            (Vector field e_r on the Euclidean plane E^2,
             Vector field e_ph on the Euclidean plane E^2)

        The orthonormal polar frame expressed in terms of the Cartesian one::

            sage: for e in E.polar_frame():
            ....:     e.display(E.cartesian_frame(), E.polar_coordinates())
            e_r = cos(ph) e_x + sin(ph) e_y
            e_ph = -sin(ph) e_x + cos(ph) e_y

        The orthonormal frame `(e_r, e_\phi)` expressed in terms of the
        coordinate frame
        `\left(\frac{\partial}{\partial r},
        \frac{\partial}{\partial\phi}\right)`::

            sage: for e in E.polar_frame():
            ....:     e.display(E.polar_coordinates())
            e_r = ∂/∂r
            e_ph = 1/r ∂/∂ph
        """
        if self._polar_frame is None:
            # create the polar chart and the associated orthonormal frame
            self.polar_coordinates()
        return self._polar_frame


###############################################################################

class Euclidean3dimSpace(EuclideanSpace):
    r"""
    3-dimensional Euclidean space.

    A *3-dimensional Euclidean space* is an affine space `E`, whose associated
    vector space is a 3-dimensional vector space over `\RR` and is equipped
    with a positive definite symmetric bilinear form, called the *scalar
    product* or *dot product*.

    The class :class:`Euclidean3dimSpace` inherits from
    :class:`~sage.manifolds.differentiable.pseudo_riemannian.PseudoRiemannianManifold`
    (via :class:`EuclideanSpace`) since a 3-dimensional Euclidean space
    can be viewed as a Riemannian manifold that is diffeomorphic to `\RR^3` and
    that has a flat metric `g`. The Euclidean scalar product is the one defined
    by the Riemannian metric `g`.

    INPUT:

    - ``name`` -- (default: ``None``) string; name (symbol) given to the
      Euclidean 3-space; if ``None``, the name will be set to ``'E^3'``
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
      Euclidean 3-space; if ``None``, it is set to ``'\mathbb{E}^{3}'`` if
      ``name`` is  ``None`` and to ``name`` otherwise
    - ``coordinates`` -- (default: ``'Cartesian'``) string describing the type
      of coordinates to be initialized at the Euclidean 3-space creation;
      allowed values are ``'Cartesian'`` (see :meth:`cartesian_coordinates`),
      ``'spherical'`` (see :meth:`spherical_coordinates`) and ``'cylindrical'``
      (see :meth:`cylindrical_coordinates`)
    - ``symbols`` -- (default: ``None``) string defining the coordinate text
      symbols and LaTeX symbols, with the same conventions as the argument
      ``coordinates`` in
      :class:`~sage.manifolds.differentiable.chart.RealDiffChart`, namely
      ``symbols`` is a string of coordinate fields separated by a blank space,
      where each field contains the coordinate's text symbol and possibly the
      coordinate's LaTeX symbol (when the latter is different from the text
      symbol), both symbols being separated by a colon (``:``); if ``None``,
      the symbols will be automatically generated according to the value of
      ``coordinates``
    - ``metric_name`` -- (default: ``'g'``) string; name (symbol) given to the
      Euclidean metric tensor
    - ``metric_latex_name`` -- (default: ``None``) string; LaTeX symbol to
      denote the Euclidean metric tensor; if none is provided, it is set to
      ``metric_name``
    - ``start_index`` -- (default: 1) integer; lower value of the range of
      indices used for "indexed objects" in the Euclidean 3-space, e.g.
      coordinates of a chart
    - ``base_manifold`` -- (default: ``None``) if not ``None``, must be a
      Euclidean 3-space; the created object is then an open subset of
      ``base_manifold``
    - ``category`` -- (default: ``None``) to specify the category; if ``None``,
      ``Manifolds(RR).Smooth() & MetricSpaces().Complete()`` is assumed
    - ``names`` -- (default: ``None``) unused argument, except if
      ``symbols`` is not provided; it must then be a tuple containing
      the coordinate symbols (this is guaranteed if the shortcut operator
      ``<,>`` is used)
    - ``init_coord_methods`` -- (default: ``None``) dictionary of methods
      to initialize the various type of coordinates, with each key being a
      string describing the type of coordinates; to be used by derived classes
      only
    - ``unique_tag`` -- (default: ``None``) tag used to force the construction
      of a new object when all the other arguments have been used previously
      (without ``unique_tag``, the
      :class:`~sage.structure.unique_representation.UniqueRepresentation`
      behavior inherited from
      :class:`~sage.manifolds.differentiable.pseudo_riemannian.PseudoRiemannianManifold`
      would return the previously constructed object corresponding to these
      arguments)

    EXAMPLES:

    A 3-dimensional Euclidean space::

        sage: E = EuclideanSpace(3); E
        Euclidean space E^3
        sage: latex(E)
        \mathbb{E}^{3}

    ``E`` belongs to the class :class:`Euclidean3dimSpace` (actually to
    a dynamically generated subclass of it via SageMath's category framework)::

        sage: type(E)
        <class 'sage.manifolds.differentiable.examples.euclidean.Euclidean3dimSpace_with_category'>

    ``E`` is both a real smooth manifold of dimension `3` and a complete metric
    space::

        sage: E.category()
        Join of Category of smooth manifolds over Real Field with 53 bits of
         precision and Category of connected manifolds over Real Field with
         53 bits of precision and Category of complete metric spaces
        sage: dim(E)
        3

    It is endowed with a default coordinate chart, which is that of Cartesian
    coordinates `(x,y,z)`::

        sage: E.atlas()
        [Chart (E^3, (x, y, z))]
        sage: E.default_chart()
        Chart (E^3, (x, y, z))
        sage: cartesian = E.cartesian_coordinates()
        sage: cartesian is E.default_chart()
        True

    A point of ``E``::

        sage: p = E((3,-2,1)); p
        Point on the Euclidean space E^3
        sage: cartesian(p)
        (3, -2, 1)
        sage: p in E
        True
        sage: p.parent() is E
        True

    ``E`` is endowed with a default metric tensor, which defines the
    Euclidean scalar product::

        sage: g = E.metric(); g
        Riemannian metric g on the Euclidean space E^3
        sage: g.display()
        g = dx⊗dx + dy⊗dy + dz⊗dz

    Curvilinear coordinates can be introduced on ``E``: see
    :meth:`spherical_coordinates` and :meth:`cylindrical_coordinates`.

    .. SEEALSO::

        :ref:`EuclideanSpace_example2`
    """
    def __init__(self, name=None, latex_name=None, coordinates='Cartesian',
                 symbols=None, metric_name='g', metric_latex_name=None,
                 start_index=1, base_manifold=None, category=None, unique_tag=None):
        r"""
        Construct a Euclidean 3-space.

        TESTS::

            sage: E = EuclideanSpace(3); E
            Euclidean space E^3
            sage: E.metric()
            Riemannian metric g on the Euclidean space E^3
            sage: TestSuite(E).run()
        """
        if coordinates not in ['Cartesian', 'spherical', 'cylindrical']:
            raise TypeError("unknown coordinate type")
        if symbols is None:
            if coordinates == 'Cartesian':
                symbols = 'x y z'
            elif coordinates == 'spherical':
                symbols = 'r th:\\theta ph:\\phi'
            elif coordinates == 'cylindrical':
                symbols = 'r ph:\\phi z'
        self._spherical_chart = None    # to be constructed later if necessary
        self._spherical_frame = None    # orthonormal frame
        self._cylindrical_chart = None
        self._cylindrical_frame = None  # orthonormal frame
        init_coord_methods = {'Cartesian': self._init_cartesian,
                              'spherical': self._init_spherical,
                              'cylindrical': self._init_cylindrical}
        EuclideanSpace.__init__(self, 3, name=name,
                                latex_name=latex_name,
                                coordinates=coordinates,
                                symbols=symbols,
                                metric_name=metric_name,
                                metric_latex_name=metric_latex_name,
                                start_index=start_index,
                                base_manifold=base_manifold, category=category,
                                init_coord_methods=init_coord_methods)
        if coordinates == 'spherical':
            # The default frame is the spherical coordinate frame; we change it
            # to the orthonormal spherical frame
            self.set_default_frame(self.spherical_frame())
        if coordinates == 'cylindrical':
            # The default frame is the cylindrical coordinate frame; we change
            # it to the orthonormal cylindrical frame
            self.set_default_frame(self.cylindrical_frame())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        TESTS::

            sage: E = EuclideanSpace(3)
            sage: E._repr_()
            'Euclidean space E^3'
            sage: E  # indirect doctest
            Euclidean space E^3
            sage: E = EuclideanSpace(3, name='E')
            sage: E._repr_()
            'Euclidean space E'
        """
        return "Euclidean space {}".format(self._name)

    def _init_spherical(self, symbols):
        r"""
        Construct the chart of spherical coordinates and initialize the
        components of the metric tensor in it.

        TESTS::

            sage: E = EuclideanSpace(3)
            sage: E.atlas()
            [Chart (E^3, (x, y, z))]
            sage: E._init_spherical(r"R Th:\Theta Ph:\Phi")
            sage: E.atlas()
            [Chart (E^3, (x, y, z)), Chart (E^3, (R, Th, Ph))]
        """
        coords = symbols.split()  # list of strings, one per coordinate
        # Adding the coordinate ranges:
        coordinates = (coords[0] + ':(0,+oo) ' + coords[1] + ':(0,pi) '
                       + coords[2] + ':(0,2*pi):periodic')
        chart = self.chart(coordinates=coordinates)
        self._spherical_chart = chart
        frame = chart.frame()
        # Initialization of the metric components and the associated
        # Christoffel symbols
        g = self.metric()
        gc = g.add_comp(frame)
        i1 = self._sindex
        i2 = i1 + 1
        i3 = i1 + 2
        r, th, ph = chart[:]
        gc[i1, i1, chart] = 1
        gc[i2, i2, chart] = r**2
        gc[i3, i3, chart] = (r*sin(th))**2
        # Orthonormal frame associated with spherical coordinates:
        to_orthonormal = self.automorphism_field()
        to_orthonormal[frame, i1, i1, chart] = 1
        to_orthonormal[frame, i2, i2, chart] = 1/r
        to_orthonormal[frame, i3, i3, chart] = 1/(r*sin(th))
        oframe = frame.new_frame(to_orthonormal, 'e',
                                 indices=(str(r), str(th), str(ph)),
                                 latex_indices=(latex(r), latex(th), latex(ph)))
        self._spherical_frame = oframe
        g.comp(oframe)

    def _init_cylindrical(self, symbols):
        r"""
        Construct the chart of cylindrical coordinates and initialize the
        components of the metric tensor in it.

        TESTS::

            sage: E = EuclideanSpace(3)
            sage: E.atlas()
            [Chart (E^3, (x, y, z))]
            sage: E._init_cylindrical(r"r ph:\phi z")
            sage: E.atlas()
            [Chart (E^3, (x, y, z)), Chart (E^3, (r, ph, z))]
        """
        coords = symbols.split()  # list of strings, one per coordinate
        # Adding the coordinate ranges:
        coordinates = (coords[0] + ':(0,+oo) ' + coords[1]
                       + ':(0,2*pi):periodic ' + coords[2])
        chart = self.chart(coordinates=coordinates)
        self._cylindrical_chart = chart
        frame = chart.frame()
        # Initialization of the metric components and the associated
        # Christoffel symbols
        g = self.metric()
        gc = g.add_comp(frame)
        i1 = self._sindex
        i2 = i1 + 1
        i3 = i1 + 2
        rh, ph, z = chart[:]
        gc[i1, i1, chart] = 1
        gc[i2, i2, chart] = rh**2
        gc[i3, i3, chart] = 1
        # Orthonormal frame associated with cylindrical coordinates:
        to_orthonormal = self.automorphism_field()
        to_orthonormal[frame, i1, i1, chart] = 1
        to_orthonormal[frame, i2, i2, chart] = 1 / rh
        to_orthonormal[frame, i3, i3, chart] = 1
        oframe = frame.new_frame(to_orthonormal, 'e',
                                 indices=(str(rh), str(ph), str(z)),
                                 latex_indices=(latex(rh), latex(ph), latex(z)))
        self._cylindrical_frame = oframe
        g.comp(oframe)

    def _transition_spherical_cartesian(self):
        r"""
        Transitions between spherical and Cartesian coordinates.

        TESTS::

            sage: E = EuclideanSpace(3)
            sage: E._init_spherical(r"r th:\theta ph:\phi")
            sage: E.atlas()
            [Chart (E^3, (x, y, z)), Chart (E^3, (r, th, ph))]
            sage: E.coord_changes()  # no transition map has been defined yet
            {}
            sage: E._transition_spherical_cartesian()  # long time

        2 transition maps have been created::

            sage: len(E.coord_changes())  # long time
            2

        Tests of the change-of-frame formulas::

            sage: # long time
            sage: spher = E.spherical_coordinates()
            sage: spher_f = E.spherical_frame()
            sage: cart_f = E.cartesian_frame()
            sage: E.change_of_frame(cart_f, spher_f)[:,spher]
            [cos(ph)*sin(th) cos(ph)*cos(th)        -sin(ph)]
            [sin(ph)*sin(th) cos(th)*sin(ph)         cos(ph)]
            [        cos(th)        -sin(th)               0]
            sage: E.change_of_frame(spher_f, cart_f)[:,spher]
            [cos(ph)*sin(th) sin(ph)*sin(th)         cos(th)]
            [cos(ph)*cos(th) cos(th)*sin(ph)        -sin(th)]
            [       -sin(ph)         cos(ph)               0]
        """
        # Transition maps spherical chart <-> Cartesian chart
        chart_cart = self._cartesian_chart
        chart_spher = self._spherical_chart
        x, y, z = chart_cart[:]
        r, th, ph = chart_spher[:]
        spher_to_cart = chart_spher.transition_map(chart_cart,
                             [r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th)])
        spher_to_cart.set_inverse(sqrt(x**2+y**2+z**2),
                                  atan2(sqrt(x**2+y**2),z), atan2(y, x),
                                  check=False)
        # Automorphism Cartesian frame → orthonormal spherical frame:
        oframe = self._spherical_frame
        cframe = chart_cart.frame()
        sframe = chart_spher.frame()
        chg = self._frame_changes
        cframe_to_oframe = chg[(sframe, oframe)] * chg[(cframe, sframe)]
        # cframe_to_oframe has been computed only in sframe;
        # its components in oframe are computed by
        cmp_of = cframe_to_oframe.comp(oframe)
        # while its components in cframe are obtained by identifying the
        # matrices in cframe and oframe:
        cmp_cf = cframe_to_oframe.add_comp(cframe)
        for i in self.irange():
            for j in self.irange():
                cmp_cf[[i,j]] = cmp_of[[i,j]]
        # Automorphism orthonormal spherical frame → Cartesian frame:
        oframe_to_cframe = chg[(sframe, cframe)] * chg[(oframe, sframe)]
        # oframe_to_cframe has been computed only in sframe;
        # its components in oframe are computed by
        cmp_of = oframe_to_cframe.comp(oframe)
        # while its components in cframe are obtained by identifying the
        # matrices in cframe and oframe:
        cmp_cf = oframe_to_cframe.add_comp(cframe)
        for i in self.irange():
            for j in self.irange():
                cmp_cf[[i,j]] = cmp_of[[i,j]]
        # Storage of the results:
        chg[(cframe, oframe)] = cframe_to_oframe
        chg[(oframe, cframe)] = oframe_to_cframe
        vmodule = self.vector_field_module()
        vmodule.set_change_of_basis(cframe, oframe, cframe_to_oframe,
                                    compute_inverse=False)
        vmodule.set_change_of_basis(oframe, cframe, oframe_to_cframe,
                                    compute_inverse=False)

    def _transition_cylindrical_cartesian(self):
        r"""
        Transitions between cylindrical and Cartesian coordinates.

        TESTS::

            sage: E = EuclideanSpace(3)
            sage: E._init_cylindrical(r"r ph:\phi z")
            sage: E.atlas()
            [Chart (E^3, (x, y, z)), Chart (E^3, (r, ph, z))]
            sage: E.coord_changes() # no transition map has been defined yet
            {}
            sage: E._transition_cylindrical_cartesian()  # long time

        2 transition maps have been created::

            sage: len(E.coord_changes())  # long time
            2

        Tests of the change-of-frame formulas::

            sage: # long time
            sage: cylind = E.cylindrical_coordinates()
            sage: cylind_f = E.cylindrical_frame()
            sage: cart_f= E.cartesian_frame()
            sage: E.change_of_frame(cart_f, cylind_f)[:,cylind]
            [ cos(ph) -sin(ph)        0]
            [ sin(ph)  cos(ph)        0]
            [       0        0        1]
            sage: E.change_of_frame(cylind_f, cart_f)[:,cylind]
            [ cos(ph)  sin(ph)        0]
            [-sin(ph)  cos(ph)        0]
            [       0        0        1]
        """
        # Transition maps cylindrical chart <-> Cartesian chart
        chart_cart = self._cartesian_chart
        chart_cylind = self._cylindrical_chart
        x, y, z = chart_cart[:]
        rh, ph, z = chart_cylind[:]
        cylind_to_cart = chart_cylind.transition_map(chart_cart,
                                                     [rh*cos(ph), rh*sin(ph), z])
        cylind_to_cart.set_inverse(sqrt(x**2+y**2), atan2(y, x), z, check=False)
        # Automorphism Cartesian frame → orthonormal cylindrical frame:
        oframe = self._cylindrical_frame
        cframe = chart_cart.frame()
        sframe = chart_cylind.frame()
        chg = self._frame_changes
        cframe_to_oframe = chg[(sframe, oframe)] * chg[(cframe, sframe)]
        # cframe_to_oframe has been computed only in sframe;
        # its components in oframe are computed by
        cmp_of = cframe_to_oframe.comp(oframe)
        # while its components in cframe are obtained by identifying the
        # matrices in cframe and oframe:
        cmp_cf = cframe_to_oframe.add_comp(cframe)
        for i in self.irange():
            for j in self.irange():
                cmp_cf[[i,j]] = cmp_of[[i,j]]
        # Automorphism orthonormal cylindrical frame → Cartesian frame:
        oframe_to_cframe = chg[(sframe, cframe)] * chg[(oframe, sframe)]
        # oframe_to_cframe has been computed only in sframe;
        # its components in oframe are computed by
        cmp_of = oframe_to_cframe.comp(oframe)
        # while its components in cframe are obtained by identifying the
        # matrices in cframe and oframe:
        cmp_cf = oframe_to_cframe.add_comp(cframe)
        for i in self.irange():
            for j in self.irange():
                cmp_cf[[i,j]] = cmp_of[[i,j]]
        # Storage of the results:
        chg[(cframe, oframe)] = cframe_to_oframe
        chg[(oframe, cframe)] = oframe_to_cframe
        vmodule = self.vector_field_module()
        vmodule.set_change_of_basis(cframe, oframe, cframe_to_oframe,
                                    compute_inverse=False)
        vmodule.set_change_of_basis(oframe, cframe, oframe_to_cframe,
                                    compute_inverse=False)

    def _transition_spherical_cylindrical(self):
        r"""
        Transitions between spherical and cylindrical coordinates.

        TESTS::

            sage: E = EuclideanSpace(3, coordinates='cylindrical')
            sage: E._init_spherical(r"r th:\theta ph:\phi")
            sage: E.atlas()
            [Chart (E^3, (r, ph, z)), Chart (E^3, (r, th, ph))]
            sage: E.coord_changes()  # no transition map has been defined yet
            {}
            sage: E._transition_spherical_cylindrical()  # long time

        2 transition maps have been created::

            sage: len(E.coord_changes())  # long time
            2

        Tests of the change-of-frame formulas::

            sage: # long time
            sage: spher = E.spherical_coordinates()
            sage: spher_f = E.spherical_frame()
            sage: cylind_f = E.cylindrical_frame()
            sage: E.change_of_frame(cylind_f, spher_f)[:, spher]
            [ sin(th)  cos(th)        0]
            [       0        0        1]
            [ cos(th) -sin(th)        0]
            sage: E.change_of_frame(spher_f, cylind_f)[:, spher]
            [ sin(th)        0  cos(th)]
            [ cos(th)        0 -sin(th)]
            [       0        1        0]
        """
        # Transition maps spherical chart <-> cylindrical chart
        cylind = self._cylindrical_chart
        spher = self._spherical_chart
        rh, ph, z = cylind[:]
        r, th, ph = spher[:]
        spher_to_cylind = spher.transition_map(cylind,
                                               [r*sin(th), ph, r*cos(th)])
        spher_to_cylind.set_inverse(sqrt(rh**2 + z**2), atan2(rh,z), ph,
                                    check=False)
        # Automorphism orthon. cylindrical frame -> orthon. spherical frame
        cf = cylind.frame() # coordinate cylindrical frame
        sf = spher.frame()  # coordinate spherical frame
        ocf = self._cylindrical_frame # orthonormal cylindrical frame
        osf = self._spherical_frame   # orthonormal spherical frame
        chg = self._frame_changes
        oc_to_os = chg[(sf, osf)] * chg[(cf, sf)] * chg[(ocf, cf)]
        # oc_to_os has been computed only in sf frame; its components in osf
        # frame are computed by:
        cmp_osf = oc_to_os.comp(osf)
        # while its components in ocf frame are obtained by identifying the
        # matrices in ocf frame and osf oframe:
        cmp_ocf = oc_to_os.add_comp(ocf)
        for i in self.irange():
            for j in self.irange():
                cmp_ocf[[i,j]] = cmp_osf[[i,j]]
        # Automorphism orthon. spherical frame -> orthon. cylindrical frame
        os_to_oc = chg[(cf, ocf)] * chg[(sf, cf)] * chg[(osf, sf)]
        # oc_to_os has been computed only in cf frame; its components in ocf
        # frame are computed by
        cmp_ocf = os_to_oc.comp(ocf)
        # while its components in osf frame are obtained by identifying the
        # matrices in osf frame and ocf oframe:
        cmp_osf = os_to_oc.add_comp(osf)
        for i in self.irange():
            for j in self.irange():
                cmp_osf[[i,j]] = cmp_ocf[[i,j]]
        # Storage of the results:
        chg[(ocf, osf)] = oc_to_os
        chg[(osf, ocf)] = os_to_oc
        vmodule = self.vector_field_module()
        vmodule.set_change_of_basis(ocf, osf, oc_to_os, compute_inverse=False)
        vmodule.set_change_of_basis(osf, ocf, os_to_oc, compute_inverse=False)

    def cartesian_coordinates(self, symbols=None, names=None):
        r"""
        Return the chart of Cartesian coordinates, possibly creating it if it
        does not already exist.

        INPUT:

        - ``symbols`` -- (default: ``None``) string defining the coordinate
          text symbols and LaTeX symbols, with the same conventions as the
          argument ``coordinates`` in
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`; this is
          used only if the Cartesian chart has not been already defined; if
          ``None`` the symbols are generated as `(x,y,z)`.
        - ``names`` -- (default: ``None``) unused argument, except if
          ``symbols`` is not provided; it must be a tuple containing
          the coordinate symbols (this is guaranteed if the shortcut operator
          ``<,>`` is used)

        OUTPUT:

        - the chart of Cartesian coordinates, as an instance of
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`

        EXAMPLES::

            sage: E = EuclideanSpace(3)
            sage: E.cartesian_coordinates()
            Chart (E^3, (x, y, z))
            sage: E.cartesian_coordinates().coord_range()
            x: (-oo, +oo); y: (-oo, +oo); z: (-oo, +oo)

        An example where the Cartesian coordinates have not been previously
        created::

            sage: E = EuclideanSpace(3, coordinates='spherical')
            sage: E.atlas()  # only spherical coordinates have been initialized
            [Chart (E^3, (r, th, ph))]
            sage: E.cartesian_coordinates(symbols='X Y Z')
            Chart (E^3, (X, Y, Z))
            sage: E.atlas()  # the Cartesian chart has been added to the atlas
            [Chart (E^3, (r, th, ph)), Chart (E^3, (X, Y, Z))]

        The coordinate variables are returned by the square bracket operator::

            sage: E.cartesian_coordinates()[1]
            X
            sage: E.cartesian_coordinates()[3]
            Z
            sage: E.cartesian_coordinates()[:]
            (X, Y, Z)

        It is also possible to use the operator ``<,>`` to set symbolic
        variable containing the coordinates::

            sage: E = EuclideanSpace(3, coordinates='spherical')
            sage: cartesian.<u,v,w> = E.cartesian_coordinates()
            sage: cartesian
            Chart (E^3, (u, v, w))
            sage: u, v, w
            (u, v, w)

        The command ``cartesian.<u,v,w> = E.cartesian_coordinates()``
        is actually a shortcut for::

            sage: cartesian = E.cartesian_coordinates(symbols='u v w')
            sage: u, v, w = cartesian[:]
        """
        if self._cartesian_chart is None:
            if symbols is None:
                if names is None:
                    symbols = 'x y z'
                else:
                    symbols = names[0] + ' ' + names[1] + ' ' + names[2]
            self._init_cartesian(symbols)
            if self._spherical_chart:
                self._transition_spherical_cartesian()
            if self._cylindrical_chart:
                self._transition_cylindrical_cartesian()
        return self._cartesian_chart

    def spherical_coordinates(self, symbols=None, names=None):
        r"""
        Return the chart of spherical coordinates, possibly creating it if it
        does not already exist.

        INPUT:

        - ``symbols`` -- (default: ``None``) string defining the coordinate
          text symbols and LaTeX symbols, with the same conventions as the
          argument ``coordinates`` in
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`; this is
          used only if the spherical chart has not been already defined; if
          ``None`` the symbols are generated as `(r,\theta,\phi)`.
        - ``names`` -- (default: ``None``) unused argument, except if
          ``symbols`` is not provided; it must be a tuple containing
          the coordinate symbols (this is guaranteed if the shortcut operator
          ``<,>`` is used)

        OUTPUT:

        - the chart of spherical coordinates, as an instance of
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`

        EXAMPLES::

            sage: E = EuclideanSpace(3)
            sage: E.spherical_coordinates()
            Chart (E^3, (r, th, ph))
            sage: latex(_)
            \left(\mathbb{E}^{3},(r, {\theta}, {\phi})\right)
            sage: E.spherical_coordinates().coord_range()
            r: (0, +oo); th: (0, pi); ph: [0, 2*pi] (periodic)

        The relation to Cartesian coordinates is::

            sage: E.coord_change(E.spherical_coordinates(),
            ....:                E.cartesian_coordinates()).display()
            x = r*cos(ph)*sin(th)
            y = r*sin(ph)*sin(th)
            z = r*cos(th)
            sage: E.coord_change(E.cartesian_coordinates(),
            ....:                E.spherical_coordinates()).display()
            r = sqrt(x^2 + y^2 + z^2)
            th = arctan2(sqrt(x^2 + y^2), z)
            ph = arctan2(y, x)

        The coordinate variables are returned by the square bracket operator::

            sage: E.spherical_coordinates()[1]
            r
            sage: E.spherical_coordinates()[3]
            ph
            sage: E.spherical_coordinates()[:]
            (r, th, ph)

        They can also be obtained via the operator ``<,>``::

            sage: spherical.<r,th,ph> = E.spherical_coordinates()
            sage: spherical
            Chart (E^3, (r, th, ph))
            sage: r, th, ph
            (r, th, ph)

        Actually, ``spherical.<r,th,ph> = E.spherical_coordinates()`` is a
        shortcut for::

            sage: spherical = E.spherical_coordinates()
            sage: r, th, ph = spherical[:]

        The coordinate symbols can be customized::

            sage: E = EuclideanSpace(3)
            sage: E.spherical_coordinates(symbols=r"R T:\Theta F:\Phi")
            Chart (E^3, (R, T, F))
            sage: latex(E.spherical_coordinates())
            \left(\mathbb{E}^{3},(R, {\Theta}, {\Phi})\right)

        Note that if the spherical coordinates have been already initialized,
        the argument ``symbols`` has no effect::

            sage: E.spherical_coordinates(symbols=r"r th:\theta ph:\phi")
            Chart (E^3, (R, T, F))
        """
        if self._spherical_chart is None:
            if symbols is None:
                if names is None:
                    symbols = 'r th:\\theta ph:\\phi'
                else:
                    names = list(names)
                    if names[1] in ['t', 'th', 'theta']:
                        names[1] = names[1] + ':\\theta'
                    if names[2] in ['p', 'ph', 'phi']:
                        names[2] = names[2] + ':\\phi'
                    symbols = names[0] + ' ' + names[1] + ' ' + names[2]
            self._init_spherical(symbols)
            if self._cartesian_chart:
                self._transition_spherical_cartesian()
            if self._cylindrical_chart:
                self._transition_spherical_cylindrical()
        return self._spherical_chart

    def spherical_frame(self):
        r"""
        Return the orthonormal vector frame associated with spherical
        coordinates.

        OUTPUT: :class:`~sage.manifolds.differentiable.vectorframe.VectorFrame`

        EXAMPLES::

            sage: E = EuclideanSpace(3)
            sage: E.spherical_frame()
            Vector frame (E^3, (e_r,e_th,e_ph))
            sage: E.spherical_frame()[1]
            Vector field e_r on the Euclidean space E^3
            sage: E.spherical_frame()[:]
            (Vector field e_r on the Euclidean space E^3,
             Vector field e_th on the Euclidean space E^3,
             Vector field e_ph on the Euclidean space E^3)

        The spherical frame expressed in terms of the Cartesian one::

            sage: for e in E.spherical_frame():
            ....:     e.display(E.cartesian_frame(), E.spherical_coordinates())
            e_r = cos(ph)*sin(th) e_x + sin(ph)*sin(th) e_y + cos(th) e_z
            e_th = cos(ph)*cos(th) e_x + cos(th)*sin(ph) e_y - sin(th) e_z
            e_ph = -sin(ph) e_x + cos(ph) e_y

        The orthonormal frame `(e_r, e_\theta, e_\phi)` expressed in terms of
        the coordinate frame
        `\left(\frac{\partial}{\partial r}, \frac{\partial}{\partial\theta},
        \frac{\partial}{\partial\phi}\right)`::

            sage: for e in E.spherical_frame():
            ....:     e.display(E.spherical_coordinates())
            e_r = ∂/∂r
            e_th = 1/r ∂/∂th
            e_ph = 1/(r*sin(th)) ∂/∂ph
        """
        if self._spherical_frame is None:
            # create the spherical chart and the associated orthonormal frame
            self.spherical_coordinates()
        return self._spherical_frame

    def cylindrical_coordinates(self, symbols=None, names=None):
        r"""
        Return the chart of cylindrical coordinates, possibly creating it if it
        does not already exist.

        INPUT:

        - ``symbols`` -- (default: ``None``) string defining the coordinate
          text symbols and LaTeX symbols, with the same conventions as the
          argument ``coordinates`` in
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`; this is
          used only if the cylindrical chart has not been already defined; if
          ``None`` the symbols are generated as `(\rho,\phi,z)`.
        - ``names`` -- (default: ``None``) unused argument, except if
          ``symbols`` is not provided; it must be a tuple containing
          the coordinate symbols (this is guaranteed if the shortcut operator
          ``<,>`` is used)

        OUTPUT:

        - the chart of cylindrical coordinates, as an instance of
          :class:`~sage.manifolds.differentiable.chart.RealDiffChart`

        EXAMPLES::

            sage: E = EuclideanSpace(3)
            sage: E.cylindrical_coordinates()
            Chart (E^3, (rh, ph, z))
            sage: latex(_)
            \left(\mathbb{E}^{3},({\rho}, {\phi}, z)\right)
            sage: E.cylindrical_coordinates().coord_range()
            rh: (0, +oo); ph: [0, 2*pi] (periodic); z: (-oo, +oo)

        The relation to Cartesian coordinates is::

            sage: E.coord_change(E.cylindrical_coordinates(),
            ....:                E.cartesian_coordinates()).display()
            x = rh*cos(ph)
            y = rh*sin(ph)
            z = z
            sage: E.coord_change(E.cartesian_coordinates(),
            ....:                E.cylindrical_coordinates()).display()
            rh = sqrt(x^2 + y^2)
            ph = arctan2(y, x)
            z = z

        The coordinate variables are returned by the square bracket operator::

            sage: E.cylindrical_coordinates()[1]
            rh
            sage: E.cylindrical_coordinates()[3]
            z
            sage: E.cylindrical_coordinates()[:]
            (rh, ph, z)

        They can also be obtained via the operator ``<,>``::

            sage: cylindrical.<rh,ph,z> = E.cylindrical_coordinates()
            sage: cylindrical
            Chart (E^3, (rh, ph, z))
            sage: rh, ph, z
            (rh, ph, z)

        Actually, ``cylindrical.<rh,ph,z> = E.cylindrical_coordinates()`` is a
        shortcut for::

            sage: cylindrical = E.cylindrical_coordinates()
            sage: rh, ph, z = cylindrical[:]

        The coordinate symbols can be customized::

            sage: E = EuclideanSpace(3)
            sage: E.cylindrical_coordinates(symbols=r"R Phi:\Phi Z")
            Chart (E^3, (R, Phi, Z))
            sage: latex(E.cylindrical_coordinates())
            \left(\mathbb{E}^{3},(R, {\Phi}, Z)\right)

        Note that if the cylindrical coordinates have been already initialized,
        the argument ``symbols`` has no effect::

            sage: E.cylindrical_coordinates(symbols=r"rh:\rho ph:\phi z")
            Chart (E^3, (R, Phi, Z))
        """
        if self._cylindrical_chart is None:
            if symbols is None:
                if names is None:
                    symbols = 'rh:\\rho ph:\\phi z'
                else:
                    names = list(names)
                    if names[0] in ['rh', 'rho']:
                        names[0] = names[0] + ':\\rho'
                    if names[1] in ['p', 'ph', 'phi']:
                        names[1] = names[1] + ':\\phi'
                    symbols = names[0] + ' ' + names[1] + ' ' + names[2]
            self._init_cylindrical(symbols)
            if self._cartesian_chart:
                self._transition_cylindrical_cartesian()
            if self._spherical_chart:
                self._transition_spherical_cylindrical()
        return self._cylindrical_chart

    def cylindrical_frame(self):
        r"""
        Return the orthonormal vector frame associated with cylindrical
        coordinates.

        OUTPUT: :class:`~sage.manifolds.differentiable.vectorframe.VectorFrame`

        EXAMPLES::

            sage: E = EuclideanSpace(3)
            sage: E.cylindrical_frame()
            Vector frame (E^3, (e_rh,e_ph,e_z))
            sage: E.cylindrical_frame()[1]
            Vector field e_rh on the Euclidean space E^3
            sage: E.cylindrical_frame()[:]
            (Vector field e_rh on the Euclidean space E^3,
             Vector field e_ph on the Euclidean space E^3,
             Vector field e_z on the Euclidean space E^3)

        The cylindrical frame expressed in terms of the Cartesian one::

            sage: for e in E.cylindrical_frame():
            ....:     e.display(E.cartesian_frame(), E.cylindrical_coordinates())
            e_rh = cos(ph) e_x + sin(ph) e_y
            e_ph = -sin(ph) e_x + cos(ph) e_y
            e_z = e_z

        The orthonormal frame `(e_r, e_\phi, e_z)` expressed in terms of
        the coordinate frame
        `\left(\frac{\partial}{\partial r}, \frac{\partial}{\partial\phi},
        \frac{\partial}{\partial z}\right)`::

            sage: for e in E.cylindrical_frame():
            ....:     e.display(E.cylindrical_coordinates())
            e_rh = ∂/∂rh
            e_ph = 1/rh ∂/∂ph
            e_z = ∂/∂z
        """
        if self._cylindrical_frame is None:
            # create the cylindrical chart and the associated orthonormal frame
            self.cylindrical_coordinates()
        return self._cylindrical_frame

    def scalar_triple_product(self, name=None, latex_name=None):
        r"""
        Return the scalar triple product operator, as a 3-form.

        The *scalar triple product* (also called *mixed product*) of three
        vector fields `u`, `v` and `w` defined on a Euclidean space `E`
        is the scalar field

        .. MATH::

            \epsilon(u,v,w) = u \cdot (v \times w).

        The scalar triple product operator `\epsilon` is a *3-form*, i.e. a
        field of fully antisymmetric trilinear forms; it is also called the
        *volume form* of `E` or the *Levi-Civita tensor* of `E`.

        INPUT:

        - ``name`` -- (default: ``None``) string; name given to the scalar
          triple product operator; if ``None``, ``'epsilon'`` is used
        - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote
          the scalar triple product; if ``None``, it is set to ``r'\epsilon'``
          if ``name`` is ``None`` and to ``name`` otherwise.

        OUTPUT:

        - the scalar triple product operator `\epsilon`, as an instance of
          :class:`~sage.manifolds.differentiable.diff_form.DiffFormParal`

        EXAMPLES::

            sage: E.<x,y,z> = EuclideanSpace()
            sage: triple_product = E.scalar_triple_product()
            sage: triple_product
            3-form epsilon on the Euclidean space E^3
            sage: latex(triple_product)
            \epsilon
            sage: u = E.vector_field(x, y, z, name='u')
            sage: v = E.vector_field(-y, x, 0, name='v')
            sage: w = E.vector_field(y*z, x*z, x*y, name='w')
            sage: s = triple_product(u, v, w); s
            Scalar field epsilon(u,v,w) on the Euclidean space E^3
            sage: s.display()
            epsilon(u,v,w): E^3 → ℝ
               (x, y, z) ↦ x^3*y + x*y^3 - 2*x*y*z^2
            sage: s.expr()
            x^3*y + x*y^3 - 2*x*y*z^2
            sage: latex(s)
            \epsilon\left(u,v,w\right)
            sage: s == - triple_product(w, v, u)
            True

        Check of the identity `\epsilon(u,v,w) = u\cdot(v\times w)`::

            sage: s == u.dot(v.cross(w))
            True

        Customizing the name::

            sage: E.scalar_triple_product(name='S')
            3-form S on the Euclidean space E^3
            sage: latex(_)
            S
            sage: E.scalar_triple_product(name='Omega', latex_name=r'\Omega')
            3-form Omega on the Euclidean space E^3
            sage: latex(_)
            \Omega
        """
        eps = self.volume_form()
        if latex_name is None:
            if name is None:
                latex_name = r'\epsilon'
            else:
                latex_name = name
        if name is None:
            name = 'epsilon'
        eps.set_name(name=name, latex_name=latex_name)
        return eps
