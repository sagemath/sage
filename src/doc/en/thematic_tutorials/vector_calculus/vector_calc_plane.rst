.. sage-doctest: needs sage.plot sage.symbolic

.. linkall

Vector calculus in the Euclidean plane
======================================

This tutorial introduces some vector calculus capabilities of SageMath in the
framework of the 2-dimensional Euclidean space. The corresponding tools have
been developed via the `SageManifolds <https://sagemanifolds.obspm.fr>`__
project.

The tutorial is also available as a Jupyter notebook, either
`passive <https://nbviewer.jupyter.org/github/sagemanifolds/SageManifolds/blob/master/Notebooks/VectorCalculus/vector_calc_plane.ipynb>`__ (``nbviewer``)
or `interactive <https://mybinder.org/v2/gh/sagemanifolds/SageManifolds/master?filepath=Notebooks/VectorCalculus/vector_calc_plane.ipynb>`__ (``binder``).


1. Defining the Euclidean plane
-------------------------------

We define the Euclidean plane `\mathbb{E}^2` as a 2-dimensional Euclidean
space, with Cartesian coordinates `(x,y)`, via the function
:class:`EuclideanSpace`::

    sage: E.<x,y> = EuclideanSpace()
    sage: E
    Euclidean plane E^2

Thanks to the use of ``<x,y>`` in the above command, the Python variables
``x`` and ``y`` are assigned to the symbolic variables `x` and `y` describing
the Cartesian coordinates (there is no need to declare them via :func:`var`,
i.e. to type ``x, y = var('x y')``)::

    sage: type(y)
    <class 'sage.symbolic.expression.Expression'>

Instead of using the variables ``x`` and ``y``, one may also access to the
coordinates by their indices in the chart of Cartesian coordinates::

    sage: cartesian = E.cartesian_coordinates()
    sage: cartesian
    Chart (E^2, (x, y))

::

    sage: cartesian[1]
    x
    sage: cartesian[2]
    y
    sage: y is cartesian[2]
    True

Each of the Cartesian coordinates spans the entire real line::

    sage: cartesian.coord_range()
    x: (-oo, +oo); y: (-oo, +oo)


2. Vector fields
----------------

The Euclidean plane `\mathbb{E}^2` is canonically endowed with the vector
frame associated with Cartesian coordinates::

    sage: E.default_frame()
    Coordinate frame (E^2, (e_x,e_y))

Vector fields on `\mathbb{E}^2` are then defined from their components in that
frame::

    sage: v = E.vector_field(-y, x, name='v')
    sage: v.display()
    v = -y e_x + x e_y

The access to individual components is performed by the square bracket
operator::

    sage: v[1]
    -y
    sage: v[:]
    [-y, x]

A plot of the vector field `v` (this is with the default parameters, see the
documentation of
:meth:`~sage.manifolds.differentiable.vectorfield.VectorField.plot` for the
various options)::

    sage: v.plot()
    Graphics object consisting of 80 graphics primitives

.. PLOT::

    E = EuclideanSpace(2)
    X = E.default_chart(); x, y = X[:]
    v = E.vector_field(-y, x, name='v')
    g = v.plot()
    sphinx_plot(g)

One may also define a vector field by setting the components in a second
stage::

    sage: w = E.vector_field(name='w')
    sage: w[1] = function('w_x')(x,y)
    sage: w[2] = function('w_y')(x,y)
    sage: w.display()
    w = w_x(x, y) e_x + w_y(x, y) e_y

Note that in the above example the components of `w` are unspecified functions
of `(x,y)`, contrary to the components of `v`.

Standard linear algebra operations can be performed on vector fields::

    sage: s = 2*v + x*w
    sage: s.display()
    (x*w_x(x, y) - 2*y) e_x + (x*w_y(x, y) + 2*x) e_y


Scalar product and norm
~~~~~~~~~~~~~~~~~~~~~~~

The dot (or scalar) product `u\cdot v` of the vector fields `u` and `v` is
obtained by the operator
:meth:`~sage.manifolds.differentiable.vectorfield.VectorField.dot_product`; it
gives rise to a scalar field on `\mathbb{E}^2`::

    sage: s = v.dot_product(w)
    sage: s
    Scalar field v.w on the Euclidean plane E^2

A shortcut alias of
:meth:`~sage.manifolds.differentiable.vectorfield.VectorField.dot_product` is
``dot``::

    sage: s == v.dot(w)
    True

::

    sage: s.display()
    v.w: E^2 → ℝ
       (x, y) ↦ -y*w_x(x, y) + x*w_y(x, y)

The symbolic expression representing the scalar field `v\cdot w` is obtained
by means of the method :meth:`~sage.manifolds.scalarfield.ScalarField.expr`::

    sage: s.expr()
    -y*w_x(x, y) + x*w_y(x, y)

The Euclidean norm of the vector field `v` is a scalar field on
`\mathbb{E}^2`::

    sage: s = norm(v)
    sage: s.display()
    |v|: E^2 → ℝ
       (x, y) ↦ sqrt(x^2 + y^2)

Again, the corresponding symbolic expression is obtained via
:meth:`~sage.manifolds.scalarfield.ScalarField.expr`::

    sage: s.expr()
    sqrt(x^2 + y^2)

::

    sage: norm(w).expr()
    sqrt(w_x(x, y)^2 + w_y(x, y)^2)

We have of course `\|v\|^2 = v\cdot v` ::

    sage: norm(v)^2 == v.dot(v)
    True


Values at a given point
~~~~~~~~~~~~~~~~~~~~~~~

We introduce a point `p\in \mathbb{E}^2` via the generic SageMath syntax for
creating an element from its parent (here `\mathbb{E}^2`), i.e. the call
operator ``()``, with the Cartesian coordinates of the point as the first
argument::

    sage: p = E((-2,3), name='p')
    sage: p
    Point p on the Euclidean plane E^2

The coordinates of `p` are returned by the method
:meth:`~sage.manifolds.point.ManifoldPoint.coord`::

    sage: p.coord()
    (-2, 3)

or by letting the chart ``cartesian`` act on the point::

    sage: cartesian(p)
    (-2, 3)

The value of the scalar field ``s = norm(v)`` at `p` is::

    sage: s(p)
    sqrt(13)

The value of a vector field at `p` is obtained by the method
:meth:`~sage.manifolds.differentiable.tensorfield_paral.TensorFieldParal.at`
(since the call operator ``()`` is reserved for the action on scalar fields,
see :ref:`vector_fields_as_derivations_plane` below)::

    sage: vp = v.at(p)
    sage: vp
    Vector v at Point p on the Euclidean plane E^2
    sage: vp.display()
    v = -3 e_x - 2 e_y
    sage: wp = w.at(p)
    sage: wp.display()
    w = w_x(-2, 3) e_x + w_y(-2, 3) e_y
    sage: s = v.at(p) + pi*w.at(p)
    sage: s.display()
    (pi*w_x(-2, 3) - 3) e_x + (pi*w_y(-2, 3) - 2) e_y


3. Differential operators
-------------------------

The standard operators `\mathrm{grad}`, `\mathrm{div}`, etc. involved in
vector calculus are accessible as methods on scalar fields and vector fields
(e.g. ``v.div()``). However, to use standard mathematical notations (e.g.
``div(v)``), let us import the functions
:func:`~sage.manifolds.operators.grad`, :func:`~sage.manifolds.operators.div`,
and :func:`~sage.manifolds.operators.laplacian` in the global namespace::

    sage: from sage.manifolds.operators import *


Divergence
~~~~~~~~~~

The divergence of a vector field is returned by the function
:func:`~sage.manifolds.operators.div`; the output is a scalar field on
`\mathbb{E}^2`::

    sage: div(v)
    Scalar field div(v) on the Euclidean plane E^2
    sage: div(v).display()
    div(v): E^2 → ℝ
       (x, y) ↦ 0

In the present case, `\mathrm{div}\, v` vanishes identically::

    sage: div(v) == 0
    True

On the contrary, the divergence of `w` is::

    sage: div(w).display()
    div(w): E^2 → ℝ
       (x, y) ↦ d(w_x)/dx + d(w_y)/dy
    sage: div(w).expr()
    diff(w_x(x, y), x) + diff(w_y(x, y), y)


Gradient
~~~~~~~~

The gradient of a scalar field, e.g. ``s = norm(v)``, is returned by the
function :func:`~sage.manifolds.operators.grad`; the output is a vector field::

    sage: s = norm(v)
    sage: grad(s)
    Vector field grad(|v|) on the Euclidean plane E^2
    sage: grad(s).display()
    grad(|v|) = x/sqrt(x^2 + y^2) e_x + y/sqrt(x^2 + y^2) e_y
    sage: grad(s)[2]
    y/sqrt(x^2 + y^2)

For a generic scalar field, like::

    sage: F = E.scalar_field(function('f')(x,y), name='F')

we have::

    sage: grad(F).display()
    grad(F) = d(f)/dx e_x + d(f)/dy e_y
    sage: grad(F)[:]
    [d(f)/dx, d(f)/dy]

Of course, we may combine :func:`~sage.manifolds.operators.grad` and
:func:`~sage.manifolds.operators.div`::

    sage: grad(div(w)).display()
    grad(div(w)) = (d^2(w_x)/dx^2 + d^2(w_y)/dxdy) e_x + (d^2(w_x)/dxdy + d^2(w_y)/dy^2) e_y


Laplace operator
~~~~~~~~~~~~~~~~

The Laplace operator `\Delta` is obtained by the function
:func:`~sage.manifolds.operators.laplacian`; it acts on scalar fields::

    sage: laplacian(F).display()
    Delta(F): E^2 → ℝ
       (x, y) ↦ d^2(f)/dx^2 + d^2(f)/dy^2

as well as on vector fields::

    sage: laplacian(w).display()
    Delta(w) = (d^2(w_x)/dx^2 + d^2(w_x)/dy^2) e_x + (d^2(w_y)/dx^2 + d^2(w_y)/dy^2) e_y

For a scalar field, we have the identity

.. MATH::

    \Delta F = \mathrm{div}\left(\mathrm{grad}\, F\right),

as we can check::

    sage: laplacian(F) == div(grad(F))
    True


4. Polar coordinates
--------------------

Polar coordinates `(r,\phi)` are introduced on `\mathbb{E}^2` by::

    sage: polar.<r,ph> = E.polar_coordinates()
    sage: polar
    Chart (E^2, (r, ph))
    sage: polar.coord_range()
    r: (0, +oo); ph: [0, 2*pi] (periodic)

They are related to Cartesian coordinates by the following transformations::

    sage: E.coord_change(polar, cartesian).display()
    x = r*cos(ph)
    y = r*sin(ph)
    sage: E.coord_change(cartesian, polar).display()
    r = sqrt(x^2 + y^2)
    ph = arctan2(y, x)

The orthonormal vector frame `(e_r, e_\phi)` associated with polar coordinates
is returned by the method
:meth:`~sage.manifolds.differentiable.euclidean.EuclideanPlane.polar_frame`::

    sage: polar_frame = E.polar_frame()
    sage: polar_frame
    Vector frame (E^2, (e_r,e_ph))

::

    sage: er = polar_frame[1]
    sage: er.display()
    e_r = x/sqrt(x^2 + y^2) e_x + y/sqrt(x^2 + y^2) e_y

The above display is in the default frame (Cartesian frame) with the default
coordinates (Cartesian). Let us ask for the display in the same frame, but
with the components expressed in polar coordinates::

    sage: er.display(cartesian.frame(), polar)
    e_r = cos(ph) e_x + sin(ph) e_y

Similarly::

    sage: eph = polar_frame[2]
    sage: eph.display()
    e_ph = -y/sqrt(x^2 + y^2) e_x + x/sqrt(x^2 + y^2) e_y
    sage: eph.display(cartesian.frame(), polar)
    e_ph = -sin(ph) e_x + cos(ph) e_y

We may check that `(e_r, e_\phi)` is an orthonormal frame::

    sage: all([er.dot(er) == 1, er.dot(eph) == 0, eph.dot(eph) == 1])
    True

Scalar fields can be expressed in terms of polar coordinates::

    sage: F.display()
    F: E^2 → ℝ
       (x, y) ↦ f(x, y)
       (r, ph) ↦ f(r*cos(ph), r*sin(ph))
    sage: F.display(polar)
    F: E^2 → ℝ
       (r, ph) ↦ f(r*cos(ph), r*sin(ph))

and we may ask for the components of vector fields in terms of the polar
frame::

    sage: v.display()  # default frame and default coordinates (both Cartesian ones)
    v = -y e_x + x e_y
    sage: v.display(polar_frame)  # polar frame and default coordinates
    v = sqrt(x^2 + y^2) e_ph
    sage: v.display(polar_frame, polar)  # polar frame and polar coordinates
    v = r e_ph

::

    sage: w.display()
    w = w_x(x, y) e_x + w_y(x, y) e_y
    sage: w.display(polar_frame, polar)
    w = (cos(ph)*w_x(r*cos(ph), r*sin(ph)) + sin(ph)*w_y(r*cos(ph), r*sin(ph))) e_r
    + (-sin(ph)*w_x(r*cos(ph), r*sin(ph)) + cos(ph)*w_y(r*cos(ph), r*sin(ph))) e_ph


Gradient in polar coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let us define a generic scalar field in terms of polar coordinates::

    sage: H = E.scalar_field({polar: function('h')(r,ph)}, name='H')
    sage: H.display(polar)
    H: E^2 → ℝ
       (r, ph) ↦ h(r, ph)

The gradient of `H` is then::

    sage: grad(H).display(polar_frame, polar)
    grad(H) = d(h)/dr e_r + d(h)/dph/r e_ph

The access to individual components is achieved via the square bracket
operator, where, in addition to the index, one has to specify the vector frame
and the coordinates if they are not the default ones::

    sage: grad(H).display(cartesian.frame(), polar)
    grad(H) = (r*cos(ph)*d(h)/dr - sin(ph)*d(h)/dph)/r e_x + (r*sin(ph)*d(h)/dr
     + cos(ph)*d(h)/dph)/r e_y
    sage: grad(H)[polar_frame, 2, polar]
    d(h)/dph/r


Divergence in polar coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let us define a generic vector field in terms of polar coordinates::

    sage: u = E.vector_field(function('u_r')(r,ph),
    ....:                    function('u_ph', latex_name=r'u_\phi')(r,ph),
    ....:                    frame=polar_frame, chart=polar, name='u')
    sage: u.display(polar_frame, polar)
    u = u_r(r, ph) e_r + u_ph(r, ph) e_ph

Its divergence is::

    sage: div(u).display(polar)
    div(u): E^2 → ℝ
       (r, ph) ↦ (r*d(u_r)/dr + u_r(r, ph) + d(u_ph)/dph)/r
    sage: div(u).expr(polar)
    (r*diff(u_r(r, ph), r) + u_r(r, ph) + diff(u_ph(r, ph), ph))/r
    sage: div(u).expr(polar).expand()
    u_r(r, ph)/r + diff(u_ph(r, ph), ph)/r + diff(u_r(r, ph), r)


Using polar coordinates by default:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to avoid specifying the arguments ``polar_frame`` and ``polar`` in
``display()``, ``expr()`` and ``[]``, we may change the default values by
means of
:meth:`~sage.manifolds.manifold.TopologicalManifold.set_default_chart` and
:meth:`~sage.manifolds.differentiable.manifold.DifferentiableManifold.set_default_frame`::

    sage: E.set_default_chart(polar)
    sage: E.set_default_frame(polar_frame)

Then we have::

    sage: u.display()
    u = u_r(r, ph) e_r + u_ph(r, ph) e_ph
    sage: u[1]
    u_r(r, ph)

::

    sage: v.display()
    v = r e_ph
    sage: v[2]
    r

::

    sage: w.display()
    w = (cos(ph)*w_x(r*cos(ph), r*sin(ph)) + sin(ph)*w_y(r*cos(ph), r*sin(ph))) e_r + (-sin(ph)*w_x(r*cos(ph), r*sin(ph)) + cos(ph)*w_y(r*cos(ph), r*sin(ph))) e_ph
    sage: div(u).expr()
    (r*diff(u_r(r, ph), r) + u_r(r, ph) + diff(u_ph(r, ph), ph))/r


5. Advanced topics: the Euclidean plane as a Riemannian manifold
----------------------------------------------------------------

`\mathbb{E}^2` is actually a *Riemannian manifold* (see
:mod:`~sage.manifolds.differentiable.pseudo_riemannian`), i.e. a smooth real
manifold endowed with a positive definite metric tensor::

    sage: E.category()
    Join of
     Category of smooth manifolds over Real Field with 53 bits of precision and
     Category of connected manifolds over Real Field with 53 bits of precision and
     Category of complete metric spaces
    sage: E.base_field() is RR
    True

Actually ``RR`` is used here as a proxy for the real field (this should
be replaced in the future, see the discussion at
:issue:`24456`) and the 53 bits of
precision play of course no role for the symbolic computations.

The user atlas of `\mathbb{E}^2` has two charts::

    sage: E.atlas()
    [Chart (E^2, (x, y)), Chart (E^2, (r, ph))]

while there are three vector frames defined on `\mathbb{E}^2`::

    sage: E.frames()
    [Coordinate frame (E^2, (e_x,e_y)),
     Coordinate frame (E^2, (∂/∂r,∂/∂ph)),
     Vector frame (E^2, (e_r,e_ph))]

Indeed, there are two frames associated with polar coordinates: the coordinate
frame `(\frac{\partial}{\partial r}, \frac{\partial}{\partial \phi})` and the
orthonormal frame `(e_r, e_\phi)`.

Riemannian metric
~~~~~~~~~~~~~~~~~

The default metric tensor of `\mathbb{E}^2` is::

    sage: g = E.metric()
    sage: g
    Riemannian metric g on the Euclidean plane E^2
    sage: g.display()
    g = e^r⊗e^r + e^ph⊗e^ph

In the above display, ``e^r`` = `e^r` and ``e^ph`` = `e^\phi` are the 1-forms
defining the coframe dual to the orthonormal polar frame `(e_r, e_\phi)`,
which is the default vector frame on `\mathbb{E}^2`::

    sage: polar_frame.coframe()
    Coframe (E^2, (e^r,e^ph))

Of course, we may ask for some display with respect to frames different from
the default one::

    sage: g.display(cartesian.frame())
    g = dx⊗dx + dy⊗dy
    sage: g.display(polar.frame())
    g = dr⊗dr + r^2 dph⊗dph
    sage: g[:]
    [1 0]
    [0 1]
    sage: g[polar.frame(),:]
    [  1   0]
    [  0 r^2]

`g` is a *flat* metric: its (Riemann) curvature tensor (see
:meth:`~sage.manifolds.differentiable.metric.PseudoRiemannianMetric.riemann`)
is zero::

    sage: g.riemann()
    Tensor field Riem(g) of type (1,3) on the Euclidean plane E^2
    sage: g.riemann().display()
    Riem(g) = 0

The metric `g` is defining the dot product on `\mathbb{E}^2`::

    sage: v.dot(w) == g(v,w)
    True
    sage: norm(v) == sqrt(g(v,v))
    True

.. _vector_fields_as_derivations_plane:

Vector fields as derivations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Vector fields act as derivations on scalar fields::

    sage: v(F)
    Scalar field v(F) on the Euclidean plane E^2
    sage: v(F).display()
    v(F): E^2 → ℝ
       (x, y) ↦ -y*d(f)/dx + x*d(f)/dy
       (r, ph) ↦ -r*sin(ph)*d(f)/d(r*cos(ph)) + r*cos(ph)*d(f)/d(r*sin(ph))
    sage: v(F) == v.dot(grad(F))
    True

::

    sage: dF = F.differential()
    sage: dF
    1-form dF on the Euclidean plane E^2
    sage: v(F) == dF(v)
    True

The set `\mathfrak{X}(\mathbb{E}^2)` of all vector fields on `\mathbb{E}^2` is
a free module of rank 2 over the commutative algebra of smooth scalar fields
on `\mathbb{E}^2`, `C^\infty(\mathbb{E}^2)`::

    sage: XE = v.parent()
    sage: XE
    Free module X(E^2) of vector fields on the Euclidean plane E^2
    sage: XE.category()
    Category of finite dimensional modules over Algebra of differentiable
     scalar fields on the Euclidean plane E^2
    sage: XE.base_ring()
    Algebra of differentiable scalar fields on the Euclidean plane E^2

::

    sage: CE = F.parent()
    sage: CE
    Algebra of differentiable scalar fields on the Euclidean plane E^2
    sage: CE is XE.base_ring()
    True
    sage: CE.category()
    Join of Category of commutative algebras over Symbolic Ring and Category of homsets of topological spaces
    sage: rank(XE)
    2

The bases of the free module `\mathfrak{X}(\mathbb{E}^2)` are nothing but the
vector frames defined on `\mathbb{E}^2`::

    sage: XE.bases()
    [Coordinate frame (E^2, (e_x,e_y)),
     Coordinate frame (E^2, (∂/∂r,∂/∂ph)),
     Vector frame (E^2, (e_r,e_ph))]


Tangent spaces
~~~~~~~~~~~~~~

A vector field evaluated at a point $p$ is a vector in the tangent space
`T_p\mathbb{E}^2`::

    sage: vp = v.at(p)
    sage: vp.display()
    v = -3 e_x - 2 e_y

::

    sage: Tp = vp.parent()
    sage: Tp
    Tangent space at Point p on the Euclidean plane E^2
    sage: Tp.category()
    Category of finite dimensional vector spaces over Symbolic Ring
    sage: dim(Tp)
    2
    sage: isinstance(Tp, FiniteRankFreeModule)
    True
    sage: sorted(Tp.bases(), key=str)
    [Basis (e_r,e_ph) on the Tangent space at Point p on the Euclidean plane E^2,
     Basis (e_x,e_y) on the Tangent space at Point p on the Euclidean plane E^2]


Levi-Civita connection
~~~~~~~~~~~~~~~~~~~~~~

The Levi-Civita connection associated to the Euclidean metric `g` is::

    sage: nabla = g.connection()
    sage: nabla
    Levi-Civita connection nabla_g associated with the Riemannian metric g on the Euclidean plane E^2

The corresponding Christoffel symbols with respect to the polar coordinates
are::

    sage: g.christoffel_symbols_display()
    Gam^r_ph,ph = -r
    Gam^ph_r,ph = 1/r

By default, only nonzero and nonredundant values are displayed (for instance
`\Gamma^\phi_{\ \, \phi r}` is skipped, since it can be deduced from
`\Gamma^\phi_{\ \, r \phi}` by symmetry on the last two indices).

The Christoffel symbols with respect to the Cartesian coordinates are all
zero::

    sage: g.christoffel_symbols_display(chart=cartesian, only_nonzero=False)
    Gam^x_xx = 0
    Gam^x_xy = 0
    Gam^x_yy = 0
    Gam^y_xx = 0
    Gam^y_xy = 0
    Gam^y_yy = 0

`\nabla_g` is the connection involved in differential operators::

    sage: grad(F) == nabla(F).up(g)
    True
    sage: nabla(F) == grad(F).down(g)
    True
    sage: div(v) == nabla(v).trace()
    True
    sage: div(w) == nabla(w).trace()
    True
    sage: laplacian(F) == nabla(nabla(F).up(g)).trace()
    True
    sage: laplacian(w) == nabla(nabla(w).up(g)).trace(1,2)
    True
