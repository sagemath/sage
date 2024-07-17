.. .. sage-doctest: needs sage.plot sage.symbolic

.. linkall

.. _change_coord_euclidean:

How to change coordinates
=========================

This tutorial introduces some vector calculus capabilities of SageMath within
the 3-dimensional Euclidean space. The corresponding tools have been developed
via the `SageManifolds <https://sagemanifolds.obspm.fr>`__ project.

The tutorial is also available as a Jupyter notebook, either
`passive <https://nbviewer.jupyter.org/github/sagemanifolds/SageManifolds/blob/master/Notebooks/VectorCalculus/vector_calc_change.ipynb>`__ (``nbviewer``)
or `interactive <https://mybinder.org/v2/gh/sagemanifolds/SageManifolds/master?filepath=Notebooks/VectorCalculus/vector_calc_change.ipynb>`__ (``binder``).


Starting from Cartesian coordinates
-----------------------------------

In this tutorial, we choose to start from the Cartesian coordinates `(x,y,z)`.
Hence, we declare the 3-dimensional Euclidean space `\mathbb{E}^3` as::

    sage: E.<x,y,z> = EuclideanSpace()
    sage: E
    Euclidean space E^3

By default, i.e. without the optional argument ``coordinates`` in
:class:`EuclideanSpace`, `\mathbb{E}^3` is initialized with the chart of
Cartesian coordinates::

    sage: E.atlas()
    [Chart (E^3, (x, y, z))]

See the tutorial :ref:`vector_calc_curvilinear` for examples of initialization
of the Euclidean space with spherical coordinates or cylindrical coordinates
instead of the Cartesian ones.

Let us denote by ``cartesian`` the chart of Cartesian coordinates::

    sage: cartesian = E.cartesian_coordinates()
    sage: cartesian
    Chart (E^3, (x, y, z))

The access to the individual coordinates is performed via the square
bracket operator::

    sage: cartesian[1]
    x
    sage: cartesian[:]
    (x, y, z)

Thanks to use of ``<x,y,z>`` when declaring ``E``, the Python variables ``x``,
``y`` and ``z`` have been created to store the coordinates `(x,y,z)` as
symbolic expressions. There is no need to declare them via :func:`var`, i.e. to
type ``x, y, z = var('x y z')``; they are immediately available::

    sage: y is cartesian[2]
    True
    sage: type(y)
    <class 'sage.symbolic.expression.Expression'>

Each of the Cartesian coordinates spans the entire real line::

    sage: cartesian.coord_range()
    x: (-oo, +oo); y: (-oo, +oo); z: (-oo, +oo)

Being the only coordinate chart created so far, ``cartesian`` is the default
chart on ``E``::

    sage: cartesian is E.default_chart()
    True

`\mathbb{E}^3` is endowed with the orthonormal vector frame `(e_x, e_y, e_z)`
associated with Cartesian coordinates::

    sage: E.frames()
    [Coordinate frame (E^3, (e_x,e_y,e_z))]

Let us denote it by ``cartesian_frame``::

    sage: cartesian_frame = E.cartesian_frame()
    sage: cartesian_frame
    Coordinate frame (E^3, (e_x,e_y,e_z))
    sage: cartesian_frame is E.default_frame()
    True

Each element of this frame is a unit vector field; for instance, we have
`e_x\cdot e_x = 1`::

    sage: e_x = cartesian_frame[1]
    sage: e_x
    Vector field e_x on the Euclidean space E^3
    sage: e_x.dot(e_x).expr()
    1

as well as `e_x\cdot e_y = 0`::

    sage: e_y = cartesian_frame[2]
    sage: e_x.dot(e_y).expr()
    0


Introducing spherical coordinates
---------------------------------

Spherical coordinates are introduced by::

    sage: spherical.<r,th,ph> = E.spherical_coordinates()
    sage: spherical
    Chart (E^3, (r, th, ph))

We have::

    sage: spherical[:]
    (r, th, ph)
    sage: spherical.coord_range()
    r: (0, +oo); th: (0, pi); ph: [0, 2*pi] (periodic)

`\mathbb{E}^3` is now endowed with two coordinate charts::

    sage: E.atlas()
    [Chart (E^3, (x, y, z)), Chart (E^3, (r, th, ph))]

The change-of-coordinate formulas have been automatically implemented during
the above call ``E.spherical_coordinates()``::

    sage: E.coord_change(spherical, cartesian).display()
    x = r*cos(ph)*sin(th)
    y = r*sin(ph)*sin(th)
    z = r*cos(th)
    sage: E.coord_change(cartesian, spherical).display()
    r = sqrt(x^2 + y^2 + z^2)
    th = arctan2(sqrt(x^2 + y^2), z)
    ph = arctan2(y, x)

These formulas are automatically used if we ask to plot the grid of spherical
coordinates in terms of Cartesian coordinates::

    sage: spherical.plot(cartesian, color={r:'red', th:'green', ph:'orange'})
    Graphics3d Object

.. PLOT::

    E = EuclideanSpace(3)
    cartesian = E.cartesian_coordinates()
    spherical = E.spherical_coordinates()
    r, th, ph = spherical[:]
    g = spherical.plot(cartesian, color={r:'red', th:'green', ph:'orange'})
    sphinx_plot(g)

Note that

- the red lines are those along which `r` varies, while
  `(\theta,\phi)` are kept fixed;
- the grid lines are those along which `\theta` varies, while
  `(r,\phi)` are kept fixed;
- the orange lines are those along which `\phi` varies, while
  `(r,\theta)` are kept fixed.

For customizing the plot, see the list of options in the documentation of
:meth:`~sage.manifolds.chart.RealChart.plot`. For instance, we may draw the
spherical coordinates in the plane `\theta=\pi/2` in terms of the coordinates
`(x, y)`::

    sage: spherical.plot(cartesian, fixed_coords={th: pi/2}, ambient_coords=(x,y),
    ....:                color={r:'red', th:'green', ph:'orange'})
    Graphics object consisting of 18 graphics primitives

.. PLOT::

    E = EuclideanSpace(3)
    cartesian = E.cartesian_coordinates()
    spherical = E.spherical_coordinates()
    x, y, z = cartesian[:]
    r, th, ph = spherical[:]
    g = spherical.plot(cartesian, fixed_coords={th: pi/2}, ambient_coords=(x,y),
                       color={r:'red', th:'green', ph:'orange'})
    sphinx_plot(g)

Similarly the grid of spherical coordinates in the half-plane `\phi=0`
drawn in terms of the coordinates `(x, z)` is obtained via::

    sage: spherical.plot(cartesian, fixed_coords={ph: 0}, ambient_coords=(x,z),
    ....:                color={r:'red', th:'green', ph:'orange'})
    Graphics object consisting of 18 graphics primitives

.. PLOT::

    E = EuclideanSpace(3)
    cartesian = E.cartesian_coordinates()
    spherical = E.spherical_coordinates()
    x, y, z = cartesian[:]
    r, th, ph = spherical[:]
    g = spherical.plot(cartesian, fixed_coords={ph: 0}, ambient_coords=(x,z),
                       color={r:'red', th:'green', ph:'orange'})
    sphinx_plot(g)

Relations between the Cartesian and spherical vector frames
-----------------------------------------------------------

At this stage, `\mathbb{E}^3` is endowed with three vector frames::

    sage: E.frames()
    [Coordinate frame (E^3, (e_x,e_y,e_z)),
     Coordinate frame (E^3, (∂/∂r,∂/∂th,∂/∂ph)),
     Vector frame (E^3, (e_r,e_th,e_ph))]

The second one is the *coordinate* frame `\left(\frac{\partial}{\partial r},
\frac{\partial}{\partial\theta}, \frac{\partial}{\partial \phi}\right)` of
spherical coordinates, while the third one is the standard *orthonormal* frame
`(e_r,e_\theta,e_\phi)` associated with spherical coordinates. For Cartesian
coordinates, the coordinate frame and the orthonormal frame coincide: it is
`(e_x,e_y,e_z)`. For spherical coordinates, the orthonormal frame is returned
by the method
:meth:`~sage.manifolds.differentiable.euclidean.Euclidean3dimSpace.spherical_frame`::

    sage: spherical_frame = E.spherical_frame()
    sage: spherical_frame
    Vector frame (E^3, (e_r,e_th,e_ph))

We may check that it is an orthonormal frame, i.e. that it obeys
`e_i\cdot e_j = \delta_{ij}`::

    sage: es = spherical_frame
    sage: [[es[i].dot(es[j]).expr() for j in E.irange()] for i in E.irange()]
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

Via the method ``display``, we may express the orthonormal spherical frame in
terms of the Cartesian one::

    sage: for vec in spherical_frame:
    ....:     vec.display(cartesian_frame, spherical)
    e_r = cos(ph)*sin(th) e_x + sin(ph)*sin(th) e_y + cos(th) e_z
    e_th = cos(ph)*cos(th) e_x + cos(th)*sin(ph) e_y - sin(th) e_z
    e_ph = -sin(ph) e_x + cos(ph) e_y


The reverse is::

    sage: for vec in cartesian_frame:
    ....:     vec.display(spherical_frame, spherical)
    e_x = cos(ph)*sin(th) e_r + cos(ph)*cos(th) e_th - sin(ph) e_ph
    e_y = sin(ph)*sin(th) e_r + cos(th)*sin(ph) e_th + cos(ph) e_ph
    e_z = cos(th) e_r - sin(th) e_th

We may also express the orthonormal frame `(e_r,e_\theta,e_\phi)` in terms on
the coordinate frame `\left(\frac{\partial}{\partial r},
\frac{\partial}{\partial\theta}, \frac{\partial}{\partial \phi}\right)` (the
latter being returned by the method
:meth:`~sage.manifolds.differentiable.chart.DiffChart.frame` acting on the
chart ``spherical``)::

    sage: for vec in spherical_frame:
    ....:     vec.display(spherical.frame(), spherical)
    e_r = ∂/∂r
    e_th = 1/r ∂/∂th
    e_ph = 1/(r*sin(th)) ∂/∂ph


Introducing cylindrical coordinates
-----------------------------------

Cylindrical coordinates are introduced in a way similar to spherical
coordinates::

    sage: cylindrical.<rh,ph,z> = E.cylindrical_coordinates()
    sage: cylindrical
    Chart (E^3, (rh, ph, z))

We have::

    sage: cylindrical[:]
    (rh, ph, z)
    sage: rh is cylindrical[1]
    True
    sage: cylindrical.coord_range()
    rh: (0, +oo); ph: [0, 2*pi] (periodic); z: (-oo, +oo)

`\mathbb{E}^3` is now endowed with three coordinate charts::

    sage: E.atlas()
    [Chart (E^3, (x, y, z)), Chart (E^3, (r, th, ph)), Chart (E^3, (rh, ph, z))]

The transformations linking the cylindrical coordinates to the Cartesian ones
are::

    sage: E.coord_change(cylindrical, cartesian).display()
    x = rh*cos(ph)
    y = rh*sin(ph)
    z = z
    sage: E.coord_change(cartesian, cylindrical).display()
    rh = sqrt(x^2 + y^2)
    ph = arctan2(y, x)
    z = z

There are now five vector frames defined on `\mathbb{E}^3`::

    sage: E.frames()
    [Coordinate frame (E^3, (e_x,e_y,e_z)),
     Coordinate frame (E^3, (∂/∂r,∂/∂th,∂/∂ph)),
     Vector frame (E^3, (e_r,e_th,e_ph)),
     Coordinate frame (E^3, (∂/∂rh,∂/∂ph,∂/∂z)),
     Vector frame (E^3, (e_rh,e_ph,e_z))]

The orthonormal frame associated with cylindrical coordinates is
`(e_\rho, e_\phi, e_z)`::

    sage: cylindrical_frame = E.cylindrical_frame()
    sage: cylindrical_frame
    Vector frame (E^3, (e_rh,e_ph,e_z))

We may check that it is an orthonormal frame::

    sage: ec = cylindrical_frame
    sage: [[ec[i].dot(ec[j]).expr() for j in E.irange()] for i in E.irange()]
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

and express it in terms of the Cartesian frame::

    sage: for vec in cylindrical_frame:
    ....:     vec.display(cartesian_frame, cylindrical)
    e_rh = cos(ph) e_x + sin(ph) e_y
    e_ph = -sin(ph) e_x + cos(ph) e_y
    e_z = e_z

The reverse is::

    sage: for vec in cartesian_frame:
    ....:     vec.display(cylindrical_frame, cylindrical)
    e_x = cos(ph) e_rh - sin(ph) e_ph
    e_y = sin(ph) e_rh + cos(ph) e_ph
    e_z = e_z

Of course, we may express the orthonormal cylindrical frame in terms of the
spherical one::

    sage: for vec in cylindrical_frame:
    ....:     vec.display(spherical_frame, spherical)
    e_rh = sin(th) e_r + cos(th) e_th
    e_ph = e_ph
    e_z = cos(th) e_r - sin(th) e_th

along with the reverse transformation::

    sage: for vec in spherical_frame:
    ....:     vec.display(cylindrical_frame, spherical)
    e_r = sin(th) e_rh + cos(th) e_z
    e_th = cos(th) e_rh - sin(th) e_z
    e_ph = e_ph

The orthonormal frame `(e_\rho,e_\phi,e_z)` can be expressed in terms on the
coordinate frame `\left(\frac{\partial}{\partial\rho},
\frac{\partial}{\partial\phi}, \frac{\partial}{\partial z}\right)` (the latter
being returned by the method
:meth:`~sage.manifolds.differentiable.chart.DiffChart.frame` acting on the
chart ``cylindrical``)::

    sage: for vec in cylindrical_frame:
    ....:     vec.display(cylindrical.frame(), cylindrical)
    e_rh = ∂/∂rh
    e_ph = 1/rh ∂/∂ph
    e_z = ∂/∂z


How to evaluate the coordinates of a point in various systems
-------------------------------------------------------------

Let us introduce a point `p\in \mathbb{E}^3` via the generic SageMath
syntax for creating an element from its parent (here
`\mathbb{E}^3`), i.e. the call operator ``()``, with the
coordinates of the point as the first argument::

    sage: p = E((-1, 1,0), chart=cartesian, name='p')
    sage: p
    Point p on the Euclidean space E^3

Actually, since the Cartesian coordinates are the default ones, the argument
``chart=cartesian`` can be omitted::

    sage: p = E((-1, 1,0), name='p')
    sage: p
    Point p on the Euclidean space E^3

The coordinates of `p` in a given coordinate chart are obtained by letting the
corresponding chart act on `p`::

    sage: cartesian(p)
    (-1, 1, 0)
    sage: spherical(p)
    (sqrt(2), 1/2*pi, 3/4*pi)
    sage: cylindrical(p)
    (sqrt(2), 3/4*pi, 0)

Here some example of a point defined from its spherical coordinates::

    sage: q = E((4,pi/3,pi), chart=spherical, name='q')
    sage: q
    Point q on the Euclidean space E^3

We have then::

    sage: spherical(q)
    (4, 1/3*pi, pi)
    sage: cartesian(q)
    (-2*sqrt(3), 0, 2)
    sage: cylindrical(q)
    (2*sqrt(3), pi, 2)


How to express a scalar field in various coordinate systems
-----------------------------------------------------------

Let us define a scalar field on `\mathbb{E}^3` from its expression in
Cartesian coordinates::

    sage: f = E.scalar_field(x^2+y^2 - z^2, name='f')

Note that since the Cartesian coordinates are the default ones, we have not
specified them in the above definition. Thanks to the known coordinate
transformations, the expression of `f` in terms of other coordinates is
automatically computed::

    sage: f.display()
    f: E^3 → ℝ
       (x, y, z) ↦ x^2 + y^2 - z^2
       (r, th, ph) ↦ -2*r^2*cos(th)^2 + r^2
       (rh, ph, z) ↦ rh^2 - z^2

We can limit the output to a single coordinate system::

    sage: f.display(cartesian)
    f: E^3 → ℝ
       (x, y, z) ↦ x^2 + y^2 - z^2
    sage: f.display(cylindrical)
    f: E^3 → ℝ
       (rh, ph, z) ↦ rh^2 - z^2

The coordinate expression in a given coordinate system is obtained via the
method :meth:`~sage.manifolds.scalarfield.ScalarField.expr`::

    sage: f.expr()  # expression in the default chart (Cartesian coordinates)
    x^2 + y^2 - z^2
    sage: f.expr(spherical)
    -2*r^2*cos(th)^2 + r^2
    sage: f.expr(cylindrical)
    rh^2 - z^2

The values of `f` at points `p` and `q` are::

    sage: f(p)
    2
    sage: f(q)
    8

Of course, we may define a scalar field from its coordinate expression in a
chart that is not the default one::

    sage: g = E.scalar_field(r^2, chart=spherical, name='g')

Instead of using the keyword argument ``chart``, one can pass a dictionary as
the first argument, with the chart as key::

    sage: g = E.scalar_field({spherical: r^2}, name='g')

The computation of the expressions of `g` in the other coordinate systems is
triggered by the method ``display()``::

    sage: g.display()
    g: E^3 → ℝ
       (x, y, z) ↦ x^2 + y^2 + z^2
       (r, th, ph) ↦ r^2
       (rh, ph, z) ↦ rh^2 + z^2


How to express a vector field in various frames
-----------------------------------------------

Let us introduce a vector field on `\mathbb{E}^3` by its components in the
Cartesian frame. Since the latter is the default vector frame on
`\mathbb{E}^3`, it suffices to write::

    sage: v = E.vector_field(-y, x, z^2, name='v')
    sage: v.display()
    v = -y e_x + x e_y + z^2 e_z

Equivalently, a vector field can be defined directly from its expansion on the
Cartesian frame::

    sage: ex, ey, ez = cartesian_frame[:]
    sage: v = -y*ex + x*ey + z^2*ez
    sage: v.display()
    -y e_x + x e_y + z^2 e_z

Let us provide ``v`` with some name, as above::

    sage: v.set_name('v')
    sage: v.display()
    v = -y e_x + x e_y + z^2 e_z

The components of `v` are returned by the square bracket operator::

    sage: v[1]
    -y
    sage: v[:]
    [-y, x, z^2]

The computation of the expression of `v` in terms of the orthonormal
spherical frame is triggered by the method ``display()``::

    sage: v.display(spherical_frame)
    v = z^3/sqrt(x^2 + y^2 + z^2) e_r
     - sqrt(x^2 + y^2)*z^2/sqrt(x^2 + y^2 + z^2) e_th + sqrt(x^2 + y^2) e_ph

We note that the components are still expressed in the default chart
(Cartesian coordinates). To have them expressed in the spherical chart, it
suffices to pass the latter as a second argument to ``display()``::

    sage: v.display(spherical_frame, spherical)
    v = r^2*cos(th)^3 e_r - r^2*cos(th)^2*sin(th) e_th + r*sin(th) e_ph

Again, the components of `v` are obtained by means of the square bracket
operator, by specifying the vector frame as first argument and the coordinate
chart as the last one::

    sage: v[spherical_frame, 1]
    z^3/sqrt(x^2 + y^2 + z^2)
    sage: v[spherical_frame, 1, spherical]
    r^2*cos(th)^3
    sage: v[spherical_frame, :, spherical]
    [r^2*cos(th)^3, -r^2*cos(th)^2*sin(th), r*sin(th)]

Similarly, the expression of `v` in terms of the cylindrical frame is::

    sage: v.display(cylindrical_frame, cylindrical)
    v = rh e_ph + z^2 e_z
    sage: v[cylindrical_frame, :, cylindrical]
    [0, rh, z^2]

The value of the vector field `v` at point `p` is::

    sage: vp = v.at(p)
    sage: vp
    Vector v at Point p on the Euclidean space E^3
    sage: vp.display()
    v = -e_x - e_y
    sage: vp.display(spherical_frame.at(p))
    v = sqrt(2) e_ph
    sage: vp.display(cylindrical_frame.at(p))
    v = sqrt(2) e_ph

The value of the vector field `v` at point `q` is::

    sage: vq = v.at(q)
    sage: vq
    Vector v at Point q on the Euclidean space E^3
    sage: vq.display()
    v = -2*sqrt(3) e_y + 4 e_z
    sage: vq.display(spherical_frame.at(q))
    v = 2 e_r - 2*sqrt(3) e_th + 2*sqrt(3) e_ph
    sage: vq.display(cylindrical_frame.at(q))
    v = 2*sqrt(3) e_ph + 4 e_z


How to change the default coordinates and default vector frame
--------------------------------------------------------------

At any time, one may change the default coordinates by the method
:meth:`~sage.manifolds.manifold.TopologicalManifold.set_default_chart`::

    sage: E.set_default_chart(spherical)

Then::

    sage: f.expr()
    -2*r^2*cos(th)^2 + r^2
    sage: v.display()
    v = -r*sin(ph)*sin(th) e_x + r*cos(ph)*sin(th) e_y + r^2*cos(th)^2 e_z

Note that the default vector frame is still the Cartesian one; to change to
the orthonormal spherical frame, use
:meth:`~sage.manifolds.differentiable.manifold.DifferentiableManifold.set_default_frame`::

    sage: E.set_default_frame(spherical_frame)

Then::

    sage: v.display()
    v = r^2*cos(th)^3 e_r - r^2*cos(th)^2*sin(th) e_th + r*sin(th) e_ph
    sage: v.display(cartesian_frame, cartesian)
    v = -y e_x + x e_y + z^2 e_z
